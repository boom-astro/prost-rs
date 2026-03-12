//! Likelihood functions for host galaxy association.
//!
//! The primary likelihood is based on the fractional offset (separation / DLR)
//! following a Gamma distribution, as described in Prost (Gagliano et al.).

/// Γ(0.75) — the gamma function evaluated at 0.75.
const GAMMA_0_75: f64 = 1.2254167024651776;

/// Compute the offset likelihood using a Gamma(a=0.75) distribution.
///
/// PDF: f(x; a=0.75) = x^(a-1) * exp(-x) / Γ(a)
///                    = x^(-0.25) * exp(-x) / Γ(0.75)
///
/// This models the distribution of fractional offsets (separation/DLR)
/// for true host galaxies.
pub fn offset_likelihood(fractional_offset: f64) -> f64 {
    if fractional_offset < 0.0 {
        return 0.0;
    }
    if fractional_offset == 0.0 {
        // Gamma(a<1) PDF diverges at x=0; use a large but finite value
        // consistent with BOOM's host.rs behavior
        return 1e6;
    }
    fractional_offset.powf(-0.25) * (-fractional_offset).exp() / GAMMA_0_75
}

/// Compute the offset likelihood with Monte Carlo position uncertainty.
///
/// Draws `n_samples` positions from a 2D Gaussian centered on the transient
/// position with the given RA/Dec uncertainties, computes the fractional
/// offset for each, and returns the mean likelihood.
pub fn offset_likelihood_mc(
    fractional_offsets: &[f64],
) -> f64 {
    if fractional_offsets.is_empty() {
        return 0.0;
    }
    let sum: f64 = fractional_offsets.iter().map(|&fo| offset_likelihood(fo)).sum();
    sum / fractional_offsets.len() as f64
}

/// Compute the redshift likelihood.
///
/// If both transient and galaxy have redshifts, use a Gaussian likelihood
/// centered on the transient's redshift. If only the galaxy has a redshift,
/// use a broad prior. If neither has a redshift, return 1.0 (uninformative).
pub fn redshift_likelihood(
    galaxy_z: Option<f64>,
    galaxy_z_err: Option<f64>,
    transient_z: Option<f64>,
    transient_z_err: Option<f64>,
) -> f64 {
    match (galaxy_z, transient_z) {
        (Some(gz), Some(tz)) => {
            let gz_err = galaxy_z_err.unwrap_or(0.01);
            let tz_err = transient_z_err.unwrap_or(0.01);
            let sigma2 = gz_err * gz_err + tz_err * tz_err;
            let dz = gz - tz;
            (-0.5 * dz * dz / sigma2).exp()
        }
        (Some(_gz), None) => {
            // Galaxy has redshift but transient doesn't — mildly informative
            // Use a flat prior (all redshifts equally likely)
            1.0
        }
        _ => 1.0,
    }
}

/// Compute the absolute magnitude likelihood.
///
/// Given apparent magnitude and redshift, compute absolute magnitude
/// and evaluate against a Schechter-like luminosity function.
/// For now, returns 1.0 (uninformative) — to be implemented with
/// proper luminosity function and k-corrections.
pub fn absmag_likelihood(
    _mag: Option<f64>,
    _mag_err: Option<f64>,
    _redshift: Option<f64>,
) -> f64 {
    // TODO: Implement Schechter luminosity function
    // M = m - 5*log10(d_L/10pc) - K(z)
    // L(M) ∝ 10^(0.4*(M*-M)*(α+1)) * exp(-10^(0.4*(M*-M)))
    1.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_offset_likelihood_zero() {
        // At FO=0, Gamma(a<1) diverges; we return a large finite value
        assert_eq!(offset_likelihood(0.0), 1e6);
    }

    #[test]
    fn test_offset_likelihood_peak() {
        // Gamma(0.75) peaks near x = a-1 = -0.25 → but since a<1,
        // PDF diverges as x→0+. So it should be monotonically decreasing
        // after a small x.
        let l1 = offset_likelihood(0.1);
        let l2 = offset_likelihood(1.0);
        let l3 = offset_likelihood(5.0);
        assert!(l1 > l2);
        assert!(l2 > l3);
    }

    #[test]
    fn test_offset_likelihood_at_one() {
        // f(1) = 1^(-0.25) * exp(-1) / Γ(0.75) = exp(-1) / Γ(0.75)
        let expected = (-1.0_f64).exp() / GAMMA_0_75;
        assert_relative_eq!(offset_likelihood(1.0), expected, epsilon = 1e-10);
    }

    #[test]
    fn test_redshift_likelihood_matching() {
        let l = redshift_likelihood(Some(0.05), Some(0.001), Some(0.05), Some(0.001));
        assert_relative_eq!(l, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_redshift_likelihood_discrepant() {
        let l = redshift_likelihood(Some(0.05), Some(0.001), Some(0.5), Some(0.001));
        assert!(l < 1e-10); // Very different redshifts → near-zero likelihood
    }

    #[test]
    fn test_redshift_likelihood_no_info() {
        assert_relative_eq!(
            redshift_likelihood(None, None, None, None),
            1.0
        );
    }

    #[test]
    fn test_mc_likelihood() {
        let offsets = vec![0.5, 1.0, 1.5, 2.0];
        let mc = offset_likelihood_mc(&offsets);
        let expected: f64 = offsets.iter().map(|&fo| offset_likelihood(fo)).sum::<f64>() / 4.0;
        assert_relative_eq!(mc, expected, epsilon = 1e-10);
    }
}
