use serde::{Deserialize, Serialize};

use crate::dlr::{compute_dlr, DlrResult};
use crate::ellipse::Ellipse;
use crate::errors::ProstError;
use crate::likelihood::{absmag_likelihood, offset_likelihood, redshift_likelihood};
use crate::prior;
use crate::types::{GalaxyCandidate, HostCandidate, Transient};

/// Configuration for the host association algorithm.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssociationConfig {
    /// Maximum fractional offset to consider (default: 10.0)
    pub max_fractional_offset: f64,
    /// Minimum semi-minor axis in arcsec (floor for tiny galaxies)
    pub min_b_arcsec: f64,
    /// Maximum number of candidates to return
    pub max_candidates: usize,
    /// Whether to use redshift information in scoring
    pub use_redshift: bool,
    /// Whether to use absolute magnitude in scoring
    pub use_absmag: bool,
}

impl Default for AssociationConfig {
    fn default() -> Self {
        Self {
            max_fractional_offset: 10.0,
            min_b_arcsec: 0.05,
            max_candidates: 10,
            use_redshift: true,
            use_absmag: false, // Not yet implemented
        }
    }
}

/// Full result of a host galaxy association.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssociationResult {
    /// Ranked list of host candidates (best first)
    pub candidates: Vec<HostCandidate>,
    /// Probability that none of the candidates is the true host
    pub p_none: f64,
    /// Number of input galaxies considered
    pub n_considered: usize,
}

impl AssociationResult {
    /// Returns the best host candidate, if any.
    pub fn best_host(&self) -> Option<&HostCandidate> {
        self.candidates.first()
    }
}

/// Perform probabilistic host galaxy association.
///
/// Given a transient and a list of galaxy candidates, compute DLR-based
/// fractional offsets and Bayesian posterior probabilities for each candidate.
pub fn associate_host(
    transient: &Transient,
    candidates: &[GalaxyCandidate],
    config: &AssociationConfig,
) -> Result<AssociationResult, ProstError> {
    if candidates.is_empty() {
        return Err(ProstError::NoCandidates);
    }

    // Build ellipses and compute DLR for each candidate
    let mut scored: Vec<(usize, DlrResult, Ellipse)> = Vec::new();

    for (i, galaxy) in candidates.iter().enumerate() {
        let ellipse = match Ellipse::from_candidate(galaxy, config.min_b_arcsec) {
            Ok(e) => e,
            Err(_) => continue, // Skip galaxies with invalid shapes
        };

        let dlr = compute_dlr(
            transient.ra,
            transient.dec,
            galaxy.ra,
            galaxy.dec,
            &ellipse,
        );

        if dlr.fractional_offset <= config.max_fractional_offset {
            scored.push((i, dlr, ellipse));
        }
    }

    if scored.is_empty() {
        return Ok(AssociationResult {
            candidates: Vec::new(),
            p_none: 1.0,
            n_considered: candidates.len(),
        });
    }

    // Sort by fractional offset
    scored.sort_by(|a, b| a.1.fractional_offset.total_cmp(&b.1.fractional_offset));

    // Compute unnormalized posteriors
    let mut unnormalized: Vec<(usize, f64, f64, f64, f64, DlrResult)> = Vec::new();

    for (i, dlr, _ellipse) in &scored {
        let galaxy = &candidates[*i];

        // Offset component
        let l_offset = offset_likelihood(dlr.fractional_offset);
        let p_offset = prior::offset_prior(dlr.fractional_offset, config.max_fractional_offset);
        let post_offset = l_offset * p_offset;

        // Redshift component
        let post_redshift = if config.use_redshift {
            redshift_likelihood(
                galaxy.redshift,
                galaxy.redshift_err,
                transient.redshift,
                transient.redshift_err,
            )
        } else {
            1.0
        };

        // Absolute magnitude component
        let post_absmag = if config.use_absmag {
            absmag_likelihood(galaxy.mag, galaxy.mag_err, galaxy.redshift)
        } else {
            1.0
        };

        let posterior = post_offset * post_redshift * post_absmag;
        unnormalized.push((*i, posterior, post_offset, post_redshift, post_absmag, dlr.clone()));
    }

    // Null hypothesis probabilities
    let p_out = prior::p_outside(scored.len());
    let p_unobs = prior::p_unobserved();
    let p_hostless = prior::p_hostless();
    let p_null = p_out + p_unobs + p_hostless;

    // Normalize
    let total: f64 = unnormalized.iter().map(|x| x.1).sum::<f64>() + p_null;

    let p_none = if total > 0.0 {
        (p_null / total).clamp(0.0, 1.0)
    } else {
        1.0
    };

    // Build ranked candidates
    // Re-sort by posterior (descending)
    let mut results: Vec<(usize, f64, f64, f64, f64, DlrResult)> = unnormalized;
    results.sort_by(|a, b| b.1.total_cmp(&a.1));

    let mut host_candidates: Vec<HostCandidate> = Vec::new();
    for (rank, (idx, post, post_off, post_z, post_mag, dlr)) in
        results.iter().enumerate().take(config.max_candidates)
    {
        let galaxy = candidates[*idx].clone();
        let normalized_posterior = if total > 0.0 { post / total } else { 0.0 };

        host_candidates.push(HostCandidate {
            galaxy,
            separation_arcsec: dlr.separation_arcsec,
            dlr: dlr.directional_radius,
            fractional_offset: dlr.fractional_offset,
            dlr_rank: (rank + 1) as u32,
            posterior: normalized_posterior,
            posterior_offset: *post_off,
            posterior_redshift: *post_z,
            posterior_absmag: *post_mag,
        });
    }

    Ok(AssociationResult {
        candidates: host_candidates,
        p_none,
        n_considered: candidates.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Transient;
    use approx::assert_relative_eq;

    fn make_galaxy(ra: f64, dec: f64, a: f64, b: f64, pa: f64) -> GalaxyCandidate {
        GalaxyCandidate {
            ra,
            dec,
            a_arcsec: a,
            b_arcsec: b,
            pa_deg: pa,
            redshift: None,
            redshift_err: None,
            mag: None,
            mag_err: None,
            objtype: None,
            objname: None,
            catalog: None,
            shape_from_image: false,
        }
    }

    #[test]
    fn test_associate_single_nearby() {
        let transient = Transient::new(180.0, 45.0);
        let galaxy = make_galaxy(180.0, 45.0 + 1.0 / 3600.0, 5.0, 3.0, 0.0);

        let config = AssociationConfig::default();
        let result = associate_host(&transient, &[galaxy], &config).unwrap();

        assert_eq!(result.candidates.len(), 1);
        assert!(result.candidates[0].posterior > 0.5);
        assert!(result.p_none < 0.5);
    }

    #[test]
    fn test_associate_empty() {
        let transient = Transient::new(180.0, 45.0);
        let config = AssociationConfig::default();
        assert!(associate_host(&transient, &[], &config).is_err());
    }

    #[test]
    fn test_associate_ranking() {
        let transient = Transient::new(180.0, 45.0);
        // Close galaxy
        let g1 = make_galaxy(180.0, 45.0 + 0.5 / 3600.0, 5.0, 3.0, 0.0);
        // Far galaxy
        let g2 = make_galaxy(180.0, 45.0 + 10.0 / 3600.0, 5.0, 3.0, 0.0);

        let config = AssociationConfig::default();
        let result = associate_host(&transient, &[g1, g2], &config).unwrap();

        assert_eq!(result.candidates.len(), 2);
        assert!(result.candidates[0].posterior > result.candidates[1].posterior);
        assert_eq!(result.candidates[0].dlr_rank, 1);
        assert_eq!(result.candidates[1].dlr_rank, 2);
    }

    #[test]
    fn test_associate_with_redshift() {
        let transient = Transient::new(180.0, 45.0).with_redshift(0.05, 0.001);

        // Galaxy at matching redshift
        let mut g1 = make_galaxy(180.0, 45.0 + 2.0 / 3600.0, 5.0, 3.0, 0.0);
        g1.redshift = Some(0.05);
        g1.redshift_err = Some(0.001);

        // Same position, wrong redshift
        let mut g2 = make_galaxy(180.0, 45.0 + 2.0 / 3600.0, 5.0, 3.0, 0.0);
        g2.redshift = Some(0.5);
        g2.redshift_err = Some(0.001);

        let config = AssociationConfig {
            use_redshift: true,
            ..Default::default()
        };
        let result = associate_host(&transient, &[g1, g2], &config).unwrap();

        // The redshift-matching galaxy should be ranked higher
        assert!(result.candidates[0].posterior > result.candidates[1].posterior);
        assert_relative_eq!(
            result.candidates[0].galaxy.redshift.unwrap(),
            0.05,
            epsilon = 0.001
        );
    }

    #[test]
    fn test_posteriors_sum_to_one() {
        let transient = Transient::new(180.0, 45.0);
        let galaxies: Vec<GalaxyCandidate> = (1..=5)
            .map(|i| {
                make_galaxy(
                    180.0,
                    45.0 + (i as f64) * 1.0 / 3600.0,
                    4.0,
                    2.0,
                    0.0,
                )
            })
            .collect();

        let config = AssociationConfig::default();
        let result = associate_host(&transient, &galaxies, &config).unwrap();

        let sum: f64 = result.candidates.iter().map(|c| c.posterior).sum::<f64>() + result.p_none;
        assert_relative_eq!(sum, 1.0, epsilon = 1e-6);
    }
}
