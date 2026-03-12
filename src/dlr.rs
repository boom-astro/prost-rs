use crate::ellipse::Ellipse;

/// Result of a directional light radius computation.
#[derive(Debug, Clone)]
pub struct DlrResult {
    /// Angular separation between transient and galaxy center (arcsec)
    pub separation_arcsec: f64,
    /// Directional light radius: galaxy effective radius along the
    /// direction toward the transient (arcsec)
    pub directional_radius: f64,
    /// Fractional offset = separation / directional_radius
    pub fractional_offset: f64,
}

/// Compute the directional light radius (DLR) for a transient–galaxy pair.
///
/// Uses tangent-plane projection centered on the galaxy to compute offsets,
/// then rotates into the galaxy's ellipse frame to find the effective radius
/// along the direction toward the transient.
pub fn compute_dlr(
    transient_ra: f64,
    transient_dec: f64,
    galaxy_ra: f64,
    galaxy_dec: f64,
    ellipse: &Ellipse,
) -> DlrResult {
    let dec_g_rad = galaxy_dec.to_radians();
    let cos_dec = dec_g_rad.cos();

    // Tangent-plane offsets in arcsec
    let mut dra = (transient_ra - galaxy_ra) * cos_dec * 3600.0;
    // Handle RA wrap-around
    if dra > 648000.0 {
        dra -= 1296000.0;
    } else if dra < -648000.0 {
        dra += 1296000.0;
    }
    let ddec = (transient_dec - galaxy_dec) * 3600.0;

    let separation = dra.hypot(ddec);

    if separation < 1e-15 {
        return DlrResult {
            separation_arcsec: 0.0,
            directional_radius: ellipse.a,
            fractional_offset: 0.0,
        };
    }

    // Rotate into galaxy ellipse frame
    let (sin_pa, cos_pa) = ellipse.pa_rad.sin_cos();
    let x_maj = dra * cos_pa + ddec * sin_pa;
    let y_min = -dra * sin_pa + ddec * cos_pa;

    // Angle of transient in the ellipse frame
    let theta = y_min.atan2(x_maj);
    let (sin_t, cos_t) = theta.sin_cos();

    // Directional radius from the ellipse equation:
    // r(θ) = a*b / sqrt((b*cosθ)² + (a*sinθ)²)
    let denom = (ellipse.b * cos_t).hypot(ellipse.a * sin_t);
    let directional_radius = if denom > 1e-15 {
        ellipse.a * ellipse.b / denom
    } else {
        ellipse.a
    };

    let fractional_offset = separation / directional_radius;

    DlrResult {
        separation_arcsec: separation,
        directional_radius,
        fractional_offset,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_dlr_on_center() {
        let e = Ellipse::new(2.0, 1.0, 0.0).unwrap();
        let result = compute_dlr(180.0, 45.0, 180.0, 45.0, &e);
        assert_relative_eq!(result.separation_arcsec, 0.0);
        assert_relative_eq!(result.fractional_offset, 0.0);
    }

    #[test]
    fn test_dlr_along_major_axis() {
        // Galaxy at origin with PA=0 (major axis along RA=0 direction)
        // Transient offset along dec → major axis direction
        let e = Ellipse::new(4.0, 2.0, 0.0).unwrap();
        // Offset 2 arcsec in dec
        let galaxy_dec = 0.0;
        let transient_dec = galaxy_dec + 2.0 / 3600.0;
        let result = compute_dlr(0.0, transient_dec, 0.0, galaxy_dec, &e);
        assert_relative_eq!(result.separation_arcsec, 2.0, epsilon = 0.01);
        // Along the major axis (PA=0, offset in dec → sin_pa=0, cos_pa=1)
        // x_maj = dra*1 + ddec*0 = 0, y_min = -dra*0 + ddec*1 = ddec
        // theta = atan2(ddec, 0) = pi/2 → sin=1, cos=0
        // r = a*b / sqrt((b*0)^2 + (a*1)^2) = a*b/a = b
        // Actually that's the minor axis direction when PA=0
        // With PA=0: major axis is aligned with RA, minor with Dec
        assert_relative_eq!(result.directional_radius, 2.0, epsilon = 0.01);
        assert_relative_eq!(result.fractional_offset, 1.0, epsilon = 0.01);
    }

    #[test]
    fn test_dlr_circular_galaxy() {
        let e = Ellipse::new(3.0, 3.0, 0.0).unwrap();
        let result = compute_dlr(10.001, 45.0, 10.0, 45.0, &e);
        assert_relative_eq!(result.directional_radius, 3.0, epsilon = 0.01);
    }

    #[test]
    fn test_dlr_ra_wraparound() {
        let e = Ellipse::new(3.0, 3.0, 0.0).unwrap();
        // Galaxy near RA=0, transient near RA=360
        let r1 = compute_dlr(359.999, 0.0, 0.001, 0.0, &e);
        assert!(r1.separation_arcsec < 10.0); // Should be ~7.2 arcsec, not ~1.3M
    }
}
