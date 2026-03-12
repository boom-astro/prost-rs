//! Integration tests using known transient–host associations.
//!
//! These tests verify that prost-rs correctly identifies host galaxies for
//! well-studied transients where the host association is unambiguous.
//!
//! Tests marked `#[ignore]` require network access (BOOM API, Legacy Survey).
//! Run them with: cargo test -- --ignored

use approx::assert_relative_eq;
use prost_rs::associate::{associate_host, AssociationConfig, AssociationResult};
use prost_rs::dlr::compute_dlr;
use prost_rs::ellipse::Ellipse;
use prost_rs::likelihood::{offset_likelihood, redshift_likelihood};
use prost_rs::types::{GalaxyCandidate, Transient};

// ---------------------------------------------------------------------------
// Well-known transient–host pairs (literature values)
// ---------------------------------------------------------------------------

/// SN 2011fe in M101 (NGC 5457)
/// One of the closest Type Ia SNe, unambiguous host.
/// SN: RA=210.774, Dec=+54.274
/// M101: RA=210.802, Dec=+54.349, a~700", b~700" (nearly face-on)
fn sn2011fe() -> (Transient, GalaxyCandidate) {
    let transient = Transient::new(210.7742, 54.2739)
        .with_redshift(0.000804, 0.000010);
    let galaxy = GalaxyCandidate {
        ra: 210.8024,
        dec: 54.3488,
        a_arcsec: 840.0, // ~14 arcmin semi-major
        b_arcsec: 840.0, // nearly circular
        pa_deg: 0.0,
        redshift: Some(0.000804),
        redshift_err: Some(0.000010),
        mag: Some(7.86),
        mag_err: Some(0.02),
        objtype: Some("G".to_string()),
        objname: Some("NGC 5457".to_string()),
        catalog: Some("NED".to_string()),
        shape_from_image: false,
    };
    (transient, galaxy)
}

/// SN 1987A in the LMC
/// SN: RA=83.8667, Dec=-69.2697
/// LMC center: RA=80.894, Dec=-69.756
fn sn1987a() -> (Transient, GalaxyCandidate) {
    let transient = Transient::new(83.8667, -69.2697);
    let galaxy = GalaxyCandidate {
        ra: 80.894,
        dec: -69.756,
        a_arcsec: 18000.0, // ~5 degrees semi-major
        b_arcsec: 9000.0,
        pa_deg: 170.0,
        redshift: Some(0.000927),
        redshift_err: Some(0.000030),
        mag: Some(0.9),
        mag_err: Some(0.1),
        objtype: Some("G".to_string()),
        objname: Some("LMC".to_string()),
        catalog: Some("NED".to_string()),
        shape_from_image: false,
    };
    (transient, galaxy)
}

/// AT2017gfo (GW170817 kilonova) in NGC 4993
/// SN: RA=197.4504, Dec=-23.3815
/// NGC 4993: RA=197.4487, Dec=-23.3839, r_eff~20"
fn at2017gfo() -> (Transient, GalaxyCandidate) {
    let transient = Transient::new(197.4504, -23.3815)
        .with_redshift(0.009787, 0.000150);
    let galaxy = GalaxyCandidate {
        ra: 197.4487,
        dec: -23.3839,
        a_arcsec: 22.0,
        b_arcsec: 16.0,
        pa_deg: 75.0,
        redshift: Some(0.009787),
        redshift_err: Some(0.000150),
        mag: Some(13.05),
        mag_err: Some(0.05),
        objtype: Some("E/S0".to_string()),
        objname: Some("NGC 4993".to_string()),
        catalog: Some("NED".to_string()),
        shape_from_image: false,
    };
    (transient, galaxy)
}

/// SN 2014J in M82
/// SN: RA=148.9256, Dec=+69.6739
/// M82: RA=148.9685, Dec=+69.6797, a~335", b~155", PA=65°
fn sn2014j() -> (Transient, GalaxyCandidate) {
    let transient = Transient::new(148.9256, 69.6739)
        .with_redshift(0.000677, 0.000010);
    let galaxy = GalaxyCandidate {
        ra: 148.9685,
        dec: 69.6797,
        a_arcsec: 335.0,
        b_arcsec: 155.0,
        pa_deg: 65.0,
        redshift: Some(0.000677),
        redshift_err: Some(0.000040),
        mag: Some(8.41),
        mag_err: Some(0.05),
        objtype: Some("Irr".to_string()),
        objname: Some("M82".to_string()),
        catalog: Some("NED".to_string()),
        shape_from_image: false,
    };
    (transient, galaxy)
}

// ---------------------------------------------------------------------------
// DLR tests with known geometry
// ---------------------------------------------------------------------------

#[test]
fn test_dlr_sn2011fe_inside_m101() {
    let (transient, galaxy) = sn2011fe();
    let ellipse = Ellipse::new(galaxy.a_arcsec, galaxy.b_arcsec, galaxy.pa_deg).unwrap();
    let result = compute_dlr(transient.ra, transient.dec, galaxy.ra, galaxy.dec, &ellipse);
    // SN 2011fe is well within M101's disk (~280" from center, r_eff~840")
    assert!(
        result.fractional_offset < 1.0,
        "SN 2011fe should be within M101's effective radius, got FO={}",
        result.fractional_offset
    );
    assert!(result.separation_arcsec > 200.0 && result.separation_arcsec < 400.0);
}

#[test]
fn test_dlr_at2017gfo_near_ngc4993() {
    let (transient, galaxy) = at2017gfo();
    let ellipse = Ellipse::new(galaxy.a_arcsec, galaxy.b_arcsec, galaxy.pa_deg).unwrap();
    let result = compute_dlr(transient.ra, transient.dec, galaxy.ra, galaxy.dec, &ellipse);
    // AT2017gfo is ~10" from NGC 4993 center, r_eff~22"
    assert!(
        result.fractional_offset < 1.0,
        "AT2017gfo should be within NGC 4993, got FO={}",
        result.fractional_offset
    );
    assert!(result.separation_arcsec > 5.0 && result.separation_arcsec < 15.0);
}

#[test]
fn test_dlr_sn2014j_inside_m82() {
    let (transient, galaxy) = sn2014j();
    let ellipse = Ellipse::new(galaxy.a_arcsec, galaxy.b_arcsec, galaxy.pa_deg).unwrap();
    let result = compute_dlr(transient.ra, transient.dec, galaxy.ra, galaxy.dec, &ellipse);
    // SN 2014J is ~54" from M82 center, well within the disk
    assert!(
        result.fractional_offset < 1.0,
        "SN 2014J should be within M82, got FO={}",
        result.fractional_offset
    );
}

#[test]
fn test_dlr_sn1987a_inside_lmc() {
    let (transient, galaxy) = sn1987a();
    let ellipse = Ellipse::new(galaxy.a_arcsec, galaxy.b_arcsec, galaxy.pa_deg).unwrap();
    let result = compute_dlr(transient.ra, transient.dec, galaxy.ra, galaxy.dec, &ellipse);
    // SN 1987A is inside the LMC
    assert!(
        result.fractional_offset < 1.0,
        "SN 1987A should be within the LMC, got FO={}",
        result.fractional_offset
    );
}

// ---------------------------------------------------------------------------
// Full association tests
// ---------------------------------------------------------------------------

#[test]
fn test_associate_at2017gfo_correct_host() {
    let (transient, host) = at2017gfo();

    // Add some decoy galaxies at various distances
    let decoy1 = GalaxyCandidate {
        ra: 197.50,
        dec: -23.40,
        a_arcsec: 5.0,
        b_arcsec: 3.0,
        pa_deg: 0.0,
        redshift: Some(0.15),
        redshift_err: Some(0.01),
        mag: Some(19.0),
        mag_err: Some(0.1),
        objtype: Some("G".to_string()),
        objname: Some("decoy1".to_string()),
        catalog: None,
        shape_from_image: false,
    };
    let decoy2 = GalaxyCandidate {
        ra: 197.46,
        dec: -23.37,
        a_arcsec: 3.0,
        b_arcsec: 2.0,
        pa_deg: 30.0,
        redshift: Some(0.08),
        redshift_err: Some(0.02),
        mag: Some(20.5),
        mag_err: Some(0.2),
        objtype: Some("G".to_string()),
        objname: Some("decoy2".to_string()),
        catalog: None,
        shape_from_image: false,
    };

    let config = AssociationConfig {
        use_redshift: true,
        ..Default::default()
    };
    let result = associate_host(&transient, &[host.clone(), decoy1, decoy2], &config).unwrap();

    assert!(!result.candidates.is_empty());
    let best = result.best_host().unwrap();
    assert_eq!(
        best.galaxy.objname.as_deref(),
        Some("NGC 4993"),
        "AT2017gfo should be associated with NGC 4993"
    );
    assert!(
        best.posterior > 0.5,
        "NGC 4993 should have high posterior, got {}",
        best.posterior
    );
}

#[test]
fn test_associate_sn2014j_correct_host() {
    let (transient, host) = sn2014j();

    let decoy = GalaxyCandidate {
        ra: 149.0,
        dec: 69.7,
        a_arcsec: 10.0,
        b_arcsec: 8.0,
        pa_deg: 0.0,
        redshift: None,
        redshift_err: None,
        mag: Some(17.0),
        mag_err: Some(0.1),
        objtype: Some("G".to_string()),
        objname: Some("field galaxy".to_string()),
        catalog: None,
        shape_from_image: false,
    };

    let config = AssociationConfig::default();
    let result = associate_host(&transient, &[host.clone(), decoy], &config).unwrap();

    let best = result.best_host().unwrap();
    assert_eq!(best.galaxy.objname.as_deref(), Some("M82"));
    assert!(best.posterior > 0.5);
}

#[test]
fn test_associate_redshift_breaks_degeneracy() {
    // Two galaxies at similar angular separation, but only one has matching redshift
    let transient = Transient::new(180.0, 45.0).with_redshift(0.03, 0.001);

    let matching_z = GalaxyCandidate {
        ra: 180.0 + 3.0 / 3600.0,
        dec: 45.0,
        a_arcsec: 8.0,
        b_arcsec: 6.0,
        pa_deg: 0.0,
        redshift: Some(0.03),
        redshift_err: Some(0.002),
        mag: Some(16.0),
        mag_err: Some(0.05),
        objtype: Some("G".to_string()),
        objname: Some("match_z".to_string()),
        catalog: None,
        shape_from_image: false,
    };

    let wrong_z = GalaxyCandidate {
        ra: 180.0 + 2.5 / 3600.0,
        dec: 45.0,
        a_arcsec: 8.0,
        b_arcsec: 6.0,
        pa_deg: 0.0,
        redshift: Some(0.5),
        redshift_err: Some(0.002),
        mag: Some(16.0),
        mag_err: Some(0.05),
        objtype: Some("G".to_string()),
        objname: Some("wrong_z".to_string()),
        catalog: None,
        shape_from_image: false,
    };

    let config = AssociationConfig {
        use_redshift: true,
        ..Default::default()
    };

    let result = associate_host(&transient, &[matching_z, wrong_z], &config).unwrap();
    let best = result.best_host().unwrap();
    assert_eq!(
        best.galaxy.objname.as_deref(),
        Some("match_z"),
        "Redshift match should win even if slightly farther"
    );
}

// ---------------------------------------------------------------------------
// Edge cases and robustness
// ---------------------------------------------------------------------------

#[test]
fn test_associate_all_far_away() {
    // All galaxies beyond max_fractional_offset
    let transient = Transient::new(180.0, 45.0);
    let far = GalaxyCandidate {
        ra: 180.1,
        dec: 45.1,
        a_arcsec: 1.0,
        b_arcsec: 0.5,
        pa_deg: 0.0,
        redshift: None,
        redshift_err: None,
        mag: None,
        mag_err: None,
        objtype: None,
        objname: None,
        catalog: None,
        shape_from_image: false,
    };

    let config = AssociationConfig::default();
    let result = associate_host(&transient, &[far], &config).unwrap();
    assert!(result.candidates.is_empty());
    assert_relative_eq!(result.p_none, 1.0, epsilon = 1e-6);
}

#[test]
fn test_associate_single_galaxy_on_top() {
    // Galaxy center coincides with transient — FO ≈ 0
    let transient = Transient::new(200.0, -30.0);
    let galaxy = GalaxyCandidate {
        ra: 200.0,
        dec: -30.0,
        a_arcsec: 10.0,
        b_arcsec: 5.0,
        pa_deg: 45.0,
        redshift: None,
        redshift_err: None,
        mag: None,
        mag_err: None,
        objtype: None,
        objname: Some("on_top".to_string()),
        catalog: None,
        shape_from_image: false,
    };

    let config = AssociationConfig::default();
    let result = associate_host(&transient, &[galaxy], &config).unwrap();
    let best = result.best_host().unwrap();
    assert_relative_eq!(best.fractional_offset, 0.0, epsilon = 1e-10);
    assert!(best.posterior > 0.9, "Galaxy on top should dominate, got {}", best.posterior);
}

#[test]
fn test_associate_many_candidates_posteriors_sum() {
    let transient = Transient::new(180.0, 0.0);
    let galaxies: Vec<GalaxyCandidate> = (1..=20)
        .map(|i| GalaxyCandidate {
            ra: 180.0 + (i as f64) * 0.5 / 3600.0,
            dec: 0.0 + (i as f64) * 0.3 / 3600.0,
            a_arcsec: 5.0,
            b_arcsec: 3.0,
            pa_deg: 0.0,
            redshift: None,
            redshift_err: None,
            mag: None,
            mag_err: None,
            objtype: None,
            objname: Some(format!("galaxy_{i}")),
            catalog: None,
            shape_from_image: false,
        })
        .collect();

    let config = AssociationConfig {
        max_candidates: 20,
        ..Default::default()
    };
    let result = associate_host(&transient, &galaxies, &config).unwrap();

    let total: f64 = result.candidates.iter().map(|c| c.posterior).sum::<f64>() + result.p_none;
    assert_relative_eq!(total, 1.0, epsilon = 1e-6);
}

#[test]
fn test_associate_highly_elongated_galaxy() {
    // Edge-on galaxy: test that DLR works correctly with extreme axis ratios
    let transient = Transient::new(100.0, 20.0);
    let galaxy = GalaxyCandidate {
        ra: 100.0,
        dec: 20.0 + 5.0 / 3600.0, // 5 arcsec N
        a_arcsec: 30.0,
        b_arcsec: 2.0,
        pa_deg: 0.0, // major axis along RA
        redshift: None,
        redshift_err: None,
        mag: None,
        mag_err: None,
        objtype: None,
        objname: Some("edge_on".to_string()),
        catalog: None,
        shape_from_image: false,
    };

    let config = AssociationConfig::default();
    let result = associate_host(&transient, &[galaxy], &config).unwrap();
    // Transient is 5" north, along the minor axis (b=2")
    // So FO along minor axis should be > 1.0
    let best = result.best_host().unwrap();
    assert!(
        best.fractional_offset > 1.0,
        "5\" offset along minor axis (b=2\") should give FO > 1, got {}",
        best.fractional_offset
    );
}

#[test]
fn test_associate_transient_along_major_axis() {
    // Transient along the major axis of an elongated galaxy
    let transient = Transient::new(100.0 + 10.0 / 3600.0, 20.0);
    let galaxy = GalaxyCandidate {
        ra: 100.0,
        dec: 20.0,
        a_arcsec: 30.0,
        b_arcsec: 2.0,
        pa_deg: 0.0, // major axis along RA
        redshift: None,
        redshift_err: None,
        mag: None,
        mag_err: None,
        objtype: None,
        objname: Some("edge_on_along".to_string()),
        catalog: None,
        shape_from_image: false,
    };

    let config = AssociationConfig::default();
    let result = associate_host(&transient, &[galaxy], &config).unwrap();
    let best = result.best_host().unwrap();
    // 10" offset along major axis (a=30") → FO ~ 0.33
    assert!(
        best.fractional_offset < 1.0,
        "10\" along major axis (a=30\") should give FO < 1, got {}",
        best.fractional_offset
    );
}

// ---------------------------------------------------------------------------
// Likelihood function tests
// ---------------------------------------------------------------------------

#[test]
fn test_offset_likelihood_monotonic_decrease() {
    // Gamma(0.75) has a=0.75 < 1, so PDF is monotonically decreasing for x > 0
    let mut prev = offset_likelihood(0.01);
    for i in 1..100 {
        let x = i as f64 * 0.1;
        let l = offset_likelihood(x);
        assert!(
            l <= prev,
            "Offset likelihood should be monotonically decreasing: f({}) = {} > f({}) = {}",
            x - 0.1,
            prev,
            x,
            l
        );
        prev = l;
    }
}

#[test]
fn test_offset_likelihood_integrates_to_one() {
    // Numerical integration of Gamma(0.75) PDF should be ~1
    let n = 100_000;
    let dx = 50.0 / n as f64;
    let integral: f64 = (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) * dx;
            offset_likelihood(x) * dx
        })
        .sum();
    assert_relative_eq!(integral, 1.0, epsilon = 0.01);
}

#[test]
fn test_redshift_likelihood_symmetric() {
    let l1 = redshift_likelihood(Some(0.05), Some(0.01), Some(0.06), Some(0.01));
    let l2 = redshift_likelihood(Some(0.06), Some(0.01), Some(0.05), Some(0.01));
    assert_relative_eq!(l1, l2, epsilon = 1e-10);
}

#[test]
fn test_redshift_likelihood_wider_errors_more_forgiving() {
    // With tighter errors, a given dz should produce lower likelihood
    let l_tight = redshift_likelihood(Some(0.05), Some(0.001), Some(0.07), Some(0.001));
    let l_loose = redshift_likelihood(Some(0.05), Some(0.05), Some(0.07), Some(0.05));
    assert!(
        l_loose > l_tight,
        "Wider error bars should produce higher likelihood for same dz"
    );
}

// ---------------------------------------------------------------------------
// Ellipse / Tractor shape tests
// ---------------------------------------------------------------------------

#[test]
fn test_tractor_known_shapes() {
    // A round DEV source with shape_r=2.0, e1=0, e2=0 → circular
    let e = Ellipse::from_tractor(2.0, 0.0, 0.0, 0.05).unwrap();
    assert_relative_eq!(e.a, 2.0, epsilon = 1e-10);
    assert_relative_eq!(e.b, 2.0, epsilon = 1e-10);
    assert_relative_eq!(e.axis_ratio, 1.0, epsilon = 1e-10);

    // An elongated source with e1=0.6, e2=0 → q = 0.4/1.6 = 0.25
    let e = Ellipse::from_tractor(3.0, 0.6, 0.0, 0.05).unwrap();
    assert_relative_eq!(e.axis_ratio, 0.25, epsilon = 1e-10);
    assert_relative_eq!(e.b, 3.0 * 0.25, epsilon = 1e-10);
    assert_relative_eq!(e.pa_rad, 0.0, epsilon = 1e-10);

    // PA from e2: e1=0, e2=0.5 → PA = 0.5*atan2(0.5, 0) = pi/4
    let e = Ellipse::from_tractor(2.0, 0.0, 0.5, 0.05).unwrap();
    assert_relative_eq!(
        e.pa_rad,
        std::f64::consts::FRAC_PI_4,
        epsilon = 1e-10
    );
}

// ---------------------------------------------------------------------------
// DLR geometry tests
// ---------------------------------------------------------------------------

#[test]
fn test_dlr_diagonal_offset() {
    // Galaxy at origin, PA=45°, a=10", b=5"
    // Transient offset diagonally
    let e = Ellipse::new(10.0, 5.0, 45.0).unwrap();
    let result = compute_dlr(0.0 + 5.0 / 3600.0, 0.0 + 5.0 / 3600.0, 0.0, 0.0, &e);
    assert!(result.separation_arcsec > 6.0 && result.separation_arcsec < 8.0);
    assert!(result.fractional_offset > 0.0);
}

#[test]
fn test_dlr_high_declination() {
    // Verify cos(dec) correction at high declination
    let e = Ellipse::new(5.0, 5.0, 0.0).unwrap();

    // At dec=80°, 1° of RA = cos(80°) * 3600" ≈ 625"
    let result = compute_dlr(100.001, 80.0, 100.0, 80.0, &e);
    // 0.001° RA at dec=80 ≈ 0.001 * cos(80°) * 3600 ≈ 0.625"
    assert!(
        result.separation_arcsec < 1.0,
        "RA offset at high dec should be small, got {}",
        result.separation_arcsec
    );
}

#[test]
fn test_dlr_near_south_pole() {
    // Test near dec = -89°
    let e = Ellipse::new(10.0, 10.0, 0.0).unwrap();
    let result = compute_dlr(45.0, -89.0 + 5.0 / 3600.0, 45.0, -89.0, &e);
    assert_relative_eq!(result.separation_arcsec, 5.0, epsilon = 0.1);
}

// ---------------------------------------------------------------------------
// Serialization round-trip
// ---------------------------------------------------------------------------

#[test]
fn test_transient_serde_roundtrip() {
    let t = Transient::new(180.0, 45.0)
        .with_redshift(0.05, 0.001)
        .with_position_err(0.1, 0.1);
    let json = serde_json::to_string(&t).unwrap();
    let t2: Transient = serde_json::from_str(&json).unwrap();
    assert_relative_eq!(t2.ra, t.ra);
    assert_relative_eq!(t2.dec, t.dec);
    assert_eq!(t2.redshift, t.redshift);
}

#[test]
fn test_association_result_serde() {
    let (transient, galaxy) = at2017gfo();
    let config = AssociationConfig::default();
    let result = associate_host(&transient, &[galaxy], &config).unwrap();

    let json = serde_json::to_string(&result).unwrap();
    let result2: AssociationResult = serde_json::from_str(&json).unwrap();
    assert_eq!(result2.candidates.len(), result.candidates.len());
    assert_relative_eq!(result2.p_none, result.p_none, epsilon = 1e-10);
}
