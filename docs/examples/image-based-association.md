# Image-Based Host Association

This example demonstrates the full image-based pipeline: loading a FITS
cutout into a `Cutout`, extracting sources, fitting Sersic profiles, and
running host association — all without a pre-existing catalog.

We use **SN 2020jfo** (ZTF20abjbgjj), a Type II supernova in the barred
spiral galaxy **M61** (NGC 4303) at z=0.00522. This is a well-studied
BOOM transient with an unambiguous host.

## Step 1: Download and load a FITS cutout

The FITS download itself uses HTTP (via `reqwest` in Rust or `requests`
in Python). Once downloaded, we parse the pixel data and wrap it in a
prost `Cutout`.

=== "Rust"

    ```rust
    use prost_rs::cutout::{Cutout, ImagingSurvey};

    // SN 2020jfo coordinates (from ZTF alerts / TNS)
    let transient_ra = 185.4788;
    let transient_dec = 4.4700;
    let size_pix: u32 = 256;
    let pixel_scale = 0.262; // arcsec/pixel for Legacy Survey

    // Download Legacy Survey DR10 r-band cutout
    let url = format!(
        "https://www.legacysurvey.org/viewer/cutout.fits\
         ?ra={}&dec={}&size={}&layer=ls-dr10&pixscale={}&bands=r",
        transient_ra, transient_dec, size_pix, pixel_scale,
    );
    let resp = reqwest::blocking::get(&url)
        .expect("cutout download failed");
    let fits_bytes = resp.bytes().unwrap();

    // Parse FITS with fitsio or a simple header skip:
    // Legacy Survey cutouts have a single 2D HDU, 256x256 f32 pixels.
    // Skip the FITS header (2880-byte blocks) to reach the data.
    let data: Vec<f64> = parse_fits_image_f64(&fits_bytes, size_pix as usize);

    let cutout = Cutout::new(
        data,
        size_pix as usize,  // width
        size_pix as usize,  // height
        pixel_scale,
        transient_ra,
        transient_dec,
        ImagingSurvey::LegacySurvey,
    ).expect("invalid cutout dimensions");

    let (w_arcsec, h_arcsec) = cutout.size_arcsec();
    println!("Cutout: {}x{} px, {:.1}\"x{:.1}\"",
        cutout.width, cutout.height, w_arcsec, h_arcsec);
    ```

=== "Python"

    ```python
    import numpy as np
    from astropy.io import fits
    import requests
    from io import BytesIO
    from prost import Cutout

    transient_ra  = 185.4788
    transient_dec = 4.4700
    size_pix = 256
    pixel_scale = 0.262  # arcsec/pixel for Legacy Survey

    resp = requests.get(
        "https://www.legacysurvey.org/viewer/cutout.fits",
        params={
            "ra": transient_ra, "dec": transient_dec,
            "size": size_pix, "layer": "ls-dr10",
            "pixscale": pixel_scale, "bands": "r",
        },
    )
    resp.raise_for_status()

    hdu = fits.open(BytesIO(resp.content))
    image_data = hdu[0].data.astype(np.float64)
    if image_data.ndim == 3:
        image_data = image_data[0]

    cutout = Cutout(
        data=image_data,
        pixel_scale=pixel_scale,
        center_ra=transient_ra,
        center_dec=transient_dec,
        survey="legacy",
    )
    print(cutout)
    ```

## Step 2: Estimate background and extract sources

=== "Rust"

    ```rust
    use prost_rs::source::{estimate_background, extract_sources, ExtractionConfig};

    let config = ExtractionConfig {
        thresh_sigma: 1.5,
        min_pixels: 5,
        back_size: 64,
        ..Default::default()
    };

    let bg = estimate_background(&cutout, &config).unwrap();
    println!("Background: mean={:.4}, rms={:.4}", bg.global_mean, bg.global_rms);

    let sources = extract_sources(&cutout, &config).unwrap();
    println!("Detected {} sources", sources.len());

    for (i, s) in sources.iter().enumerate().take(5) {
        println!(
            "  [{}] x={:.1} y={:.1} flux={:.1} SNR={:.1} a={:.1}px b={:.1}px",
            i, s.x, s.y, s.flux, s.snr, s.a_pix, s.b_pix,
        );
    }
    ```

=== "Python"

    ```python
    from prost import ExtractionConfig, extract_sources, estimate_background

    config = ExtractionConfig(thresh_sigma=1.5, min_pixels=5)

    bg = estimate_background(cutout)
    print(f"Background: mean={bg.global_mean:.4f}, rms={bg.global_rms:.4f}")

    sources = extract_sources(cutout, config)
    print(f"Detected {len(sources)} sources")

    for i, s in enumerate(sources[:5]):
        print(f"  [{i}] x={s.x:.1f} y={s.y:.1f} flux={s.flux:.1f} "
              f"SNR={s.snr:.1f} a={s.a_pix:.1f}px b={s.b_pix:.1f}px")
    ```

## Step 3: Fit Sersic profiles to resolved sources

Extended sources (galaxies) can be distinguished from point sources
(stars) by their size relative to the PSF. We fit Sersic profiles to
sources with semi-major axes larger than ~2 pixels, then convert the
fitted morphology to `GalaxyCandidate`s for association.

=== "Rust"

    ```rust
    use prost_rs::morphology::{fit_sersic, FitConfig};
    use prost_rs::types::GalaxyCandidate;

    let fit_config = FitConfig {
        max_sersic_n: 8.0,
        min_r_eff_pix: 0.5,
        ..Default::default()
    };

    let mut galaxy_candidates: Vec<GalaxyCandidate> = Vec::new();

    for (idx, src) in sources.iter().enumerate() {
        // Skip likely stars (unresolved at LS resolution)
        if src.a_pix < 2.0 {
            continue;
        }

        let fit = match fit_sersic(&cutout, src, &fit_config) {
            Ok(f) if f.converged => f,
            _ => continue,
        };

        println!(
            "  Sersic fit: n={:.2}, r_eff={:.2}\", q={:.2}, PA={:.1}°, chi2r={:.2}",
            fit.n, fit.r_eff, fit.axis_ratio, fit.pa_deg, fit.chi2_reduced,
        );

        // Convert fitted morphology to a GalaxyCandidate.
        // Use the half-light ellipse for DLR: a ~ r_eff/q, b ~ r_eff.
        let mut a_arcsec = fit.r_eff / fit.axis_ratio;
        let mut b_arcsec = fit.r_eff;
        if a_arcsec < b_arcsec {
            std::mem::swap(&mut a_arcsec, &mut b_arcsec);
        }

        galaxy_candidates.push(GalaxyCandidate {
            ra: fit.ra,
            dec: fit.dec,
            a_arcsec,
            b_arcsec,
            pa_deg: fit.pa_deg,
            redshift: None,
            redshift_err: None,
            mag: None,
            mag_err: None,
            objtype: None,
            objname: Some(format!("src_{idx}")),
            catalog: None,
            shape_from_image: true,
        });
    }

    println!("{} galaxy candidates from image", galaxy_candidates.len());
    ```

=== "Python"

    ```python
    from prost import FitConfig, fit_sersic, GalaxyCandidate

    fit_cfg = FitConfig(max_sersic_n=8.0, min_r_eff_pix=0.5)

    galaxy_candidates = []
    for i, src in enumerate(sources):
        if src.a_pix < 2.0:
            continue

        try:
            fit = fit_sersic(cutout, src, fit_cfg)
        except RuntimeError:
            continue

        if not fit.converged:
            continue

        print(f"  Sersic fit: n={fit.n:.2f}, r_eff={fit.r_eff:.2f}\", "
              f"q={fit.axis_ratio:.2f}, PA={fit.pa_deg:.1f}°, "
              f"chi2r={fit.chi2_reduced:.2f}")

        a_arcsec = fit.r_eff / fit.axis_ratio
        b_arcsec = fit.r_eff
        if a_arcsec < b_arcsec:
            a_arcsec, b_arcsec = b_arcsec, a_arcsec

        galaxy_candidates.append(GalaxyCandidate(
            ra=fit.ra, dec=fit.dec,
            a_arcsec=a_arcsec, b_arcsec=b_arcsec,
            pa_deg=fit.pa_deg,
            objname=f"src_{i}",
            shape_from_image=True,
        ))

    print(f"{len(galaxy_candidates)} galaxy candidates from image")
    ```

## Step 4: Run host association

=== "Rust"

    ```rust
    use prost_rs::associate::{associate_host, AssociationConfig};
    use prost_rs::types::Transient;

    let transient = Transient::new(transient_ra, transient_dec);

    let config = AssociationConfig::default();
    let result = associate_host(&transient, &galaxy_candidates, &config)
        .expect("association failed");

    println!("Candidates: {}", result.candidates.len());
    println!("P(none):    {:.4}", result.p_none);

    if let Some(best) = result.best_host() {
        println!("Best host:      {:?}", best.galaxy.objname);
        println!("  Separation:   {:.1}\"", best.separation_arcsec);
        println!("  DLR rank:     {}", best.dlr_rank);
        println!("  Frac offset:  {:.3}", best.fractional_offset);
        println!("  P(host):      {:.4}", best.posterior);
        println!("  a={:.2}\" b={:.2}\"", best.galaxy.a_arcsec, best.galaxy.b_arcsec);
    }
    ```

=== "Python"

    ```python
    from prost import Transient, associate_host

    transient = Transient(ra=transient_ra, dec=transient_dec, name="SN 2020jfo")
    result = associate_host(transient, galaxy_candidates)

    print(f"Candidates: {len(result.candidates)}")
    print(f"P(none):    {result.p_none:.4f}")

    if result.best_host:
        h = result.best_host
        print(f"Best host:      {h.galaxy.objname}")
        print(f"  Separation:   {h.separation_arcsec:.1f}\"")
        print(f"  DLR rank:     {h.dlr_rank}")
        print(f"  Frac offset:  {h.fractional_offset:.3f}")
        print(f"  P(host):      {h.posterior:.4f}")
        print(f"  a={h.galaxy.a_arcsec:.2f}\" b={h.galaxy.b_arcsec:.2f}\"")
    ```

Expected output: the brightest extended source (M61's core) is identified
as the most probable host, with SN 2020jfo well within its disk
(fractional offset < 1).

## Step 5: Combine with catalog redshifts

For best results, combine image-derived morphology with catalog redshifts
from BOOM or NED:

=== "Rust"

    ```rust
    // Enrich the best candidate with a known redshift
    if let Some(best) = result.best_host() {
        for g in galaxy_candidates.iter_mut() {
            if g.ra == best.galaxy.ra && g.dec == best.galaxy.dec {
                g.redshift = Some(0.00522);
                g.redshift_err = Some(0.00001);
                g.objname = Some("NGC 4303 (M61)".to_string());
            }
        }
    }

    // Re-run with redshift scoring
    let mut transient_z = Transient::new(transient_ra, transient_dec)
        .with_redshift(0.00522, 0.0001);

    let config_z = AssociationConfig {
        use_redshift: true,
        ..Default::default()
    };
    let result_z = associate_host(&transient_z, &galaxy_candidates, &config_z)
        .expect("redshift-aware association failed");

    if let Some(best) = result_z.best_host() {
        println!(
            "With redshift: P(host)={:.4} (offset={:.3}, z={:.3})",
            best.posterior, best.posterior_offset, best.posterior_redshift,
        );
    }
    ```

=== "Python"

    ```python
    from prost import AssociationConfig

    if result.best_host:
        for g in galaxy_candidates:
            if g.ra == result.best_host.galaxy.ra:
                g.redshift = 0.00522
                g.redshift_err = 0.00001
                g.objname = "NGC 4303 (M61)"

    config = AssociationConfig(use_redshift=True)
    transient.redshift = 0.00522
    transient.redshift_err = 0.0001

    result_z = associate_host(transient, galaxy_candidates, config)

    if result_z.best_host:
        h = result_z.best_host
        print(f"With redshift: P(host)={h.posterior:.4f} "
              f"(offset={h.posterior_offset:.3f}, z={h.posterior_redshift:.3f})")
    ```

!!! tip "When to use image-based vs catalog-based association"
    **Catalog-based** (BOOM/LS_DR10 cross-match) is faster and includes
    redshifts, magnitudes, and Tractor shape parameters out of the box.
    Use it when a catalog covers your field.

    **Image-based** extraction is essential when:

    - No catalog coverage exists (e.g., Southern sky gaps, faint hosts)
    - Catalog shapes are unreliable (blended sources, shredded objects)
    - You need Sersic profile parameters (n, r_eff) for science
    - The host is extended enough that Tractor PSF models fail
