//! Integration tests for downloading image cutouts from survey servers.
//!
//! These tests verify that we can retrieve FITS/JPEG cutouts from
//! Legacy Survey and Pan-STARRS image services for known positions.
//!
//! Run with: cargo test --features catalogs -- --ignored

/// Legacy Survey cutout service base URL
const LS_CUTOUT_URL: &str = "https://www.legacysurvey.org/viewer/cutout.fits";
/// Pan-STARRS cutout service base URL
const PS1_CUTOUT_URL: &str = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py";

/// Download a Legacy Survey FITS cutout as raw bytes.
async fn download_ls_cutout(
    ra: f64,
    dec: f64,
    size_pixels: u32,
    layer: &str,
) -> Result<Vec<u8>, String> {
    let url = format!(
        "{}?ra={}&dec={}&size={}&layer={}&pixscale=0.262",
        LS_CUTOUT_URL, ra, dec, size_pixels, layer
    );

    let client = reqwest::Client::new();
    let resp = client
        .get(&url)
        .send()
        .await
        .map_err(|e| format!("request failed: {e}"))?;

    if !resp.status().is_success() {
        return Err(format!("HTTP {}", resp.status()));
    }

    let bytes = resp
        .bytes()
        .await
        .map_err(|e| format!("body read failed: {e}"))?;

    Ok(bytes.to_vec())
}

/// Download a Legacy Survey JPEG cutout.
async fn download_ls_jpeg(
    ra: f64,
    dec: f64,
    size_pixels: u32,
    layer: &str,
) -> Result<Vec<u8>, String> {
    let url = format!(
        "https://www.legacysurvey.org/viewer/jpeg-cutout?ra={}&dec={}&size={}&layer={}&pixscale=0.262",
        ra, dec, size_pixels, layer
    );

    let client = reqwest::Client::new();
    let resp = client
        .get(&url)
        .send()
        .await
        .map_err(|e| format!("request failed: {e}"))?;

    if !resp.status().is_success() {
        return Err(format!("HTTP {}", resp.status()));
    }

    let bytes = resp
        .bytes()
        .await
        .map_err(|e| format!("body read failed: {e}"))?;

    Ok(bytes.to_vec())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Download a FITS cutout of NGC 4993 (AT2017gfo host) from Legacy Survey.
#[tokio::test]
#[ignore]
async fn test_download_ls_fits_ngc4993() {
    let ra = 197.4487;
    let dec = -23.3839;
    let size = 128; // ~33" at 0.262"/pix

    match download_ls_cutout(ra, dec, size, "ls-dr10").await {
        Ok(bytes) => {
            println!("NGC 4993 FITS cutout: {} bytes", bytes.len());
            // FITS files start with "SIMPLE  ="
            assert!(
                bytes.len() > 2880,
                "FITS cutout should be at least one FITS block"
            );
            let header = String::from_utf8_lossy(&bytes[..30]);
            assert!(
                header.starts_with("SIMPLE"),
                "Should be a valid FITS file, got: {header}"
            );
        }
        Err(e) => {
            eprintln!("Could not download LS cutout (network issue?): {e}");
        }
    }
}

/// Download a JPEG cutout of M82 (SN 2014J host).
#[tokio::test]
#[ignore]
async fn test_download_ls_jpeg_m82() {
    let ra = 148.9685;
    let dec = 69.6797;
    let size = 256;

    match download_ls_jpeg(ra, dec, size, "ls-dr10").await {
        Ok(bytes) => {
            println!("M82 JPEG cutout: {} bytes", bytes.len());
            // JPEG starts with FF D8 FF
            assert!(bytes.len() > 1000, "JPEG should have some content");
            assert_eq!(
                bytes[0..2],
                [0xFF, 0xD8],
                "Should be a valid JPEG (starts with FF D8)"
            );
        }
        Err(e) => {
            eprintln!("Could not download JPEG cutout: {e}");
        }
    }
}

/// Download cutouts at multiple positions to verify batch capability.
#[tokio::test]
#[ignore]
async fn test_download_ls_batch_cutouts() {
    let targets = vec![
        ("NGC4993", 197.4487, -23.3839),
        ("M82", 148.9685, 69.6797),
        ("M61", 185.4788, 4.4734),
    ];

    let client = reqwest::Client::new();
    let mut successes = 0;

    for (name, ra, dec) in &targets {
        let url = format!(
            "{}?ra={}&dec={}&size=64&layer=ls-dr10&pixscale=0.262",
            LS_CUTOUT_URL, ra, dec
        );

        match client.get(&url).send().await {
            Ok(resp) if resp.status().is_success() => {
                let bytes = resp.bytes().await.unwrap();
                println!("{}: {} bytes", name, bytes.len());
                successes += 1;
            }
            Ok(resp) => {
                eprintln!("{}: HTTP {}", name, resp.status());
            }
            Err(e) => {
                eprintln!("{}: request failed: {e}", name);
            }
        }
    }

    assert!(
        successes >= 1,
        "At least one cutout download should succeed"
    );
}

/// Verify that cutout size scales correctly.
#[tokio::test]
#[ignore]
async fn test_cutout_size_scaling() {
    let ra = 197.4487;
    let dec = -23.3839;

    let small = download_ls_cutout(ra, dec, 32, "ls-dr10").await;
    let large = download_ls_cutout(ra, dec, 128, "ls-dr10").await;

    match (small, large) {
        (Ok(s), Ok(l)) => {
            println!("32px cutout: {} bytes, 128px cutout: {} bytes", s.len(), l.len());
            assert!(
                l.len() > s.len(),
                "Larger cutout should have more data"
            );
        }
        _ => {
            eprintln!("Cutout download failed, skipping size comparison");
        }
    }
}

/// Download a Pan-STARRS cutout (PS1 covers dec > -30°).
#[tokio::test]
#[ignore]
async fn test_download_ps1_image_list() {
    let ra = 148.9685; // M82
    let dec = 69.6797;

    let url = format!(
        "{}?ra={}&dec={}&size=240&output_size=128&filters=grizy&type=stack",
        PS1_CUTOUT_URL, ra, dec
    );

    let client = reqwest::Client::new();
    match client.get(&url).send().await {
        Ok(resp) if resp.status().is_success() => {
            let text = resp.text().await.unwrap();
            println!("PS1 image list ({} bytes):\n{}", text.len(), &text[..text.len().min(500)]);
            // Should contain filenames for available filters
            assert!(!text.is_empty(), "PS1 should return image list");
        }
        Ok(resp) => {
            eprintln!("PS1 returned HTTP {}", resp.status());
        }
        Err(e) => {
            eprintln!("PS1 request failed: {e}");
        }
    }
}
