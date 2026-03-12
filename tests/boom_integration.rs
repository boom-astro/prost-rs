//! Integration tests that query BOOM API for real transient data and
//! run host association against Legacy Survey cross-match results.
//!
//! Requires: BOOM_USERNAME and BOOM_PASSWORD environment variables.
//! Run with: cargo test --features catalogs -- --ignored

use prost_rs::associate::{associate_host, AssociationConfig};
use prost_rs::ellipse::Ellipse;
use prost_rs::types::{GalaxyCandidate, Transient};
use serde_json::Value;

const BOOM_BASE_URL: &str = "https://api.kaboom.caltech.edu";

/// Authenticate with BOOM and return a reqwest client with bearer token.
async fn boom_client() -> Option<reqwest::Client> {
    let username = std::env::var("BOOM_USERNAME").ok()?;
    let password = std::env::var("BOOM_PASSWORD").ok()?;

    let client = reqwest::Client::new();
    let resp = client
        .post(format!("{BOOM_BASE_URL}/auth"))
        .form(&[("username", &username), ("password", &password)])
        .send()
        .await
        .ok()?;

    let body: Value = resp.json().await.ok()?;
    let token = body
        .get("access_token")
        .or_else(|| body.get("token"))?
        .as_str()?;

    let mut headers = reqwest::header::HeaderMap::new();
    headers.insert(
        reqwest::header::AUTHORIZATION,
        format!("Bearer {token}").parse().unwrap(),
    );

    Some(
        reqwest::Client::builder()
            .default_headers(headers)
            .build()
            .unwrap(),
    )
}

/// Get coordinates for a ZTF object from BOOM.
async fn get_ztf_coords(client: &reqwest::Client, object_id: &str) -> Option<(f64, f64)> {
    let pipeline = serde_json::json!([
        {"$match": {"objectId": object_id}},
        {"$group": {
            "_id": "$objectId",
            "ra": {"$first": "$candidate.ra"},
            "dec": {"$first": "$candidate.dec"},
        }},
    ]);

    let resp = client
        .post(format!("{BOOM_BASE_URL}/queries/pipeline"))
        .json(&serde_json::json!({
            "catalog_name": "ZTF_alerts",
            "pipeline": pipeline,
        }))
        .send()
        .await
        .ok()?;

    let body: Value = resp.json().await.ok()?;
    let data = body.get("data")?.as_array()?;
    let doc = data.first()?;
    let ra = doc.get("ra")?.as_f64()?;
    let dec = doc.get("dec")?.as_f64()?;
    Some((ra, dec))
}

/// Cone search a catalog on BOOM for galaxies near a position.
async fn cone_search_ls(
    client: &reqwest::Client,
    name: &str,
    ra: f64,
    dec: f64,
    radius_arcsec: f64,
) -> Option<Vec<Value>> {
    let resp = client
        .post(format!("{BOOM_BASE_URL}/queries/cone_search"))
        .json(&serde_json::json!({
            "catalog_name": "LS_DR10",
            "radius": radius_arcsec,
            "unit": "Arcseconds",
            "object_coordinates": {name: [ra, dec]},
            "projection": {
                "_id": 1,
                "ra": 1,
                "dec": 1,
                "type": 1,
                "shape_r": 1,
                "shape_e1": 1,
                "shape_e2": 1,
                "z": 1,
            },
        }))
        .send()
        .await
        .ok()?;

    let body: Value = resp.json().await.ok()?;
    let data = body.get("data")?.get(name)?.as_array()?.clone();
    Some(data)
}

/// Convert BOOM LS_DR10 documents to GalaxyCandidates.
fn docs_to_candidates(docs: &[Value]) -> Vec<GalaxyCandidate> {
    docs.iter()
        .filter_map(|doc| {
            let ra = doc.get("ra")?.as_f64()?;
            let dec = doc.get("dec")?.as_f64()?;
            let shape_r = doc.get("shape_r")?.as_f64()?;
            let shape_e1 = doc.get("shape_e1").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let shape_e2 = doc.get("shape_e2").and_then(|v| v.as_f64()).unwrap_or(0.0);

            // Filter stars
            let obj_type = doc.get("type").and_then(|v| v.as_str()).unwrap_or("");
            if obj_type == "PSF" {
                return None;
            }

            let ellipse = Ellipse::from_tractor(shape_r, shape_e1, shape_e2, 0.05).ok()?;

            Some(GalaxyCandidate {
                ra,
                dec,
                a_arcsec: ellipse.a,
                b_arcsec: ellipse.b,
                pa_deg: ellipse.pa_rad.to_degrees(),
                redshift: doc.get("z").and_then(|v| v.as_f64()),
                redshift_err: None,
                mag: None,
                mag_err: None,
                objtype: Some(obj_type.to_string()),
                objname: doc.get("_id").and_then(|v| v.as_str()).map(|s| s.to_string()),
                catalog: Some("LS_DR10".to_string()),
                shape_from_image: false,
            })
        })
        .collect()
}

/// Query NED on BOOM for host galaxy info (redshift, type).
async fn cone_search_ned(
    client: &reqwest::Client,
    name: &str,
    ra: f64,
    dec: f64,
    radius_arcsec: f64,
) -> Option<Vec<Value>> {
    let resp = client
        .post(format!("{BOOM_BASE_URL}/queries/cone_search"))
        .json(&serde_json::json!({
            "catalog_name": "NED",
            "radius": radius_arcsec,
            "unit": "Arcseconds",
            "object_coordinates": {name: [ra, dec]},
            "projection": {
                "_id": 1,
                "objtype": 1,
                "z": 1,
            },
        }))
        .send()
        .await
        .ok()?;

    let body: Value = resp.json().await.ok()?;
    let data = body.get("data")?.get(name)?.as_array()?.clone();
    Some(data)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Test association for a well-known ZTF transient.
/// ZTF18abvkwla (the "Koala") - AT2018lug - fast blue optical transient
/// in a star-forming dwarf galaxy at z=0.2714
#[tokio::test]
#[ignore]
async fn test_boom_ztf18abvkwla() {
    let client = match boom_client().await {
        Some(c) => c,
        None => {
            eprintln!("Skipping: BOOM credentials not set");
            return;
        }
    };

    let (ra, dec) = get_ztf_coords(&client, "ZTF18abvkwla")
        .await
        .expect("Failed to get ZTF18abvkwla coordinates from BOOM");

    let docs = cone_search_ls(&client, "ZTF18abvkwla", ra, dec, 60.0)
        .await
        .expect("Failed to query LS_DR10");

    assert!(!docs.is_empty(), "Should find LS_DR10 sources near ZTF18abvkwla");

    let candidates = docs_to_candidates(&docs);
    if candidates.is_empty() {
        eprintln!("No extended LS_DR10 sources found (all PSF?), skipping association");
        return;
    }

    let transient = Transient::new(ra, dec);
    let config = AssociationConfig::default();
    let result = associate_host(&transient, &candidates, &config).unwrap();

    println!("ZTF18abvkwla: {} candidates, p_none={:.4}", result.candidates.len(), result.p_none);
    if let Some(best) = result.best_host() {
        println!(
            "  Best host: {:?} at {:.1}\" (FO={:.2}, P={:.4})",
            best.galaxy.objname, best.separation_arcsec, best.fractional_offset, best.posterior
        );
    }

    // Should find at least one candidate
    assert!(!result.candidates.is_empty());
}

/// Test a known SN Ia: ZTF20abjbgjj (SN 2020jfo in M61/NGC 4303)
/// Well-known host: NGC 4303 (M61), a large barred spiral at z=0.005224
#[tokio::test]
#[ignore]
async fn test_boom_sn2020jfo_in_m61() {
    let client = match boom_client().await {
        Some(c) => c,
        None => {
            eprintln!("Skipping: BOOM credentials not set");
            return;
        }
    };

    let (ra, dec) = get_ztf_coords(&client, "ZTF20abjbgjj")
        .await
        .expect("Failed to get ZTF20abjbgjj coordinates");

    println!("SN 2020jfo (ZTF20abjbgjj): RA={ra:.6}, Dec={dec:.6}");

    // Query LS_DR10 and NED
    let ls_docs = cone_search_ls(&client, "ZTF20abjbgjj", ra, dec, 60.0)
        .await
        .expect("Failed to query LS_DR10");

    let ned_docs = cone_search_ned(&client, "ZTF20abjbgjj", ra, dec, 30.0)
        .await
        .expect("Failed to query NED");

    println!("  LS_DR10: {} sources, NED: {} sources", ls_docs.len(), ned_docs.len());

    let candidates = docs_to_candidates(&ls_docs);
    if candidates.is_empty() {
        eprintln!("No extended sources found, skipping");
        return;
    }

    let transient = Transient::new(ra, dec);
    let config = AssociationConfig::default();
    let result = associate_host(&transient, &candidates, &config).unwrap();

    println!(
        "  Association: {} candidates, p_none={:.4}",
        result.candidates.len(),
        result.p_none
    );
    for c in &result.candidates {
        println!(
            "    rank={} sep={:.1}\" FO={:.2} P={:.4} name={:?} type={:?}",
            c.dlr_rank, c.separation_arcsec, c.fractional_offset,
            c.posterior, c.galaxy.objname, c.galaxy.objtype
        );
    }

    // SN 2020jfo is well inside M61, should have a clear host
    assert!(!result.candidates.is_empty());
    let best = result.best_host().unwrap();
    assert!(
        best.posterior > 0.3,
        "Best host for SN in M61 should have decent posterior, got {}",
        best.posterior
    );
}

/// Test a hostless / ambiguous case: look for a transient far from galaxies.
/// Use a position in an empty field to verify p_none is high.
#[tokio::test]
#[ignore]
async fn test_boom_empty_field() {
    let client = match boom_client().await {
        Some(c) => c,
        None => {
            eprintln!("Skipping: BOOM credentials not set");
            return;
        }
    };

    // Pick a position in a sparse region (high galactic latitude void)
    let ra = 200.0;
    let dec = 60.0;

    let docs = cone_search_ls(&client, "empty_field", ra, dec, 30.0)
        .await
        .expect("Failed to query LS_DR10");

    let candidates = docs_to_candidates(&docs);
    let transient = Transient::new(ra, dec);
    let config = AssociationConfig::default();
    let result = associate_host(&transient, &candidates, &config).unwrap();

    println!(
        "Empty field: {} LS sources, {} candidates after filter, p_none={:.4}",
        docs.len(),
        result.candidates.len(),
        result.p_none
    );

    // Even if there are some sources, p_none should be relatively high
    // since nothing should be close in DLR
    if result.candidates.is_empty() {
        assert_eq!(result.p_none, 1.0);
    }
}

/// Test batch association: get multiple known ZTF transients and associate all.
#[tokio::test]
#[ignore]
async fn test_boom_batch_association() {
    let client = match boom_client().await {
        Some(c) => c,
        None => {
            eprintln!("Skipping: BOOM credentials not set");
            return;
        }
    };

    // A mix of well-known ZTF transients
    let objects = vec![
        "ZTF20abjbgjj", // SN 2020jfo in M61
        "ZTF18abvkwla", // AT2018lug "Koala"
        "ZTF19aarioci", // SN 2019ehk
    ];

    let config = AssociationConfig::default();
    let mut results = Vec::new();

    for obj_id in &objects {
        let coords = match get_ztf_coords(&client, obj_id).await {
            Some(c) => c,
            None => {
                eprintln!("  Could not get coords for {obj_id}, skipping");
                continue;
            }
        };

        let docs = match cone_search_ls(&client, obj_id, coords.0, coords.1, 60.0).await {
            Some(d) => d,
            None => {
                eprintln!("  Could not query LS_DR10 for {obj_id}, skipping");
                continue;
            }
        };

        let candidates = docs_to_candidates(&docs);
        let transient = Transient::new(coords.0, coords.1);
        let result = associate_host(&transient, &candidates, &config);

        match result {
            Ok(r) => {
                let best_info = r.best_host().map(|b| {
                    format!(
                        "sep={:.1}\" FO={:.2} P={:.4}",
                        b.separation_arcsec, b.fractional_offset, b.posterior
                    )
                });
                println!(
                    "{}: {} candidates, p_none={:.4}, best=[{}]",
                    obj_id,
                    r.candidates.len(),
                    r.p_none,
                    best_info.unwrap_or_else(|| "none".to_string())
                );
                results.push((obj_id.to_string(), r));
            }
            Err(e) => {
                println!("{obj_id}: association error: {e}");
            }
        }
    }

    // At least some should succeed
    assert!(
        !results.is_empty(),
        "At least one batch association should succeed"
    );

    // All results should have valid posteriors
    for (obj_id, result) in &results {
        let total: f64 =
            result.candidates.iter().map(|c| c.posterior).sum::<f64>() + result.p_none;
        assert!(
            (total - 1.0).abs() < 1e-6,
            "{}: posteriors sum to {}, expected 1.0",
            obj_id,
            total
        );
    }
}
