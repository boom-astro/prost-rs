//! BOOM alert broker client for querying transients and galaxy catalogs.
//!
//! Provides authenticated access to the BOOM API for retrieving transient
//! coordinates, Legacy Survey cross-matches, and NED catalog data.
//!
//! # Configuration
//!
//! `BoomConfig` reads credentials from environment variables by default
//! (`BOOM_USERNAME`, `BOOM_PASSWORD`). Without credentials, the client
//! can still be constructed but authenticated endpoints will fail.
//!
//! # Feature gate
//!
//! This module requires the `catalogs` feature:
//! ```toml
//! prost-rs = { git = "...", features = ["catalogs"] }
//! ```

use crate::ellipse::Ellipse;
use crate::errors::ProstError;
use crate::types::GalaxyCandidate;

/// Default BOOM API base URL.
pub const DEFAULT_BOOM_URL: &str = "https://api.kaboom.caltech.edu";

/// Configuration for the BOOM client.
#[derive(Debug, Clone)]
pub struct BoomConfig {
    /// BOOM API base URL.
    pub base_url: String,
    /// Username for authentication (None = anonymous / read from env).
    pub username: Option<String>,
    /// Password for authentication (None = anonymous / read from env).
    pub password: Option<String>,
    /// Default search radius for cone searches (arcsec).
    pub default_radius_arcsec: f64,
    /// Default catalog for galaxy cross-matches.
    pub default_catalog: String,
    /// Minimum semi-minor axis floor for Tractor conversion (arcsec).
    pub min_b_arcsec: f64,
}

impl Default for BoomConfig {
    fn default() -> Self {
        Self {
            base_url: DEFAULT_BOOM_URL.to_string(),
            username: std::env::var("BOOM_USERNAME").ok(),
            password: std::env::var("BOOM_PASSWORD").ok(),
            default_radius_arcsec: 60.0,
            default_catalog: "LS_DR10".to_string(),
            min_b_arcsec: 0.05,
        }
    }
}

impl BoomConfig {
    /// Create a config with explicit credentials.
    pub fn with_credentials(username: &str, password: &str) -> Self {
        Self {
            username: Some(username.to_string()),
            password: Some(password.to_string()),
            ..Default::default()
        }
    }

    /// Create a config that reads credentials from environment variables.
    /// Returns a config with `username` and `password` set to `None` if
    /// the env vars are not set — the client will fail on auth endpoints.
    pub fn from_env() -> Self {
        Self::default()
    }

    /// Check whether credentials are available.
    pub fn has_credentials(&self) -> bool {
        self.username.is_some() && self.password.is_some()
    }
}

/// Async BOOM API client.
///
/// Handles authentication, transient lookups, and catalog cone searches.
pub struct BoomClient {
    config: BoomConfig,
    http: reqwest::Client,
    token: Option<String>,
}

impl BoomClient {
    /// Create a new client. Does not authenticate until `authenticate()` is called.
    pub fn new(config: BoomConfig) -> Self {
        Self {
            config,
            http: reqwest::Client::new(),
            token: None,
        }
    }

    /// Create a client with default config (credentials from env).
    pub fn from_env() -> Self {
        Self::new(BoomConfig::from_env())
    }

    /// Authenticate with BOOM and store the bearer token.
    pub async fn authenticate(&mut self) -> Result<(), ProstError> {
        let username = self.config.username.as_deref().ok_or_else(|| {
            ProstError::BoomError("BOOM_USERNAME not set".to_string())
        })?;
        let password = self.config.password.as_deref().ok_or_else(|| {
            ProstError::BoomError("BOOM_PASSWORD not set".to_string())
        })?;

        let resp = self
            .http
            .post(format!("{}/auth", self.config.base_url))
            .form(&[("username", username), ("password", password)])
            .send()
            .await
            .map_err(|e| ProstError::BoomError(format!("auth request failed: {e}")))?;

        if !resp.status().is_success() {
            return Err(ProstError::BoomError(format!(
                "auth failed: HTTP {}",
                resp.status()
            )));
        }

        let body: serde_json::Value = resp
            .json()
            .await
            .map_err(|e| ProstError::BoomError(format!("auth response parse failed: {e}")))?;

        let token = body
            .get("access_token")
            .or_else(|| body.get("token"))
            .and_then(|v| v.as_str())
            .ok_or_else(|| ProstError::BoomError("no token in auth response".to_string()))?;

        self.token = Some(token.to_string());

        // Rebuild HTTP client with default auth header
        let mut headers = reqwest::header::HeaderMap::new();
        headers.insert(
            reqwest::header::AUTHORIZATION,
            format!("Bearer {token}")
                .parse()
                .map_err(|e| ProstError::BoomError(format!("invalid token: {e}")))?,
        );
        self.http = reqwest::Client::builder()
            .default_headers(headers)
            .build()
            .map_err(|e| ProstError::BoomError(format!("client build failed: {e}")))?;

        Ok(())
    }

    /// Check whether the client has authenticated.
    pub fn is_authenticated(&self) -> bool {
        self.token.is_some()
    }

    /// Get coordinates for a ZTF object by its objectId.
    ///
    /// Queries the ZTF_alerts catalog using an aggregation pipeline.
    pub async fn get_ztf_coordinates(&self, object_id: &str) -> Result<(f64, f64), ProstError> {
        let pipeline = serde_json::json!([
            {"$match": {"objectId": object_id}},
            {"$group": {
                "_id": "$objectId",
                "ra": {"$first": "$candidate.ra"},
                "dec": {"$first": "$candidate.dec"},
            }},
        ]);

        let resp = self
            .http
            .post(format!("{}/queries/pipeline", self.config.base_url))
            .json(&serde_json::json!({
                "catalog_name": "ZTF_alerts",
                "pipeline": pipeline,
            }))
            .send()
            .await
            .map_err(|e| ProstError::BoomError(format!("pipeline query failed: {e}")))?;

        if !resp.status().is_success() {
            return Err(ProstError::BoomError(format!(
                "pipeline query: HTTP {}",
                resp.status()
            )));
        }

        let body: serde_json::Value = resp
            .json()
            .await
            .map_err(|e| ProstError::BoomError(format!("pipeline response parse failed: {e}")))?;

        let data = body
            .get("data")
            .and_then(|d| d.as_array())
            .ok_or_else(|| ProstError::BoomError("no data in pipeline response".to_string()))?;

        let doc = data
            .first()
            .ok_or_else(|| ProstError::BoomError(format!("object {object_id} not found")))?;

        let ra = doc
            .get("ra")
            .and_then(|v| v.as_f64())
            .ok_or_else(|| ProstError::BoomError("missing ra in response".to_string()))?;
        let dec = doc
            .get("dec")
            .and_then(|v| v.as_f64())
            .ok_or_else(|| ProstError::BoomError("missing dec in response".to_string()))?;

        Ok((ra, dec))
    }

    /// Cone search a catalog on BOOM.
    ///
    /// Returns raw JSON documents from the catalog. Use `docs_to_candidates`
    /// to convert LS_DR10 results to `GalaxyCandidate`s.
    pub async fn cone_search(
        &self,
        name: &str,
        ra: f64,
        dec: f64,
        radius_arcsec: f64,
        catalog: &str,
        projection: &serde_json::Value,
    ) -> Result<Vec<serde_json::Value>, ProstError> {
        let resp = self
            .http
            .post(format!("{}/queries/cone_search", self.config.base_url))
            .json(&serde_json::json!({
                "catalog_name": catalog,
                "radius": radius_arcsec,
                "unit": "Arcseconds",
                "object_coordinates": {name: [ra, dec]},
                "projection": projection,
            }))
            .send()
            .await
            .map_err(|e| ProstError::BoomError(format!("cone search failed: {e}")))?;

        if !resp.status().is_success() {
            return Err(ProstError::BoomError(format!(
                "cone search: HTTP {}",
                resp.status()
            )));
        }

        let body: serde_json::Value = resp
            .json()
            .await
            .map_err(|e| ProstError::BoomError(format!("cone search parse failed: {e}")))?;

        let data = body
            .get("data")
            .and_then(|d| d.get(name))
            .and_then(|d| d.as_array())
            .cloned()
            .unwrap_or_default();

        Ok(data)
    }

    /// Cone search Legacy Survey (LS_DR10) for galaxies near a position.
    ///
    /// Convenience wrapper around `cone_search` with LS-specific projection.
    pub async fn cone_search_ls(
        &self,
        name: &str,
        ra: f64,
        dec: f64,
        radius_arcsec: Option<f64>,
    ) -> Result<Vec<serde_json::Value>, ProstError> {
        let radius = radius_arcsec.unwrap_or(self.config.default_radius_arcsec);
        let projection = serde_json::json!({
            "_id": 1, "ra": 1, "dec": 1, "type": 1,
            "shape_r": 1, "shape_e1": 1, "shape_e2": 1, "z": 1,
        });
        self.cone_search(name, ra, dec, radius, &self.config.default_catalog, &projection)
            .await
    }

    /// Cone search NED for galaxy redshifts and types.
    pub async fn cone_search_ned(
        &self,
        name: &str,
        ra: f64,
        dec: f64,
        radius_arcsec: Option<f64>,
    ) -> Result<Vec<serde_json::Value>, ProstError> {
        let radius = radius_arcsec.unwrap_or(30.0);
        let projection = serde_json::json!({
            "_id": 1, "objtype": 1, "z": 1,
        });
        self.cone_search(name, ra, dec, radius, "NED", &projection)
            .await
    }

    /// Query BOOM for a ZTF transient and return galaxy candidates.
    ///
    /// This is the high-level convenience method: looks up coordinates,
    /// runs a Legacy Survey cone search, and converts results to
    /// `GalaxyCandidate`s ready for `associate_host`.
    pub async fn get_candidates(
        &self,
        object_id: &str,
        radius_arcsec: Option<f64>,
    ) -> Result<((f64, f64), Vec<GalaxyCandidate>), ProstError> {
        let (ra, dec) = self.get_ztf_coordinates(object_id).await?;
        let docs = self.cone_search_ls(object_id, ra, dec, radius_arcsec).await?;
        let candidates = docs_to_candidates(&docs, self.config.min_b_arcsec);
        Ok(((ra, dec), candidates))
    }
}

/// Convert BOOM LS_DR10 JSON documents to `GalaxyCandidate`s.
///
/// Filters out PSF (star) sources and converts Tractor shape parameters
/// to ellipse parameters. Documents missing `shape_r` are skipped.
pub fn docs_to_candidates(docs: &[serde_json::Value], min_b: f64) -> Vec<GalaxyCandidate> {
    docs.iter()
        .filter_map(|doc| {
            let ra = doc.get("ra")?.as_f64()?;
            let dec = doc.get("dec")?.as_f64()?;
            let shape_r = doc.get("shape_r")?.as_f64()?;
            let shape_e1 = doc
                .get("shape_e1")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let shape_e2 = doc
                .get("shape_e2")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);

            // Filter stars
            let obj_type = doc.get("type").and_then(|v| v.as_str()).unwrap_or("");
            if obj_type == "PSF" {
                return None;
            }

            let ellipse = Ellipse::from_tractor(shape_r, shape_e1, shape_e2, min_b).ok()?;

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
                objname: doc
                    .get("_id")
                    .and_then(|v| v.as_str())
                    .map(|s| s.to_string()),
                catalog: Some("LS_DR10".to_string()),
                shape_from_image: false,
            })
        })
        .collect()
}
