use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use prost_rs::associate::{
    associate_host as rs_associate_host, AssociationConfig as RsConfig,
};
use prost_rs::cutout::{Cutout as RsCutout, ImagingSurvey};
use prost_rs::dlr::compute_dlr as rs_compute_dlr;
use prost_rs::ellipse::Ellipse as RsEllipse;
use prost_rs::likelihood::{offset_likelihood as rs_offset_likelihood, redshift_likelihood as rs_redshift_likelihood};
use prost_rs::morphology::{
    elliptical_radius as rs_elliptical_radius, fit_sersic as rs_fit_sersic,
    sersic_profile as rs_sersic_profile, FitConfig as RsFitConfig,
    SersicFit as RsSersicFit,
};
use prost_rs::source::{
    estimate_background as rs_estimate_background, extract_sources as rs_extract_sources,
    Background as RsBackground, DetectedSource as RsDetectedSource,
    ExtractionConfig as RsExtractionConfig,
};
use prost_rs::types::{
    GalaxyCandidate as RsGalaxy, HostCandidate as RsHost, Transient as RsTransient,
};

// ---------------------------------------------------------------------------
// Transient
// ---------------------------------------------------------------------------

#[pyclass(name = "Transient")]
#[derive(Clone)]
pub struct PyTransient {
    #[pyo3(get, set)]
    pub ra: f64,
    #[pyo3(get, set)]
    pub dec: f64,
    #[pyo3(get, set)]
    pub ra_err: f64,
    #[pyo3(get, set)]
    pub dec_err: f64,
    #[pyo3(get, set)]
    pub redshift: Option<f64>,
    #[pyo3(get, set)]
    pub redshift_err: Option<f64>,
    #[pyo3(get, set)]
    pub name: Option<String>,
}

#[pymethods]
impl PyTransient {
    #[new]
    #[pyo3(signature = (ra, dec, ra_err=0.0, dec_err=0.0, redshift=None, redshift_err=None, name=None))]
    fn new(
        ra: f64,
        dec: f64,
        ra_err: f64,
        dec_err: f64,
        redshift: Option<f64>,
        redshift_err: Option<f64>,
        name: Option<String>,
    ) -> Self {
        Self {
            ra,
            dec,
            ra_err,
            dec_err,
            redshift,
            redshift_err,
            name,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Transient(ra={:.6}, dec={:.6}, name={:?})",
            self.ra, self.dec, self.name
        )
    }
}

impl PyTransient {
    fn to_rust(&self) -> RsTransient {
        let mut t = RsTransient::new(self.ra, self.dec)
            .with_position_err(self.ra_err, self.dec_err);
        if let (Some(z), Some(ze)) = (self.redshift, self.redshift_err) {
            t = t.with_redshift(z, ze);
        } else if let Some(z) = self.redshift {
            t = t.with_redshift(z, 0.0);
        }
        t.name = self.name.clone();
        t
    }
}

// ---------------------------------------------------------------------------
// GalaxyCandidate
// ---------------------------------------------------------------------------

#[pyclass(name = "GalaxyCandidate")]
#[derive(Clone)]
pub struct PyGalaxyCandidate {
    #[pyo3(get, set)]
    pub ra: f64,
    #[pyo3(get, set)]
    pub dec: f64,
    #[pyo3(get, set)]
    pub a_arcsec: f64,
    #[pyo3(get, set)]
    pub b_arcsec: f64,
    #[pyo3(get, set)]
    pub pa_deg: f64,
    #[pyo3(get, set)]
    pub redshift: Option<f64>,
    #[pyo3(get, set)]
    pub redshift_err: Option<f64>,
    #[pyo3(get, set)]
    pub mag: Option<f64>,
    #[pyo3(get, set)]
    pub mag_err: Option<f64>,
    #[pyo3(get, set)]
    pub objtype: Option<String>,
    #[pyo3(get, set)]
    pub objname: Option<String>,
    #[pyo3(get, set)]
    pub catalog: Option<String>,
    #[pyo3(get, set)]
    pub shape_from_image: bool,
}

#[pymethods]
impl PyGalaxyCandidate {
    #[new]
    #[pyo3(signature = (ra, dec, a_arcsec, b_arcsec, pa_deg, redshift=None, redshift_err=None, mag=None, mag_err=None, objtype=None, objname=None, catalog=None, shape_from_image=false))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        ra: f64,
        dec: f64,
        a_arcsec: f64,
        b_arcsec: f64,
        pa_deg: f64,
        redshift: Option<f64>,
        redshift_err: Option<f64>,
        mag: Option<f64>,
        mag_err: Option<f64>,
        objtype: Option<String>,
        objname: Option<String>,
        catalog: Option<String>,
        shape_from_image: bool,
    ) -> Self {
        Self {
            ra,
            dec,
            a_arcsec,
            b_arcsec,
            pa_deg,
            redshift,
            redshift_err,
            mag,
            mag_err,
            objtype,
            objname,
            catalog,
            shape_from_image,
        }
    }

    /// Create a GalaxyCandidate from Tractor shape parameters.
    #[staticmethod]
    #[pyo3(signature = (ra, dec, shape_r, shape_e1, shape_e2, min_b=0.05, redshift=None, objtype=None, objname=None))]
    fn from_tractor(
        ra: f64,
        dec: f64,
        shape_r: f64,
        shape_e1: f64,
        shape_e2: f64,
        min_b: f64,
        redshift: Option<f64>,
        objtype: Option<String>,
        objname: Option<String>,
    ) -> PyResult<Self> {
        let ellipse = RsEllipse::from_tractor(shape_r, shape_e1, shape_e2, min_b)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        Ok(Self {
            ra,
            dec,
            a_arcsec: ellipse.a,
            b_arcsec: ellipse.b,
            pa_deg: ellipse.pa_rad.to_degrees(),
            redshift,
            redshift_err: None,
            mag: None,
            mag_err: None,
            objtype,
            objname,
            catalog: None,
            shape_from_image: false,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "GalaxyCandidate(ra={:.6}, dec={:.6}, a={:.2}\", b={:.2}\", name={:?})",
            self.ra, self.dec, self.a_arcsec, self.b_arcsec, self.objname
        )
    }
}

impl PyGalaxyCandidate {
    fn to_rust(&self) -> RsGalaxy {
        RsGalaxy {
            ra: self.ra,
            dec: self.dec,
            a_arcsec: self.a_arcsec,
            b_arcsec: self.b_arcsec,
            pa_deg: self.pa_deg,
            redshift: self.redshift,
            redshift_err: self.redshift_err,
            mag: self.mag,
            mag_err: self.mag_err,
            objtype: self.objtype.clone(),
            objname: self.objname.clone(),
            catalog: self.catalog.clone(),
            shape_from_image: self.shape_from_image,
        }
    }
}

// ---------------------------------------------------------------------------
// HostCandidate (result)
// ---------------------------------------------------------------------------

#[pyclass(name = "HostCandidate")]
#[derive(Clone)]
pub struct PyHostCandidate {
    #[pyo3(get)]
    pub galaxy: PyGalaxyCandidate,
    #[pyo3(get)]
    pub separation_arcsec: f64,
    #[pyo3(get)]
    pub dlr: f64,
    #[pyo3(get)]
    pub fractional_offset: f64,
    #[pyo3(get)]
    pub dlr_rank: u32,
    #[pyo3(get)]
    pub posterior: f64,
    #[pyo3(get)]
    pub posterior_offset: f64,
    #[pyo3(get)]
    pub posterior_redshift: f64,
    #[pyo3(get)]
    pub posterior_absmag: f64,
}

#[pymethods]
impl PyHostCandidate {
    fn __repr__(&self) -> String {
        format!(
            "HostCandidate(rank={}, sep={:.1}\", FO={:.2}, P={:.4}, name={:?})",
            self.dlr_rank,
            self.separation_arcsec,
            self.fractional_offset,
            self.posterior,
            self.galaxy.objname
        )
    }
}

impl From<RsHost> for PyHostCandidate {
    fn from(h: RsHost) -> Self {
        Self {
            galaxy: PyGalaxyCandidate {
                ra: h.galaxy.ra,
                dec: h.galaxy.dec,
                a_arcsec: h.galaxy.a_arcsec,
                b_arcsec: h.galaxy.b_arcsec,
                pa_deg: h.galaxy.pa_deg,
                redshift: h.galaxy.redshift,
                redshift_err: h.galaxy.redshift_err,
                mag: h.galaxy.mag,
                mag_err: h.galaxy.mag_err,
                objtype: h.galaxy.objtype,
                objname: h.galaxy.objname,
                catalog: h.galaxy.catalog,
                shape_from_image: h.galaxy.shape_from_image,
            },
            separation_arcsec: h.separation_arcsec,
            dlr: h.dlr,
            fractional_offset: h.fractional_offset,
            dlr_rank: h.dlr_rank,
            posterior: h.posterior,
            posterior_offset: h.posterior_offset,
            posterior_redshift: h.posterior_redshift,
            posterior_absmag: h.posterior_absmag,
        }
    }
}

// ---------------------------------------------------------------------------
// AssociationConfig
// ---------------------------------------------------------------------------

#[pyclass(name = "AssociationConfig")]
#[derive(Clone)]
pub struct PyAssociationConfig {
    #[pyo3(get, set)]
    pub max_fractional_offset: f64,
    #[pyo3(get, set)]
    pub min_b_arcsec: f64,
    #[pyo3(get, set)]
    pub max_candidates: usize,
    #[pyo3(get, set)]
    pub use_redshift: bool,
    #[pyo3(get, set)]
    pub use_absmag: bool,
}

#[pymethods]
impl PyAssociationConfig {
    #[new]
    #[pyo3(signature = (max_fractional_offset=10.0, min_b_arcsec=0.05, max_candidates=10, use_redshift=true, use_absmag=false))]
    fn new(
        max_fractional_offset: f64,
        min_b_arcsec: f64,
        max_candidates: usize,
        use_redshift: bool,
        use_absmag: bool,
    ) -> Self {
        Self {
            max_fractional_offset,
            min_b_arcsec,
            max_candidates,
            use_redshift,
            use_absmag,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "AssociationConfig(max_fo={}, min_b={}\", max_cand={}, use_z={})",
            self.max_fractional_offset, self.min_b_arcsec, self.max_candidates, self.use_redshift
        )
    }
}

impl PyAssociationConfig {
    fn to_rust(&self) -> RsConfig {
        RsConfig {
            max_fractional_offset: self.max_fractional_offset,
            min_b_arcsec: self.min_b_arcsec,
            max_candidates: self.max_candidates,
            use_redshift: self.use_redshift,
            use_absmag: self.use_absmag,
        }
    }
}

// ---------------------------------------------------------------------------
// AssociationResult
// ---------------------------------------------------------------------------

#[pyclass(name = "AssociationResult")]
pub struct PyAssociationResult {
    #[pyo3(get)]
    pub candidates: Vec<PyHostCandidate>,
    #[pyo3(get)]
    pub p_none: f64,
    #[pyo3(get)]
    pub n_considered: usize,
}

#[pymethods]
impl PyAssociationResult {
    /// The best host candidate, or None if no candidates found.
    #[getter]
    fn best_host(&self) -> Option<PyHostCandidate> {
        self.candidates.first().cloned()
    }

    fn __repr__(&self) -> String {
        format!(
            "AssociationResult(n_candidates={}, p_none={:.4}, best={:?})",
            self.candidates.len(),
            self.p_none,
            self.candidates.first().map(|c| format!(
                "rank {} P={:.3}",
                c.dlr_rank, c.posterior
            ))
        )
    }

    fn __len__(&self) -> usize {
        self.candidates.len()
    }
}

// ---------------------------------------------------------------------------
// Top-level functions
// ---------------------------------------------------------------------------

/// Associate a transient with its most likely host galaxy.
///
/// Parameters
/// ----------
/// transient : Transient
///     The transient event to associate.
/// candidates : list[GalaxyCandidate]
///     List of galaxy candidates to consider.
/// config : AssociationConfig, optional
///     Configuration parameters. Uses defaults if not provided.
///
/// Returns
/// -------
/// AssociationResult
///     Ranked candidates with posterior probabilities.
#[pyfunction]
#[pyo3(signature = (transient, candidates, config=None))]
fn associate_host(
    transient: &PyTransient,
    candidates: Vec<PyGalaxyCandidate>,
    config: Option<&PyAssociationConfig>,
) -> PyResult<PyAssociationResult> {
    let rs_transient = transient.to_rust();
    let rs_candidates: Vec<RsGalaxy> = candidates.iter().map(|c| c.to_rust()).collect();
    let rs_config = config
        .map(|c| c.to_rust())
        .unwrap_or_default();

    let result = rs_associate_host(&rs_transient, &rs_candidates, &rs_config)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    Ok(PyAssociationResult {
        candidates: result.candidates.into_iter().map(PyHostCandidate::from).collect(),
        p_none: result.p_none,
        n_considered: result.n_considered,
    })
}

/// Compute the directional light radius for a transient–galaxy pair.
///
/// Parameters
/// ----------
/// transient_ra, transient_dec : float
///     Transient position in degrees.
/// galaxy_ra, galaxy_dec : float
///     Galaxy center in degrees.
/// a_arcsec, b_arcsec : float
///     Galaxy semi-major and semi-minor axes in arcsec.
/// pa_deg : float
///     Galaxy position angle in degrees (N through E).
///
/// Returns
/// -------
/// dict
///     Keys: separation_arcsec, directional_radius, fractional_offset
#[pyfunction]
#[pyo3(signature = (transient_ra, transient_dec, galaxy_ra, galaxy_dec, a_arcsec, b_arcsec, pa_deg))]
fn compute_dlr(
    py: Python<'_>,
    transient_ra: f64,
    transient_dec: f64,
    galaxy_ra: f64,
    galaxy_dec: f64,
    a_arcsec: f64,
    b_arcsec: f64,
    pa_deg: f64,
) -> PyResult<Bound<'_, PyDict>> {
    let ellipse = RsEllipse::new(a_arcsec, b_arcsec, pa_deg)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let result = rs_compute_dlr(transient_ra, transient_dec, galaxy_ra, galaxy_dec, &ellipse);

    let dict = PyDict::new(py);
    dict.set_item("separation_arcsec", result.separation_arcsec)?;
    dict.set_item("directional_radius", result.directional_radius)?;
    dict.set_item("fractional_offset", result.fractional_offset)?;
    Ok(dict)
}

/// Compute offset likelihood using Gamma(0.75) distribution.
///
/// Parameters
/// ----------
/// fractional_offset : float or numpy array
///     Fractional offset value(s).
///
/// Returns
/// -------
/// float or numpy array
///     Likelihood value(s).
#[pyfunction]
fn offset_likelihood(py: Python<'_>, fractional_offset: PyObject) -> PyResult<PyObject> {
    // Try as array first, then scalar
    if let Ok(arr) = fractional_offset.extract::<PyReadonlyArray1<f64>>(py) {
        let slice = arr.as_slice().unwrap();
        let result: Vec<f64> = slice.iter().map(|&x| rs_offset_likelihood(x)).collect();
        Ok(PyArray1::from_vec(py, result).into_any().unbind())
    } else {
        let x: f64 = fractional_offset.extract(py)?;
        Ok(rs_offset_likelihood(x).into_pyobject(py)?.into_any().unbind())
    }
}

/// Compute redshift likelihood for a galaxy–transient pair.
#[pyfunction]
#[pyo3(signature = (galaxy_z=None, galaxy_z_err=None, transient_z=None, transient_z_err=None))]
fn redshift_likelihood(
    galaxy_z: Option<f64>,
    galaxy_z_err: Option<f64>,
    transient_z: Option<f64>,
    transient_z_err: Option<f64>,
) -> f64 {
    rs_redshift_likelihood(galaxy_z, galaxy_z_err, transient_z, transient_z_err)
}

/// Evaluate a Sérsic surface brightness profile.
///
/// Parameters
/// ----------
/// r : float
///     Radius from galaxy center (arcsec or pixels).
/// r_eff : float
///     Effective (half-light) radius (same units as r).
/// n : float
///     Sérsic index (1=exponential, 4=de Vaucouleurs).
/// i_eff : float
///     Surface brightness at r_eff.
///
/// Returns
/// -------
/// float
///     Surface brightness at radius r.
#[pyfunction]
fn sersic_profile(r: f64, r_eff: f64, n: f64, i_eff: f64) -> f64 {
    rs_sersic_profile(r, r_eff, n, i_eff)
}

/// Compute elliptical radius given offsets and shape parameters.
#[pyfunction]
fn elliptical_radius(dx: f64, dy: f64, axis_ratio: f64, pa_deg: f64) -> f64 {
    rs_elliptical_radius(dx, dy, axis_ratio, pa_deg.to_radians())
}

// ---------------------------------------------------------------------------
// Cutout
// ---------------------------------------------------------------------------

#[pyclass(name = "Cutout")]
#[derive(Clone)]
pub struct PyCutout {
    inner: RsCutout,
}

#[pymethods]
impl PyCutout {
    /// Create a new cutout from pixel data.
    ///
    /// Parameters
    /// ----------
    /// data : numpy.ndarray
    ///     2D array of pixel values (row-major).
    /// pixel_scale : float
    ///     Pixel scale in arcsec/pixel.
    /// center_ra : float
    ///     RA of image center in degrees.
    /// center_dec : float
    ///     Dec of image center in degrees.
    /// survey : str, optional
    ///     Survey name: "legacy", "panstarrs", or "skymapper". Default: "legacy".
    #[new]
    #[pyo3(signature = (data, pixel_scale, center_ra, center_dec, survey="legacy"))]
    fn new(
        data: PyReadonlyArray2<f64>,
        pixel_scale: f64,
        center_ra: f64,
        center_dec: f64,
        survey: &str,
    ) -> PyResult<Self> {
        let shape = data.shape();
        let height = shape[0];
        let width = shape[1];
        let flat: Vec<f64> = data.as_slice()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?
            .to_vec();
        let surv = match survey.to_lowercase().as_str() {
            "panstarrs" | "ps1" => ImagingSurvey::PanSTARRS,
            "skymapper" | "sm" => ImagingSurvey::SkyMapper,
            _ => ImagingSurvey::LegacySurvey,
        };
        let cutout = RsCutout::new(flat, width, height, pixel_scale, center_ra, center_dec, surv)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner: cutout })
    }

    /// Create a cutout filled with zeros.
    #[staticmethod]
    fn zeros(width: usize, height: usize, pixel_scale: f64, center_ra: f64, center_dec: f64) -> Self {
        Self {
            inner: RsCutout::zeros(width, height, pixel_scale, center_ra, center_dec),
        }
    }

    #[getter]
    fn width(&self) -> usize {
        self.inner.width
    }

    #[getter]
    fn height(&self) -> usize {
        self.inner.height
    }

    #[getter]
    fn pixel_scale(&self) -> f64 {
        self.inner.pixel_scale
    }

    #[getter]
    fn center_ra(&self) -> f64 {
        self.inner.center_ra
    }

    #[getter]
    fn center_dec(&self) -> f64 {
        self.inner.center_dec
    }

    /// Get pixel data as a 2D numpy array.
    #[getter]
    fn data<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        PyArray2::from_vec2(
            py,
            &self.inner.data
                .chunks(self.inner.width)
                .map(|row| row.to_vec())
                .collect::<Vec<_>>(),
        )
        .unwrap()
    }

    /// Get pixel value at (row, col).
    fn get(&self, row: usize, col: usize) -> Option<f64> {
        self.inner.get(row, col)
    }

    /// Set pixel value at (row, col).
    fn set(&mut self, row: usize, col: usize, value: f64) -> bool {
        self.inner.set(row, col, value)
    }

    /// Convert pixel coordinates to sky coordinates (RA, Dec).
    fn pixel_to_sky(&self, row: f64, col: f64) -> (f64, f64) {
        self.inner.pixel_to_sky(row, col)
    }

    /// Convert sky coordinates to pixel coordinates (row, col).
    fn sky_to_pixel(&self, ra: f64, dec: f64) -> (f64, f64) {
        self.inner.sky_to_pixel(ra, dec)
    }

    /// Size of cutout in arcsec (width, height).
    fn size_arcsec(&self) -> (f64, f64) {
        self.inner.size_arcsec()
    }

    /// Extract a sub-image centered on a pixel position.
    fn sub_image(&self, center_row: usize, center_col: usize, radius: usize) -> Option<PyCutout> {
        self.inner
            .sub_image(center_row, center_col, radius)
            .map(|c| PyCutout { inner: c })
    }

    /// Add a 2D Gaussian to the cutout.
    fn add_gaussian(
        &mut self,
        center_row: f64,
        center_col: f64,
        amplitude: f64,
        sigma_major: f64,
        sigma_minor: f64,
        pa_rad: f64,
    ) {
        self.inner.add_gaussian(center_row, center_col, amplitude, sigma_major, sigma_minor, pa_rad);
    }

    /// Add a 2D Sérsic profile to the cutout.
    #[pyo3(signature = (center_row, center_col, i_eff, r_eff_pix, n, axis_ratio, pa_rad))]
    fn add_sersic(
        &mut self,
        center_row: f64,
        center_col: f64,
        i_eff: f64,
        r_eff_pix: f64,
        n: f64,
        axis_ratio: f64,
        pa_rad: f64,
    ) {
        self.inner.add_sersic(center_row, center_col, i_eff, r_eff_pix, n, axis_ratio, pa_rad);
    }

    fn __repr__(&self) -> String {
        format!(
            "Cutout({}x{}, scale={:.3}\"/px, center=({:.4}, {:.4}))",
            self.inner.width, self.inner.height, self.inner.pixel_scale,
            self.inner.center_ra, self.inner.center_dec
        )
    }
}

// ---------------------------------------------------------------------------
// ExtractionConfig
// ---------------------------------------------------------------------------

#[pyclass(name = "ExtractionConfig")]
#[derive(Clone)]
pub struct PyExtractionConfig {
    #[pyo3(get, set)]
    pub thresh_sigma: f64,
    #[pyo3(get, set)]
    pub min_pixels: usize,
    #[pyo3(get, set)]
    pub back_size: usize,
    #[pyo3(get, set)]
    pub back_clip_iters: usize,
    #[pyo3(get, set)]
    pub back_clip_sigma: f64,
    #[pyo3(get, set)]
    pub deblend_contrast: f64,
}

#[pymethods]
impl PyExtractionConfig {
    #[new]
    #[pyo3(signature = (thresh_sigma=1.5, min_pixels=5, back_size=64, back_clip_iters=3, back_clip_sigma=3.0, deblend_contrast=0.005))]
    fn new(
        thresh_sigma: f64,
        min_pixels: usize,
        back_size: usize,
        back_clip_iters: usize,
        back_clip_sigma: f64,
        deblend_contrast: f64,
    ) -> Self {
        Self {
            thresh_sigma,
            min_pixels,
            back_size,
            back_clip_iters,
            back_clip_sigma,
            deblend_contrast,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "ExtractionConfig(thresh={:.1}σ, min_pix={}, back_size={})",
            self.thresh_sigma, self.min_pixels, self.back_size
        )
    }
}

impl PyExtractionConfig {
    fn to_rust(&self) -> RsExtractionConfig {
        RsExtractionConfig {
            thresh_sigma: self.thresh_sigma,
            min_pixels: self.min_pixels,
            back_size: self.back_size,
            back_clip_iters: self.back_clip_iters,
            back_clip_sigma: self.back_clip_sigma,
            deblend_contrast: self.deblend_contrast,
        }
    }
}

// ---------------------------------------------------------------------------
// DetectedSource
// ---------------------------------------------------------------------------

#[pyclass(name = "DetectedSource")]
#[derive(Clone)]
pub struct PyDetectedSource {
    #[pyo3(get)]
    pub x: f64,
    #[pyo3(get)]
    pub y: f64,
    #[pyo3(get)]
    pub ra: f64,
    #[pyo3(get)]
    pub dec: f64,
    #[pyo3(get)]
    pub a_pix: f64,
    #[pyo3(get)]
    pub b_pix: f64,
    #[pyo3(get)]
    pub pa_deg: f64,
    #[pyo3(get)]
    pub flux: f64,
    #[pyo3(get)]
    pub npix: usize,
    #[pyo3(get)]
    pub peak: f64,
    #[pyo3(get)]
    pub snr: f64,
}

#[pymethods]
impl PyDetectedSource {
    /// Semi-major axis in arcsec.
    fn a_arcsec(&self, pixel_scale: f64) -> f64 {
        self.a_pix * pixel_scale
    }

    /// Semi-minor axis in arcsec.
    fn b_arcsec(&self, pixel_scale: f64) -> f64 {
        self.b_pix * pixel_scale
    }

    /// Axis ratio b/a.
    fn axis_ratio(&self) -> f64 {
        if self.a_pix > 0.0 { self.b_pix / self.a_pix } else { 1.0 }
    }

    /// Ellipticity 1 - b/a.
    fn ellipticity(&self) -> f64 {
        1.0 - self.axis_ratio()
    }

    fn __repr__(&self) -> String {
        format!(
            "DetectedSource(x={:.1}, y={:.1}, flux={:.1}, snr={:.1}, a={:.1}px, b={:.1}px)",
            self.x, self.y, self.flux, self.snr, self.a_pix, self.b_pix
        )
    }
}

impl From<RsDetectedSource> for PyDetectedSource {
    fn from(s: RsDetectedSource) -> Self {
        Self {
            x: s.x,
            y: s.y,
            ra: s.ra,
            dec: s.dec,
            a_pix: s.a_pix,
            b_pix: s.b_pix,
            pa_deg: s.pa_deg,
            flux: s.flux,
            npix: s.npix,
            peak: s.peak,
            snr: s.snr,
        }
    }
}

// ---------------------------------------------------------------------------
// Background
// ---------------------------------------------------------------------------

#[pyclass(name = "Background")]
#[derive(Clone)]
pub struct PyBackground {
    #[pyo3(get)]
    pub global_mean: f64,
    #[pyo3(get)]
    pub global_rms: f64,
    #[pyo3(get)]
    pub width: usize,
    #[pyo3(get)]
    pub height: usize,
    level: Vec<f64>,
    rms: Vec<f64>,
}

#[pymethods]
impl PyBackground {
    /// Background level as 2D numpy array.
    #[getter]
    fn level<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        PyArray2::from_vec2(
            py,
            &self.level.chunks(self.width).map(|r| r.to_vec()).collect::<Vec<_>>(),
        )
        .unwrap()
    }

    /// Background RMS as 2D numpy array.
    #[getter]
    fn rms<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        PyArray2::from_vec2(
            py,
            &self.rms.chunks(self.width).map(|r| r.to_vec()).collect::<Vec<_>>(),
        )
        .unwrap()
    }

    fn __repr__(&self) -> String {
        format!(
            "Background({}x{}, mean={:.2}, rms={:.2})",
            self.width, self.height, self.global_mean, self.global_rms
        )
    }
}

impl From<RsBackground> for PyBackground {
    fn from(b: RsBackground) -> Self {
        Self {
            global_mean: b.global_mean,
            global_rms: b.global_rms,
            width: b.width,
            height: b.height,
            level: b.level,
            rms: b.rms,
        }
    }
}

// ---------------------------------------------------------------------------
// FitConfig
// ---------------------------------------------------------------------------

#[pyclass(name = "FitConfig")]
#[derive(Clone)]
pub struct PyFitConfig {
    #[pyo3(get, set)]
    pub max_sersic_n: f64,
    #[pyo3(get, set)]
    pub min_sersic_n: f64,
    #[pyo3(get, set)]
    pub min_r_eff_pix: f64,
    #[pyo3(get, set)]
    pub max_r_eff_pix: f64,
    #[pyo3(get, set)]
    pub max_iter: usize,
    #[pyo3(get, set)]
    pub tol: f64,
    #[pyo3(get, set)]
    pub aperture_factor: f64,
    #[pyo3(get, set)]
    pub lambda_init: f64,
}

#[pymethods]
impl PyFitConfig {
    #[new]
    #[pyo3(signature = (max_sersic_n=8.0, min_sersic_n=0.3, min_r_eff_pix=0.5, max_r_eff_pix=200.0, max_iter=200, tol=1e-6, aperture_factor=3.0, lambda_init=1e-3))]
    fn new(
        max_sersic_n: f64,
        min_sersic_n: f64,
        min_r_eff_pix: f64,
        max_r_eff_pix: f64,
        max_iter: usize,
        tol: f64,
        aperture_factor: f64,
        lambda_init: f64,
    ) -> Self {
        Self {
            max_sersic_n,
            min_sersic_n,
            min_r_eff_pix,
            max_r_eff_pix,
            max_iter,
            tol,
            aperture_factor,
            lambda_init,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "FitConfig(n=[{:.1}, {:.1}], r_eff=[{:.1}, {:.1}]px, max_iter={})",
            self.min_sersic_n, self.max_sersic_n, self.min_r_eff_pix, self.max_r_eff_pix, self.max_iter
        )
    }
}

impl PyFitConfig {
    fn to_rust(&self) -> RsFitConfig {
        RsFitConfig {
            max_sersic_n: self.max_sersic_n,
            min_sersic_n: self.min_sersic_n,
            min_r_eff_pix: self.min_r_eff_pix,
            max_r_eff_pix: self.max_r_eff_pix,
            max_iter: self.max_iter,
            tol: self.tol,
            aperture_factor: self.aperture_factor,
            lambda_init: self.lambda_init,
        }
    }
}

// ---------------------------------------------------------------------------
// SersicFit
// ---------------------------------------------------------------------------

#[pyclass(name = "SersicFit")]
#[derive(Clone)]
pub struct PySersicFit {
    #[pyo3(get)]
    pub ra: f64,
    #[pyo3(get)]
    pub dec: f64,
    #[pyo3(get)]
    pub x: f64,
    #[pyo3(get)]
    pub y: f64,
    #[pyo3(get)]
    pub r_eff: f64,
    #[pyo3(get)]
    pub r_eff_pix: f64,
    #[pyo3(get)]
    pub n: f64,
    #[pyo3(get)]
    pub axis_ratio: f64,
    #[pyo3(get)]
    pub pa_deg: f64,
    #[pyo3(get)]
    pub i_eff: f64,
    #[pyo3(get)]
    pub flux: f64,
    #[pyo3(get)]
    pub chi2_reduced: f64,
    #[pyo3(get)]
    pub ndata: usize,
    #[pyo3(get)]
    pub converged: bool,
}

#[pymethods]
impl PySersicFit {
    fn __repr__(&self) -> String {
        format!(
            "SersicFit(n={:.2}, r_eff={:.2}\", q={:.2}, PA={:.1}°, chi2r={:.2}, conv={})",
            self.n, self.r_eff, self.axis_ratio, self.pa_deg, self.chi2_reduced, self.converged
        )
    }
}

impl From<RsSersicFit> for PySersicFit {
    fn from(f: RsSersicFit) -> Self {
        Self {
            ra: f.ra,
            dec: f.dec,
            x: f.x,
            y: f.y,
            r_eff: f.r_eff,
            r_eff_pix: f.r_eff_pix,
            n: f.n,
            axis_ratio: f.axis_ratio,
            pa_deg: f.pa_deg,
            i_eff: f.i_eff,
            flux: f.flux,
            chi2_reduced: f.chi2_reduced,
            ndata: f.ndata,
            converged: f.converged,
        }
    }
}

// ---------------------------------------------------------------------------
// Image-based functions
// ---------------------------------------------------------------------------

/// Estimate the background of an image cutout.
///
/// Parameters
/// ----------
/// cutout : Cutout
///     The image cutout.
/// config : ExtractionConfig, optional
///     Configuration for background estimation.
///
/// Returns
/// -------
/// Background
///     Background level and RMS maps.
#[pyfunction]
#[pyo3(signature = (cutout, config=None))]
fn estimate_background(
    cutout: &PyCutout,
    config: Option<&PyExtractionConfig>,
) -> PyResult<PyBackground> {
    let rs_config = config
        .map(|c| c.to_rust())
        .unwrap_or_default();
    let bg = rs_estimate_background(&cutout.inner, &rs_config)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyBackground::from(bg))
}

/// Extract sources from an image cutout.
///
/// Parameters
/// ----------
/// cutout : Cutout
///     The image cutout to search.
/// config : ExtractionConfig, optional
///     Detection configuration. Uses defaults if not provided.
///
/// Returns
/// -------
/// list[DetectedSource]
///     Detected sources sorted by flux (brightest first).
#[pyfunction]
#[pyo3(signature = (cutout, config=None))]
fn extract_sources(
    cutout: &PyCutout,
    config: Option<&PyExtractionConfig>,
) -> PyResult<Vec<PyDetectedSource>> {
    let rs_config = config
        .map(|c| c.to_rust())
        .unwrap_or_default();
    let sources = rs_extract_sources(&cutout.inner, &rs_config)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(sources.into_iter().map(PyDetectedSource::from).collect())
}

/// Fit a Sérsic profile to a detected source.
///
/// Parameters
/// ----------
/// cutout : Cutout
///     The image cutout containing the source.
/// source : DetectedSource
///     The source to fit (from extract_sources).
/// config : FitConfig, optional
///     Fitting configuration. Uses defaults if not provided.
///
/// Returns
/// -------
/// SersicFit
///     Fitted Sérsic profile parameters.
#[pyfunction]
#[pyo3(signature = (cutout, source, config=None))]
fn fit_sersic(
    cutout: &PyCutout,
    source: &PyDetectedSource,
    config: Option<&PyFitConfig>,
) -> PyResult<PySersicFit> {
    let rs_config = config
        .map(|c| c.to_rust())
        .unwrap_or_default();
    let rs_source = RsDetectedSource {
        x: source.x,
        y: source.y,
        ra: source.ra,
        dec: source.dec,
        a_pix: source.a_pix,
        b_pix: source.b_pix,
        pa_deg: source.pa_deg,
        flux: source.flux,
        npix: source.npix,
        peak: source.peak,
        snr: source.snr,
    };
    let fit = rs_fit_sersic(&cutout.inner, &rs_source, &rs_config)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PySersicFit::from(fit))
}

// ---------------------------------------------------------------------------
// BOOM client
// ---------------------------------------------------------------------------

use prost_rs::boom::{BoomClient as RsBoomClient, BoomConfig as RsBoomConfig};

#[pyclass(name = "BoomConfig")]
#[derive(Clone)]
pub struct PyBoomConfig {
    #[pyo3(get, set)]
    pub base_url: String,
    #[pyo3(get, set)]
    pub username: Option<String>,
    #[pyo3(get, set)]
    pub password: Option<String>,
    #[pyo3(get, set)]
    pub default_radius_arcsec: f64,
    #[pyo3(get, set)]
    pub default_catalog: String,
    #[pyo3(get, set)]
    pub min_b_arcsec: f64,
}

#[pymethods]
impl PyBoomConfig {
    #[new]
    #[pyo3(signature = (base_url=None, username=None, password=None, default_radius_arcsec=60.0, default_catalog="LS_DR10", min_b_arcsec=0.05))]
    fn new(
        base_url: Option<String>,
        username: Option<String>,
        password: Option<String>,
        default_radius_arcsec: f64,
        default_catalog: &str,
        min_b_arcsec: f64,
    ) -> Self {
        let default = RsBoomConfig::default();
        Self {
            base_url: base_url.unwrap_or(default.base_url),
            username: username.or(default.username),
            password: password.or(default.password),
            default_radius_arcsec,
            default_catalog: default_catalog.to_string(),
            min_b_arcsec,
        }
    }

    /// Create a config from environment variables (BOOM_USERNAME, BOOM_PASSWORD).
    #[staticmethod]
    fn from_env() -> Self {
        let cfg = RsBoomConfig::from_env();
        Self {
            base_url: cfg.base_url,
            username: cfg.username,
            password: cfg.password,
            default_radius_arcsec: cfg.default_radius_arcsec,
            default_catalog: cfg.default_catalog,
            min_b_arcsec: cfg.min_b_arcsec,
        }
    }

    /// Create a config with explicit credentials.
    #[staticmethod]
    fn with_credentials(username: &str, password: &str) -> Self {
        let cfg = RsBoomConfig::with_credentials(username, password);
        Self {
            base_url: cfg.base_url,
            username: cfg.username,
            password: cfg.password,
            default_radius_arcsec: cfg.default_radius_arcsec,
            default_catalog: cfg.default_catalog,
            min_b_arcsec: cfg.min_b_arcsec,
        }
    }

    /// Check whether credentials are available.
    fn has_credentials(&self) -> bool {
        self.username.is_some() && self.password.is_some()
    }

    fn __repr__(&self) -> String {
        format!(
            "BoomConfig(url={:?}, has_creds={}, radius={:.0}\")",
            self.base_url, self.has_credentials(), self.default_radius_arcsec
        )
    }
}

impl PyBoomConfig {
    fn to_rust(&self) -> RsBoomConfig {
        RsBoomConfig {
            base_url: self.base_url.clone(),
            username: self.username.clone(),
            password: self.password.clone(),
            default_radius_arcsec: self.default_radius_arcsec,
            default_catalog: self.default_catalog.clone(),
            min_b_arcsec: self.min_b_arcsec,
        }
    }
}

#[pyclass(name = "BoomClient")]
pub struct PyBoomClient {
    client: RsBoomClient,
    rt: tokio::runtime::Runtime,
}

#[pymethods]
impl PyBoomClient {
    /// Create a new BOOM client.
    ///
    /// Parameters
    /// ----------
    /// config : BoomConfig, optional
    ///     Client configuration. Uses env-based defaults if not provided.
    #[new]
    #[pyo3(signature = (config=None))]
    fn new(config: Option<&PyBoomConfig>) -> PyResult<Self> {
        let rs_config = config
            .map(|c| c.to_rust())
            .unwrap_or_default();
        let rt = tokio::runtime::Runtime::new()
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        Ok(Self {
            client: RsBoomClient::new(rs_config),
            rt,
        })
    }

    /// Create a client with credentials from environment variables.
    #[staticmethod]
    fn from_env() -> PyResult<Self> {
        let rt = tokio::runtime::Runtime::new()
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        Ok(Self {
            client: RsBoomClient::from_env(),
            rt,
        })
    }

    /// Authenticate with the BOOM API.
    fn authenticate(&mut self) -> PyResult<()> {
        self.rt
            .block_on(self.client.authenticate())
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    /// Check whether the client has authenticated.
    fn is_authenticated(&self) -> bool {
        self.client.is_authenticated()
    }

    /// Get ZTF coordinates for an object by its objectId.
    ///
    /// Returns
    /// -------
    /// tuple[float, float]
    ///     (ra, dec) in degrees.
    fn get_ztf_coordinates(&self, object_id: &str) -> PyResult<(f64, f64)> {
        self.rt
            .block_on(self.client.get_ztf_coordinates(object_id))
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    /// Cone search Legacy Survey (LS_DR10) for galaxies near a position.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Object name (used as key in BOOM response).
    /// ra, dec : float
    ///     Search center in degrees.
    /// radius_arcsec : float, optional
    ///     Search radius. Uses config default if not provided.
    ///
    /// Returns
    /// -------
    /// list[dict]
    ///     Raw catalog documents from BOOM.
    #[pyo3(signature = (name, ra, dec, radius_arcsec=None))]
    fn cone_search_ls(
        &self,
        py: Python<'_>,
        name: &str,
        ra: f64,
        dec: f64,
        radius_arcsec: Option<f64>,
    ) -> PyResult<PyObject> {
        let docs = self.rt
            .block_on(self.client.cone_search_ls(name, ra, dec, radius_arcsec))
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        json_values_to_py(py, &docs)
    }

    /// Cone search NED for galaxy redshifts and types.
    #[pyo3(signature = (name, ra, dec, radius_arcsec=None))]
    fn cone_search_ned(
        &self,
        py: Python<'_>,
        name: &str,
        ra: f64,
        dec: f64,
        radius_arcsec: Option<f64>,
    ) -> PyResult<PyObject> {
        let docs = self.rt
            .block_on(self.client.cone_search_ned(name, ra, dec, radius_arcsec))
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        json_values_to_py(py, &docs)
    }

    /// Query BOOM for a ZTF transient and return galaxy candidates.
    ///
    /// High-level method: looks up coordinates, runs LS cone search,
    /// and converts results to GalaxyCandidate objects.
    ///
    /// Parameters
    /// ----------
    /// object_id : str
    ///     ZTF object ID (e.g., "ZTF20aajnksq").
    /// radius_arcsec : float, optional
    ///     Search radius. Uses config default if not provided.
    ///
    /// Returns
    /// -------
    /// tuple[tuple[float, float], list[GalaxyCandidate]]
    ///     ((ra, dec), candidates) — transient coordinates and galaxy candidates.
    #[pyo3(signature = (object_id, radius_arcsec=None))]
    fn get_candidates(
        &self,
        object_id: &str,
        radius_arcsec: Option<f64>,
    ) -> PyResult<((f64, f64), Vec<PyGalaxyCandidate>)> {
        let ((ra, dec), candidates) = self.rt
            .block_on(self.client.get_candidates(object_id, radius_arcsec))
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let py_candidates: Vec<PyGalaxyCandidate> = candidates.into_iter().map(|g| {
            PyGalaxyCandidate {
                ra: g.ra,
                dec: g.dec,
                a_arcsec: g.a_arcsec,
                b_arcsec: g.b_arcsec,
                pa_deg: g.pa_deg,
                redshift: g.redshift,
                redshift_err: g.redshift_err,
                mag: g.mag,
                mag_err: g.mag_err,
                objtype: g.objtype,
                objname: g.objname,
                catalog: g.catalog,
                shape_from_image: g.shape_from_image,
            }
        }).collect();

        Ok(((ra, dec), py_candidates))
    }

    fn __repr__(&self) -> String {
        format!(
            "BoomClient(authenticated={})",
            self.client.is_authenticated()
        )
    }
}

/// Convert serde_json::Value array to Python list of dicts.
fn json_values_to_py(py: Python<'_>, values: &[serde_json::Value]) -> PyResult<PyObject> {
    use pyo3::types::PyList;

    let items: Vec<PyObject> = values
        .iter()
        .map(|v| json_value_to_py(py, v))
        .collect::<PyResult<_>>()?;
    Ok(PyList::new(py, &items)?.into_any().unbind())
}

fn json_value_to_py(py: Python<'_>, value: &serde_json::Value) -> PyResult<PyObject> {
    use pyo3::types::{PyDict as PyDictType, PyList};

    match value {
        serde_json::Value::Null => Ok(py.None()),
        serde_json::Value::Bool(b) => Ok(b.into_pyobject(py).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?.to_owned().into_any().unbind()),
        serde_json::Value::Number(n) => {
            if let Some(i) = n.as_i64() {
                Ok(i.into_pyobject(py)?.into_any().unbind())
            } else if let Some(f) = n.as_f64() {
                Ok(f.into_pyobject(py)?.into_any().unbind())
            } else {
                Ok(py.None())
            }
        }
        serde_json::Value::String(s) => Ok(s.into_pyobject(py)?.into_any().unbind()),
        serde_json::Value::Array(arr) => {
            let items: Vec<PyObject> = arr
                .iter()
                .map(|v| json_value_to_py(py, v))
                .collect::<PyResult<_>>()?;
            Ok(PyList::new(py, &items)?.into_any().unbind())
        }
        serde_json::Value::Object(map) => {
            let dict = PyDictType::new(py);
            for (k, v) in map {
                dict.set_item(k, json_value_to_py(py, v)?)?;
            }
            Ok(dict.into_any().unbind())
        }
    }
}

// ---------------------------------------------------------------------------
// Module
// ---------------------------------------------------------------------------

#[pymodule]
fn prost_extension(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Catalog-based association
    m.add_class::<PyTransient>()?;
    m.add_class::<PyGalaxyCandidate>()?;
    m.add_class::<PyHostCandidate>()?;
    m.add_class::<PyAssociationConfig>()?;
    m.add_class::<PyAssociationResult>()?;
    m.add_function(wrap_pyfunction!(associate_host, m)?)?;
    m.add_function(wrap_pyfunction!(compute_dlr, m)?)?;
    m.add_function(wrap_pyfunction!(offset_likelihood, m)?)?;
    m.add_function(wrap_pyfunction!(redshift_likelihood, m)?)?;
    // Image-based analysis
    m.add_class::<PyCutout>()?;
    m.add_class::<PyExtractionConfig>()?;
    m.add_class::<PyDetectedSource>()?;
    m.add_class::<PyBackground>()?;
    m.add_class::<PyFitConfig>()?;
    m.add_class::<PySersicFit>()?;
    m.add_function(wrap_pyfunction!(estimate_background, m)?)?;
    m.add_function(wrap_pyfunction!(extract_sources, m)?)?;
    m.add_function(wrap_pyfunction!(fit_sersic, m)?)?;
    m.add_function(wrap_pyfunction!(sersic_profile, m)?)?;
    m.add_function(wrap_pyfunction!(elliptical_radius, m)?)?;
    // BOOM client
    m.add_class::<PyBoomConfig>()?;
    m.add_class::<PyBoomClient>()?;
    Ok(())
}
