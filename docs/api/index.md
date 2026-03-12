# API Reference

## Catalog-Based Association

### `Transient`

A transient event to associate with a host galaxy.

```python
Transient(ra, dec, ra_err=0.0, dec_err=0.0, redshift=None, redshift_err=None, name=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `ra` | `float` | Right ascension (degrees) |
| `dec` | `float` | Declination (degrees) |
| `ra_err` | `float` | RA uncertainty (arcsec) |
| `dec_err` | `float` | Dec uncertainty (arcsec) |
| `redshift` | `float \| None` | Spectroscopic or photometric redshift |
| `redshift_err` | `float \| None` | Redshift uncertainty |
| `name` | `str \| None` | Identifier |

=== "Rust"

    ```rust
    let t = Transient::new(197.45, -23.38)
        .with_position_err(0.1, 0.1)
        .with_redshift(0.00978, 0.00015);
    ```

---

### `GalaxyCandidate`

A galaxy from a catalog or image analysis.

```python
GalaxyCandidate(ra, dec, a_arcsec, b_arcsec, pa_deg, redshift=None,
                redshift_err=None, mag=None, mag_err=None, objtype=None,
                objname=None, catalog=None, shape_from_image=False)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `ra` | `float` | Right ascension (degrees) |
| `dec` | `float` | Declination (degrees) |
| `a_arcsec` | `float` | Semi-major axis (arcsec) |
| `b_arcsec` | `float` | Semi-minor axis (arcsec) |
| `pa_deg` | `float` | Position angle (degrees, N through E) |
| `redshift` | `float \| None` | Redshift |
| `redshift_err` | `float \| None` | Redshift uncertainty |
| `mag` | `float \| None` | Apparent magnitude |
| `mag_err` | `float \| None` | Magnitude uncertainty |
| `objtype` | `str \| None` | Morphological type (e.g., `"EXP"`, `"DEV"`, `"SER"`) |
| `objname` | `str \| None` | Object identifier |
| `catalog` | `str \| None` | Source catalog (e.g., `"LS_DR10"`) |
| `shape_from_image` | `bool` | Whether shape was derived from image analysis (default `False`) |

#### Static Methods

##### `GalaxyCandidate.from_tractor(ra, dec, shape_r, shape_e1, shape_e2, min_b=0.05, redshift=None, objtype=None, objname=None)`

Create a `GalaxyCandidate` from Legacy Survey Tractor shape parameters.
Converts Tractor ellipticity components to semi-axes and position angle
(see [Tractor Shape Parameters](../physics/association.md#tractor-shape-parameters)).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `shape_r` | `float` | required | Tractor effective radius (arcsec) |
| `shape_e1` | `float` | required | Tractor ellipticity component 1 |
| `shape_e2` | `float` | required | Tractor ellipticity component 2 |
| `min_b` | `float` | 0.05 | Floor for semi-minor axis (arcsec) |

Returns a `GalaxyCandidate` with `a_arcsec`, `b_arcsec`, and `pa_deg` derived
from the Tractor parameters. Raises `ValueError` if `shape_r <= 0`.

---

### `AssociationConfig`

Configuration for the association algorithm.

```python
AssociationConfig(max_fractional_offset=10.0, min_b_arcsec=0.05,
                  max_candidates=10, use_redshift=True, use_absmag=False)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_fractional_offset` | `float` | 10.0 | Maximum fractional offset to consider a galaxy |
| `min_b_arcsec` | `float` | 0.05 | Floor for semi-minor axis (prevents division by zero) |
| `max_candidates` | `int` | 10 | Maximum number of ranked candidates returned |
| `use_redshift` | `bool` | True | Include redshift likelihood in posterior |
| `use_absmag` | `bool` | False | Include absolute magnitude likelihood in posterior |

---

### `HostCandidate`

Result for a single galaxy candidate (read-only).

| Property | Type | Description |
|----------|------|-------------|
| `galaxy` | `GalaxyCandidate` | The galaxy's properties |
| `separation_arcsec` | `float` | Angular separation from transient (arcsec) |
| `dlr` | `float` | Directional light radius at the transient's angle (arcsec) |
| `fractional_offset` | `float` | `separation_arcsec / dlr` |
| `dlr_rank` | `int` | Rank by fractional offset (1 = closest in DLR) |
| `posterior` | `float` | Normalized posterior probability P(host) |
| `posterior_offset` | `float` | Offset likelihood component |
| `posterior_redshift` | `float` | Redshift likelihood component (1.0 if no redshift) |
| `posterior_absmag` | `float` | Absolute magnitude component (1.0 if unused) |

---

### `AssociationResult`

Full association result.

| Property | Type | Description |
|----------|------|-------------|
| `candidates` | `list[HostCandidate]` | Ranked candidates (best first by posterior) |
| `best_host` | `HostCandidate \| None` | First candidate, or `None` if empty |
| `p_none` | `float` | P(no host among candidates) |
| `n_considered` | `int` | Number of input galaxies evaluated |

`len(result)` returns the number of candidates.

---

## Association Functions

### `associate_host(transient, candidates, config=None)`

Main entry point. Ranks galaxy candidates by DLR, computes Bayesian
posteriors, and returns an `AssociationResult`.

| Parameter | Type | Description |
|-----------|------|-------------|
| `transient` | `Transient` | The transient event |
| `candidates` | `list[GalaxyCandidate]` | Galaxy candidates to evaluate |
| `config` | `AssociationConfig \| None` | Configuration (uses defaults if `None`) |

Returns `AssociationResult`. Raises `ValueError` on invalid input.

### `compute_dlr(transient_ra, transient_dec, galaxy_ra, galaxy_dec, a_arcsec, b_arcsec, pa_deg)`

Compute DLR for a single transient-galaxy pair. Returns a `dict`:

| Key | Type | Description |
|-----|------|-------------|
| `separation_arcsec` | `float` | Angular separation |
| `directional_radius` | `float` | Galaxy radius in direction of transient |
| `fractional_offset` | `float` | `separation / directional_radius` |

### `offset_likelihood(fractional_offset)`

Evaluate the Gamma(0.75) offset likelihood. Accepts a scalar `float` or
a `numpy.ndarray` and returns the same type.

### `redshift_likelihood(galaxy_z=None, galaxy_z_err=None, transient_z=None, transient_z_err=None)`

Gaussian redshift likelihood for a galaxy-transient pair. Returns `1.0`
if either redshift is `None` (agnostic).

---

## Image-Based Analysis

### `Cutout`

A 2D image cutout with pixel/sky coordinate transforms.

```python
Cutout(data, pixel_scale, center_ra, center_dec, survey="legacy")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `data` | `numpy.ndarray` | 2D array of pixel values (float64, row-major) |
| `pixel_scale` | `float` | Pixel scale (arcsec/pixel) |
| `center_ra` | `float` | RA of image center (degrees) |
| `center_dec` | `float` | Dec of image center (degrees) |
| `survey` | `str` | Survey name: `"legacy"`, `"panstarrs"`/`"ps1"`, or `"skymapper"`/`"sm"` |

=== "Rust"

    ```rust
    let cutout = Cutout::new(
        data, width, height, pixel_scale,
        center_ra, center_dec, ImagingSurvey::LegacySurvey,
    )?;
    ```

#### Properties

| Property | Type | Description |
|----------|------|-------------|
| `width` | `int` | Image width in pixels |
| `height` | `int` | Image height in pixels |
| `pixel_scale` | `float` | Pixel scale (arcsec/pixel) |
| `center_ra` | `float` | RA of image center (degrees) |
| `center_dec` | `float` | Dec of image center (degrees) |
| `data` | `numpy.ndarray` | 2D pixel array (read-only copy) |

#### Methods

##### `Cutout.zeros(width, height, pixel_scale, center_ra, center_dec)`

Static method. Create a cutout filled with zeros (useful for testing or
building synthetic images).

##### `cutout.get(row, col)`

Get pixel value at `(row, col)`. Returns `None` if out of bounds.

##### `cutout.set(row, col, value)`

Set pixel value at `(row, col)`. Returns `True` on success, `False` if
out of bounds.

##### `cutout.pixel_to_sky(row, col)`

Convert pixel coordinates to sky coordinates using a TAN projection.
Returns `(ra, dec)` in degrees.

##### `cutout.sky_to_pixel(ra, dec)`

Convert sky coordinates to pixel coordinates. Returns `(row, col)`.

##### `cutout.size_arcsec()`

Returns `(width_arcsec, height_arcsec)` — the angular size of the cutout.

##### `cutout.sub_image(center_row, center_col, radius)`

Extract a square sub-image of side `2*radius+1` centered on `(center_row,
center_col)`. Returns a new `Cutout`, or `None` if the region extends
beyond the image boundary.

##### `cutout.add_gaussian(center_row, center_col, amplitude, sigma_major, sigma_minor, pa_rad)`

Add a 2D elliptical Gaussian to the pixel data in-place. Useful for
generating synthetic test images.

##### `cutout.add_sersic(center_row, center_col, i_eff, r_eff_pix, n, axis_ratio, pa_rad)`

Add a 2D Sersic profile to the pixel data in-place.

| Parameter | Type | Description |
|-----------|------|-------------|
| `i_eff` | `float` | Surface brightness at the effective radius |
| `r_eff_pix` | `float` | Effective (half-light) radius in pixels |
| `n` | `float` | Sersic index (1 = exponential, 4 = de Vaucouleurs) |
| `axis_ratio` | `float` | Axis ratio b/a (0, 1] |
| `pa_rad` | `float` | Position angle in radians |

---

### `ExtractionConfig`

Configuration for source extraction and background estimation.

```python
ExtractionConfig(thresh_sigma=1.5, min_pixels=5, back_size=64,
                 back_clip_iters=3, back_clip_sigma=3.0,
                 deblend_contrast=0.005)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `thresh_sigma` | `float` | 1.5 | Detection threshold in sigma above background |
| `min_pixels` | `int` | 5 | Minimum connected pixels for a detection |
| `back_size` | `int` | 64 | Background mesh cell size (pixels) |
| `back_clip_iters` | `int` | 3 | Sigma-clipping iterations for background |
| `back_clip_sigma` | `float` | 3.0 | Sigma-clipping threshold |
| `deblend_contrast` | `float` | 0.005 | Deblend minimum contrast (0 = off) |

---

### `DetectedSource`

A source detected in an image cutout (read-only).

| Property | Type | Description |
|----------|------|-------------|
| `x` | `float` | Centroid column (pixels) |
| `y` | `float` | Centroid row (pixels) |
| `ra` | `float` | Centroid RA (degrees) |
| `dec` | `float` | Centroid Dec (degrees) |
| `a_pix` | `float` | Semi-major axis from second moments (pixels) |
| `b_pix` | `float` | Semi-minor axis from second moments (pixels) |
| `pa_deg` | `float` | Position angle (degrees, CCW from +x axis) |
| `flux` | `float` | Total flux above background (ADU) |
| `npix` | `int` | Number of pixels above threshold |
| `peak` | `float` | Peak pixel value (background-subtracted) |
| `snr` | `float` | Signal-to-noise ratio (flux / flux_err) |

#### Methods

##### `source.a_arcsec(pixel_scale)`

Semi-major axis converted to arcsec: `a_pix * pixel_scale`.

##### `source.b_arcsec(pixel_scale)`

Semi-minor axis converted to arcsec: `b_pix * pixel_scale`.

##### `source.axis_ratio()`

Axis ratio `b/a`. Returns `1.0` if `a_pix == 0`.

##### `source.ellipticity()`

Ellipticity `1 - b/a`.

---

### `Background`

Background estimation result from `estimate_background()`.

| Property | Type | Description |
|----------|------|-------------|
| `global_mean` | `float` | Mean background level |
| `global_rms` | `float` | Global background RMS |
| `width` | `int` | Image width |
| `height` | `int` | Image height |
| `level` | `numpy.ndarray` | 2D background level map (same size as image) |
| `rms` | `numpy.ndarray` | 2D background RMS map |

---

### `FitConfig`

Configuration for Sersic profile fitting.

```python
FitConfig(max_sersic_n=8.0, min_sersic_n=0.3, min_r_eff_pix=0.5,
          max_r_eff_pix=200.0, max_iter=200, tol=1e-6,
          aperture_factor=3.0, lambda_init=1e-3)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_sersic_n` | `float` | 8.0 | Maximum Sersic index allowed |
| `min_sersic_n` | `float` | 0.3 | Minimum Sersic index allowed |
| `min_r_eff_pix` | `float` | 0.5 | Minimum effective radius (pixels) |
| `max_r_eff_pix` | `float` | 200.0 | Maximum effective radius (pixels) |
| `max_iter` | `int` | 200 | Maximum Levenberg-Marquardt iterations |
| `tol` | `float` | 1e-6 | Convergence tolerance (relative chi-squared change) |
| `aperture_factor` | `float` | 3.0 | Fitting aperture = `aperture_factor * a_pix` |
| `lambda_init` | `float` | 1e-3 | Initial LM damping parameter |

---

### `SersicFit`

Result of fitting a Sersic profile to a detected source (read-only).

| Property | Type | Description |
|----------|------|-------------|
| `ra` | `float` | Fitted centroid RA (degrees) |
| `dec` | `float` | Fitted centroid Dec (degrees) |
| `x` | `float` | Fitted centroid column (pixels) |
| `y` | `float` | Fitted centroid row (pixels) |
| `r_eff` | `float` | Effective (half-light) radius (arcsec) |
| `r_eff_pix` | `float` | Effective radius (pixels) |
| `n` | `float` | Sersic index (1 = exponential, 4 = de Vaucouleurs) |
| `axis_ratio` | `float` | Axis ratio b/a |
| `pa_deg` | `float` | Position angle (degrees, N through E) |
| `i_eff` | `float` | Surface brightness at r_eff |
| `flux` | `float` | Total integrated flux |
| `chi2_reduced` | `float` | Reduced chi-squared of the fit |
| `ndata` | `int` | Number of pixels used in fit |
| `converged` | `bool` | Whether the optimizer converged |

---

## Image-Based Functions

### `estimate_background(cutout, config=None)`

Estimate the background and RMS of an image using sigma-clipped
statistics on a mesh grid, bilinearly interpolated to full resolution.

| Parameter | Type | Description |
|-----------|------|-------------|
| `cutout` | `Cutout` | The image cutout |
| `config` | `ExtractionConfig \| None` | Controls `back_size`, `back_clip_*` params |

Returns `Background`.

### `extract_sources(cutout, config=None)`

Detect sources in an image cutout using background subtraction,
thresholding, and connected-component labeling with 8-connectivity.
Source shapes are measured from flux-weighted second moments.

| Parameter | Type | Description |
|-----------|------|-------------|
| `cutout` | `Cutout` | The image cutout to search |
| `config` | `ExtractionConfig \| None` | Detection parameters |

Returns `list[DetectedSource]` sorted by flux (brightest first).

### `fit_sersic(cutout, source, config=None)`

Fit a 2D Sersic profile to a detected source using Levenberg-Marquardt
optimization with numerical Jacobian. Initial parameters are derived from
the source's moment-based shape measurements.

| Parameter | Type | Description |
|-----------|------|-------------|
| `cutout` | `Cutout` | The image containing the source |
| `source` | `DetectedSource` | Source to fit (from `extract_sources`) |
| `config` | `FitConfig \| None` | Fitting bounds and convergence criteria |

Returns `SersicFit`. Raises `RuntimeError` if the fit fails (e.g., too
few pixels in the aperture).

---

## Morphology Utilities

### `sersic_profile(r, r_eff, n, i_eff)`

Evaluate a Sersic surface brightness profile:

\[
I(r) = I_e \exp\left(-b_n \left[\left(\frac{r}{r_e}\right)^{1/n} - 1\right]\right)
\]

where \(b_n \approx 1.9992n - 0.3271\).

| Parameter | Type | Description |
|-----------|------|-------------|
| `r` | `float` | Radius from galaxy center |
| `r_eff` | `float` | Effective (half-light) radius (same units as `r`) |
| `n` | `float` | Sersic index |
| `i_eff` | `float` | Surface brightness at `r_eff` |

### `elliptical_radius(dx, dy, axis_ratio, pa_deg)`

Compute the elliptical radius given pixel offsets and shape parameters.
Equivalent to projecting `(dx, dy)` into the galaxy frame and computing
the distance along the ellipse.

| Parameter | Type | Description |
|-----------|------|-------------|
| `dx` | `float` | Offset in x (pixels or arcsec) |
| `dy` | `float` | Offset in y (same units as `dx`) |
| `axis_ratio` | `float` | Axis ratio b/a |
| `pa_deg` | `float` | Position angle (degrees) |
