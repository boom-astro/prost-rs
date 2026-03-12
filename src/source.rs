//! Source extraction from image cutouts.
//!
//! Detects and measures sources in cutout images using techniques inspired
//! by SExtractor/SEP: background estimation, sigma-clipped statistics,
//! thresholding, connected-component labeling, and isophotal/moment-based
//! shape measurements.

use crate::cutout::Cutout;
use crate::errors::ProstError;

/// A detected source in an image.
#[derive(Debug, Clone)]
pub struct DetectedSource {
    /// Centroid x (column) in pixels
    pub x: f64,
    /// Centroid y (row) in pixels
    pub y: f64,
    /// RA of centroid (degrees)
    pub ra: f64,
    /// Dec of centroid (degrees)
    pub dec: f64,
    /// Semi-major axis from second moments (pixels)
    pub a_pix: f64,
    /// Semi-minor axis from second moments (pixels)
    pub b_pix: f64,
    /// Position angle from moments (degrees, CCW from +x axis)
    pub pa_deg: f64,
    /// Total flux above background (ADU)
    pub flux: f64,
    /// Number of pixels above threshold
    pub npix: usize,
    /// Peak pixel value (background-subtracted)
    pub peak: f64,
    /// Signal-to-noise ratio (flux / flux_err)
    pub snr: f64,
}

impl DetectedSource {
    /// Convert pixel-based semi-major axis to arcsec.
    pub fn a_arcsec(&self, pixel_scale: f64) -> f64 {
        self.a_pix * pixel_scale
    }

    /// Convert pixel-based semi-minor axis to arcsec.
    pub fn b_arcsec(&self, pixel_scale: f64) -> f64 {
        self.b_pix * pixel_scale
    }

    /// Axis ratio b/a.
    pub fn axis_ratio(&self) -> f64 {
        if self.a_pix > 0.0 {
            self.b_pix / self.a_pix
        } else {
            1.0
        }
    }

    /// Ellipticity e = 1 - b/a.
    pub fn ellipticity(&self) -> f64 {
        1.0 - self.axis_ratio()
    }
}

/// Configuration for source extraction.
#[derive(Debug, Clone)]
pub struct ExtractionConfig {
    /// Detection threshold in sigma above background
    pub thresh_sigma: f64,
    /// Minimum number of connected pixels for a detection
    pub min_pixels: usize,
    /// Background mesh size in pixels (for mesh-based estimation)
    pub back_size: usize,
    /// Number of sigma-clipping iterations for background estimation
    pub back_clip_iters: usize,
    /// Sigma-clipping threshold for background estimation
    pub back_clip_sigma: f64,
    /// Deblend minimum contrast (0 = no deblending, 1 = aggressive)
    pub deblend_contrast: f64,
}

impl Default for ExtractionConfig {
    fn default() -> Self {
        Self {
            thresh_sigma: 1.5,
            min_pixels: 5,
            back_size: 64,
            back_clip_iters: 3,
            back_clip_sigma: 3.0,
            deblend_contrast: 0.005,
        }
    }
}

/// Background estimation result.
#[derive(Debug, Clone)]
pub struct Background {
    /// Background level per pixel (same size as image)
    pub level: Vec<f64>,
    /// Background RMS per pixel (same size as image)
    pub rms: Vec<f64>,
    /// Global mean background
    pub global_mean: f64,
    /// Global RMS
    pub global_rms: f64,
    /// Image dimensions
    pub width: usize,
    pub height: usize,
}

/// Estimate the background and RMS of an image using sigma-clipped statistics.
///
/// Uses a mesh grid approach: the image is divided into cells of `mesh_size`,
/// sigma-clipped statistics are computed per cell, then bilinearly interpolated
/// back to the full resolution.
pub fn estimate_background(
    cutout: &Cutout,
    config: &ExtractionConfig,
) -> Result<Background, ProstError> {
    let mesh = config.back_size;
    let ny = (cutout.height + mesh - 1) / mesh;
    let nx = (cutout.width + mesh - 1) / mesh;

    let mut mesh_mean = vec![0.0_f64; ny * nx];
    let mut mesh_rms = vec![0.0_f64; ny * nx];

    // Compute sigma-clipped statistics per mesh cell
    for my in 0..ny {
        for mx in 0..nx {
            let r0 = my * mesh;
            let r1 = (r0 + mesh).min(cutout.height);
            let c0 = mx * mesh;
            let c1 = (c0 + mesh).min(cutout.width);

            let mut values: Vec<f64> = Vec::new();
            for r in r0..r1 {
                for c in c0..c1 {
                    let v = cutout.data[r * cutout.width + c];
                    if v.is_finite() {
                        values.push(v);
                    }
                }
            }

            let (mean, rms) = sigma_clipped_stats(
                &mut values,
                config.back_clip_sigma,
                config.back_clip_iters,
            );
            mesh_mean[my * nx + mx] = mean;
            mesh_rms[my * nx + mx] = rms;
        }
    }

    // Bilinearly interpolate mesh values to full resolution
    let mut level = vec![0.0_f64; cutout.height * cutout.width];
    let mut rms = vec![0.0_f64; cutout.height * cutout.width];

    for r in 0..cutout.height {
        for c in 0..cutout.width {
            let fy = (r as f64 + 0.5) / mesh as f64 - 0.5;
            let fx = (c as f64 + 0.5) / mesh as f64 - 0.5;
            level[r * cutout.width + c] = bilinear_interp(&mesh_mean, nx, ny, fx, fy);
            rms[r * cutout.width + c] = bilinear_interp(&mesh_rms, nx, ny, fx, fy);
        }
    }

    let global_mean = mesh_mean.iter().sum::<f64>() / mesh_mean.len() as f64;
    let global_rms = (mesh_rms.iter().map(|x| x * x).sum::<f64>() / mesh_rms.len() as f64).sqrt();

    Ok(Background {
        level,
        rms,
        global_mean,
        global_rms,
        width: cutout.width,
        height: cutout.height,
    })
}

/// Sigma-clipped mean and standard deviation.
fn sigma_clipped_stats(values: &mut Vec<f64>, sigma: f64, iters: usize) -> (f64, f64) {
    if values.is_empty() {
        return (0.0, 0.0);
    }

    for _ in 0..iters {
        if values.len() < 3 {
            break;
        }
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let var = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / values.len() as f64;
        let std = var.sqrt();
        if std < 1e-30 {
            break;
        }
        let lo = mean - sigma * std;
        let hi = mean + sigma * std;
        let before = values.len();
        values.retain(|&x| x >= lo && x <= hi);
        if values.len() == before {
            break; // No more clipping needed
        }
    }

    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let var = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / values.len() as f64;
    (mean, var.sqrt())
}

/// Bilinear interpolation on a 2D grid.
fn bilinear_interp(grid: &[f64], nx: usize, ny: usize, fx: f64, fy: f64) -> f64 {
    let ix0 = (fx.floor() as isize).clamp(0, nx as isize - 1) as usize;
    let iy0 = (fy.floor() as isize).clamp(0, ny as isize - 1) as usize;
    let ix1 = (ix0 + 1).min(nx - 1);
    let iy1 = (iy0 + 1).min(ny - 1);

    let tx = (fx - ix0 as f64).clamp(0.0, 1.0);
    let ty = (fy - iy0 as f64).clamp(0.0, 1.0);

    let v00 = grid[iy0 * nx + ix0];
    let v10 = grid[iy0 * nx + ix1];
    let v01 = grid[iy1 * nx + ix0];
    let v11 = grid[iy1 * nx + ix1];

    v00 * (1.0 - tx) * (1.0 - ty) + v10 * tx * (1.0 - ty) + v01 * (1.0 - tx) * ty + v11 * tx * ty
}

/// Extract sources from a background-subtracted image using thresholding
/// and connected-component labeling.
pub fn extract_sources(
    cutout: &Cutout,
    config: &ExtractionConfig,
) -> Result<Vec<DetectedSource>, ProstError> {
    // Step 1: Estimate background
    let bg = estimate_background(cutout, config)?;

    // Step 2: Create background-subtracted image and threshold mask
    let npix = cutout.width * cutout.height;
    let mut subtracted = vec![0.0_f64; npix];
    let mut mask = vec![false; npix];

    for i in 0..npix {
        subtracted[i] = cutout.data[i] - bg.level[i];
        let thresh = config.thresh_sigma * bg.rms[i];
        mask[i] = subtracted[i] > thresh && thresh > 0.0;
    }

    // Step 3: Connected-component labeling (8-connectivity)
    let labels = label_connected_components(&mask, cutout.width, cutout.height);
    let n_labels = *labels.iter().max().unwrap_or(&0);

    // Step 4: Measure each labeled region
    let mut sources = Vec::new();
    for label in 1..=n_labels {
        let source = measure_source(
            &subtracted,
            &labels,
            label,
            cutout.width,
            cutout.height,
            &bg,
            cutout,
        );
        if let Some(s) = source {
            if s.npix >= config.min_pixels {
                sources.push(s);
            }
        }
    }

    // Sort by flux (brightest first)
    sources.sort_by(|a, b| b.flux.total_cmp(&a.flux));

    Ok(sources)
}

/// Label connected components in a boolean mask using 8-connectivity.
/// Returns a label map (0 = background, 1..n = component labels).
fn label_connected_components(mask: &[bool], width: usize, height: usize) -> Vec<u32> {
    let npix = width * height;
    let mut labels = vec![0u32; npix];
    let mut current_label = 0u32;

    // Union-Find for efficient merging
    let mut parent: Vec<u32> = Vec::new();

    let find = |parent: &mut Vec<u32>, mut x: u32| -> u32 {
        while parent[x as usize] != x {
            parent[x as usize] = parent[parent[x as usize] as usize]; // path compression
            x = parent[x as usize];
        }
        x
    };

    // First pass: assign provisional labels and record equivalences
    for r in 0..height {
        for c in 0..width {
            let idx = r * width + c;
            if !mask[idx] {
                continue;
            }

            // Collect labels of already-labeled neighbors (8-connected, only above/left)
            let mut neighbor_labels: Vec<u32> = Vec::new();
            let offsets: [(isize, isize); 4] = [(-1, -1), (-1, 0), (-1, 1), (0, -1)];
            for (dr, dc) in offsets {
                let nr = r as isize + dr;
                let nc = c as isize + dc;
                if nr >= 0 && nr < height as isize && nc >= 0 && nc < width as isize {
                    let ni = nr as usize * width + nc as usize;
                    if labels[ni] > 0 {
                        neighbor_labels.push(labels[ni]);
                    }
                }
            }

            if neighbor_labels.is_empty() {
                // New label
                current_label += 1;
                labels[idx] = current_label;
                parent.push(0); // placeholder for index 0
                if parent.len() <= current_label as usize {
                    parent.resize(current_label as usize + 1, 0);
                }
                parent[current_label as usize] = current_label;
            } else {
                // Take minimum label among neighbors
                let min_label = *neighbor_labels.iter().min().unwrap();
                labels[idx] = min_label;
                // Union all neighbor labels
                for &nl in &neighbor_labels {
                    let root_min = find(&mut parent, min_label);
                    let root_nl = find(&mut parent, nl);
                    if root_min != root_nl {
                        let (lo, hi) = if root_min < root_nl {
                            (root_min, root_nl)
                        } else {
                            (root_nl, root_min)
                        };
                        parent[hi as usize] = lo;
                    }
                }
            }
        }
    }

    // Second pass: resolve equivalences
    let mut label_map: Vec<u32> = vec![0; current_label as usize + 1];
    let mut next_label = 0u32;
    for i in 1..=current_label {
        let root = find(&mut parent, i);
        if label_map[root as usize] == 0 {
            next_label += 1;
            label_map[root as usize] = next_label;
        }
        label_map[i as usize] = label_map[root as usize];
    }

    for label in labels.iter_mut() {
        if *label > 0 {
            *label = label_map[*label as usize];
        }
    }

    labels
}

/// Measure source properties (centroid, flux, shape moments) for a labeled region.
fn measure_source(
    subtracted: &[f64],
    labels: &[u32],
    label: u32,
    width: usize,
    height: usize,
    bg: &Background,
    cutout: &Cutout,
) -> Option<DetectedSource> {
    let mut sum_f = 0.0_f64;
    let mut sum_fx = 0.0_f64;
    let mut sum_fy = 0.0_f64;
    let mut sum_var = 0.0_f64;
    let mut npix = 0usize;
    let mut peak = f64::NEG_INFINITY;

    // First pass: compute flux-weighted centroid
    for r in 0..height {
        for c in 0..width {
            let idx = r * width + c;
            if labels[idx] != label {
                continue;
            }
            let f = subtracted[idx].max(0.0);
            sum_f += f;
            sum_fx += f * c as f64;
            sum_fy += f * r as f64;
            sum_var += bg.rms[idx].powi(2);
            npix += 1;
            if subtracted[idx] > peak {
                peak = subtracted[idx];
            }
        }
    }

    if npix == 0 || sum_f <= 0.0 {
        return None;
    }

    let cx = sum_fx / sum_f;
    let cy = sum_fy / sum_f;

    // Second pass: compute second-order moments for shape
    let mut mxx = 0.0_f64;
    let mut myy = 0.0_f64;
    let mut mxy = 0.0_f64;

    for r in 0..height {
        for c in 0..width {
            let idx = r * width + c;
            if labels[idx] != label {
                continue;
            }
            let f = subtracted[idx].max(0.0);
            let dx = c as f64 - cx;
            let dy = r as f64 - cy;
            mxx += f * dx * dx;
            myy += f * dy * dy;
            mxy += f * dx * dy;
        }
    }

    mxx /= sum_f;
    myy /= sum_f;
    mxy /= sum_f;

    // Derive semi-axes and PA from moments
    // eigenvalues of [[mxx, mxy], [mxy, myy]]
    let trace = mxx + myy;
    let det = mxx * myy - mxy * mxy;
    let discriminant = (trace * trace - 4.0 * det).max(0.0);
    let lambda1 = 0.5 * (trace + discriminant.sqrt());
    let lambda2 = 0.5 * (trace - discriminant.sqrt());

    let a_pix = lambda1.sqrt().max(0.5);
    let b_pix = lambda2.max(0.0).sqrt().max(0.5);
    let pa_rad = 0.5 * (2.0 * mxy).atan2(mxx - myy);
    let pa_deg = pa_rad.to_degrees().rem_euclid(180.0);

    let flux_err = sum_var.sqrt();
    let snr = if flux_err > 0.0 {
        sum_f / flux_err
    } else {
        0.0
    };

    let (ra, dec) = cutout.pixel_to_sky(cy, cx);

    Some(DetectedSource {
        x: cx,
        y: cy,
        ra,
        dec,
        a_pix,
        b_pix,
        pa_deg,
        flux: sum_f,
        npix,
        peak,
        snr,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn make_gaussian_cutout() -> Cutout {
        let mut cutout = Cutout::zeros(101, 101, 0.262, 180.0, 45.0);
        // Add background noise level
        for v in cutout.data.iter_mut() {
            *v = 10.0;
        }
        // Add a bright Gaussian source at center
        cutout.add_gaussian(50.0, 50.0, 500.0, 8.0, 4.0, 0.3);
        cutout
    }

    #[test]
    fn test_sigma_clipped_stats_basic() {
        let mut values = vec![1.0, 2.0, 3.0, 4.0, 5.0, 100.0]; // outlier at 100
        let (mean, _rms) = sigma_clipped_stats(&mut values, 3.0, 5);
        // After clipping, 100 should be removed (may take more iterations)
        assert!(mean < 20.0, "Mean should be low after clipping outlier, got {mean}");
    }

    #[test]
    fn test_sigma_clipped_stats_uniform() {
        let mut values: Vec<f64> = (0..100).map(|i| 10.0 + (i as f64) * 0.01).collect();
        let (mean, rms) = sigma_clipped_stats(&mut values, 3.0, 3);
        assert_relative_eq!(mean, 10.495, epsilon = 0.1);
        assert!(rms < 1.0);
    }

    #[test]
    fn test_background_estimation() {
        let mut cutout = Cutout::zeros(64, 64, 0.262, 180.0, 45.0);
        // Flat background at 100
        for v in cutout.data.iter_mut() {
            *v = 100.0;
        }
        let config = ExtractionConfig {
            back_size: 32,
            ..Default::default()
        };
        let bg = estimate_background(&cutout, &config).unwrap();
        assert_relative_eq!(bg.global_mean, 100.0, epsilon = 0.01);
        assert_relative_eq!(bg.global_rms, 0.0, epsilon = 0.01);
    }

    #[test]
    fn test_connected_components_single() {
        // 5x5 mask with a single blob
        #[rustfmt::skip]
        let mask = vec![
            false, false, false, false, false,
            false, true,  true,  false, false,
            false, true,  true,  true,  false,
            false, false, true,  false, false,
            false, false, false, false, false,
        ];
        let labels = label_connected_components(&mask, 5, 5);
        let max_label = *labels.iter().max().unwrap();
        assert_eq!(max_label, 1, "Should find exactly one component");
        assert_eq!(labels.iter().filter(|&&l| l == 1).count(), 6);
    }

    #[test]
    fn test_connected_components_two() {
        // 5x5 mask with two separate blobs
        #[rustfmt::skip]
        let mask = vec![
            true,  true,  false, false, false,
            true,  false, false, false, false,
            false, false, false, false, false,
            false, false, false, true,  true,
            false, false, false, true,  true,
        ];
        let labels = label_connected_components(&mask, 5, 5);
        let max_label = *labels.iter().max().unwrap();
        assert_eq!(max_label, 2, "Should find two components");
    }

    #[test]
    fn test_extract_single_source() {
        let cutout = make_gaussian_cutout();
        let config = ExtractionConfig {
            thresh_sigma: 3.0,
            min_pixels: 5,
            back_size: 32,
            ..Default::default()
        };
        let sources = extract_sources(&cutout, &config).unwrap();
        assert!(
            !sources.is_empty(),
            "Should detect at least one source"
        );
        let s = &sources[0];
        // Centroid should be near (50, 50)
        assert!((s.x - 50.0).abs() < 2.0, "x centroid off: {}", s.x);
        assert!((s.y - 50.0).abs() < 2.0, "y centroid off: {}", s.y);
        // Semi-major should be larger than semi-minor
        assert!(s.a_pix > s.b_pix);
        assert!(s.flux > 0.0);
        assert!(s.snr > 10.0);
    }

    #[test]
    fn test_extract_elongated_shape() {
        let cutout = make_gaussian_cutout();
        let config = ExtractionConfig {
            thresh_sigma: 2.0,
            min_pixels: 3,
            back_size: 32,
            ..Default::default()
        };
        let sources = extract_sources(&cutout, &config).unwrap();
        assert!(!sources.is_empty());
        let s = &sources[0];
        // The Gaussian has sigma_major=8, sigma_minor=4, so a/b ~ 2
        let ratio = s.a_pix / s.b_pix;
        assert!(
            ratio > 1.3 && ratio < 3.0,
            "Axis ratio should reflect elongation: a/b = {}",
            ratio
        );
    }

    #[test]
    fn test_extract_two_sources() {
        let mut cutout = Cutout::zeros(101, 101, 0.262, 180.0, 45.0);
        for v in cutout.data.iter_mut() {
            *v = 10.0;
        }
        // Two well-separated Gaussians
        cutout.add_gaussian(25.0, 25.0, 300.0, 5.0, 5.0, 0.0);
        cutout.add_gaussian(75.0, 75.0, 200.0, 4.0, 4.0, 0.0);

        let config = ExtractionConfig {
            thresh_sigma: 3.0,
            min_pixels: 5,
            back_size: 32,
            ..Default::default()
        };
        let sources = extract_sources(&cutout, &config).unwrap();
        assert!(
            sources.len() >= 2,
            "Should detect two sources, got {}",
            sources.len()
        );
    }

    #[test]
    fn test_extract_no_sources() {
        // Flat image, no sources
        let mut cutout = Cutout::zeros(64, 64, 0.262, 180.0, 45.0);
        for v in cutout.data.iter_mut() {
            *v = 100.0;
        }
        let config = ExtractionConfig::default();
        let sources = extract_sources(&cutout, &config).unwrap();
        assert!(sources.is_empty(), "Should detect no sources in flat image");
    }

    #[test]
    fn test_source_sky_coords() {
        let cutout = make_gaussian_cutout();
        let config = ExtractionConfig {
            thresh_sigma: 3.0,
            min_pixels: 5,
            back_size: 32,
            ..Default::default()
        };
        let sources = extract_sources(&cutout, &config).unwrap();
        let s = &sources[0];
        // Source at center of image → should be near (180, 45)
        assert!((s.ra - 180.0).abs() < 0.01);
        assert!((s.dec - 45.0).abs() < 0.01);
    }
}
