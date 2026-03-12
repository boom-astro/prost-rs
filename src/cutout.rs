//! Image cutout representation and WCS transformations.
//!
//! Provides the `Cutout` struct for working with image data: pixel access,
//! coordinate transformations, and sub-image extraction. Cutout retrieval
//! from survey servers is handled by the integration tests and Python layer.

use crate::errors::ProstError;

/// Supported imaging surveys for cutout retrieval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImagingSurvey {
    /// DESI Legacy Imaging Surveys (DECaLS / MzLS / BASS)
    LegacySurvey,
    /// Pan-STARRS1 3pi survey
    PanSTARRS,
    /// SkyMapper Southern Survey
    SkyMapper,
}

/// A 2D image cutout with WCS metadata.
#[derive(Debug, Clone)]
pub struct Cutout {
    /// Pixel data (row-major, flattened)
    pub data: Vec<f64>,
    /// Image width in pixels (number of columns)
    pub width: usize,
    /// Image height in pixels (number of rows)
    pub height: usize,
    /// Pixel scale in arcsec/pixel
    pub pixel_scale: f64,
    /// RA of image center (degrees)
    pub center_ra: f64,
    /// Dec of image center (degrees)
    pub center_dec: f64,
    /// Source survey
    pub survey: ImagingSurvey,
}

impl Cutout {
    /// Create a new cutout from raw pixel data.
    pub fn new(
        data: Vec<f64>,
        width: usize,
        height: usize,
        pixel_scale: f64,
        center_ra: f64,
        center_dec: f64,
        survey: ImagingSurvey,
    ) -> Result<Self, ProstError> {
        if data.len() != width * height {
            return Err(ProstError::CutoutError(format!(
                "data length {} does not match {}x{}",
                data.len(),
                width,
                height
            )));
        }
        Ok(Self {
            data,
            width,
            height,
            pixel_scale,
            center_ra,
            center_dec,
            survey,
        })
    }

    /// Create a cutout filled with zeros (for testing).
    pub fn zeros(
        width: usize,
        height: usize,
        pixel_scale: f64,
        center_ra: f64,
        center_dec: f64,
    ) -> Self {
        Self {
            data: vec![0.0; width * height],
            width,
            height,
            pixel_scale,
            center_ra,
            center_dec,
            survey: ImagingSurvey::LegacySurvey,
        }
    }

    /// Get pixel value at (row, col). Returns None if out of bounds.
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> Option<f64> {
        if row < self.height && col < self.width {
            Some(self.data[row * self.width + col])
        } else {
            None
        }
    }

    /// Set pixel value at (row, col). Returns false if out of bounds.
    #[inline]
    pub fn set(&mut self, row: usize, col: usize, value: f64) -> bool {
        if row < self.height && col < self.width {
            self.data[row * self.width + col] = value;
            true
        } else {
            false
        }
    }

    /// Convert pixel coordinates to sky coordinates (TAN projection).
    pub fn pixel_to_sky(&self, row: f64, col: f64) -> (f64, f64) {
        let dx = (col - self.width as f64 / 2.0) * self.pixel_scale / 3600.0;
        let dy = (row - self.height as f64 / 2.0) * self.pixel_scale / 3600.0;
        let dec = self.center_dec + dy;
        let ra = self.center_ra + dx / dec.to_radians().cos();
        (ra, dec)
    }

    /// Convert sky coordinates to pixel coordinates (TAN projection).
    pub fn sky_to_pixel(&self, ra: f64, dec: f64) -> (f64, f64) {
        let dx = (ra - self.center_ra) * dec.to_radians().cos() * 3600.0 / self.pixel_scale;
        let dy = (dec - self.center_dec) * 3600.0 / self.pixel_scale;
        let col = dx + self.width as f64 / 2.0;
        let row = dy + self.height as f64 / 2.0;
        (row, col)
    }

    /// Size of the cutout in arcsec (width, height).
    pub fn size_arcsec(&self) -> (f64, f64) {
        (
            self.width as f64 * self.pixel_scale,
            self.height as f64 * self.pixel_scale,
        )
    }

    /// Extract a sub-image centered on a pixel position.
    /// Returns None if the sub-image would extend beyond the cutout.
    pub fn sub_image(&self, center_row: usize, center_col: usize, radius: usize) -> Option<Cutout> {
        let r0 = center_row.checked_sub(radius)?;
        let c0 = center_col.checked_sub(radius)?;
        let r1 = center_row + radius + 1;
        let c1 = center_col + radius + 1;
        if r1 > self.height || c1 > self.width {
            return None;
        }
        let sub_w = c1 - c0;
        let sub_h = r1 - r0;
        let mut sub_data = Vec::with_capacity(sub_w * sub_h);
        for r in r0..r1 {
            for c in c0..c1 {
                sub_data.push(self.data[r * self.width + c]);
            }
        }
        let (ra, dec) = self.pixel_to_sky(center_row as f64, center_col as f64);
        Some(Cutout {
            data: sub_data,
            width: sub_w,
            height: sub_h,
            pixel_scale: self.pixel_scale,
            center_ra: ra,
            center_dec: dec,
            survey: self.survey,
        })
    }

    /// Add a 2D Gaussian to the cutout (for generating test images).
    pub fn add_gaussian(
        &mut self,
        center_row: f64,
        center_col: f64,
        amplitude: f64,
        sigma_major: f64,
        sigma_minor: f64,
        pa_rad: f64,
    ) {
        let (sin_pa, cos_pa) = pa_rad.sin_cos();
        for r in 0..self.height {
            for c in 0..self.width {
                let dy = r as f64 - center_row;
                let dx = c as f64 - center_col;
                let x_rot = dx * cos_pa + dy * sin_pa;
                let y_rot = -dx * sin_pa + dy * cos_pa;
                let exponent =
                    -0.5 * ((x_rot / sigma_major).powi(2) + (y_rot / sigma_minor).powi(2));
                self.data[r * self.width + c] += amplitude * exponent.exp();
            }
        }
    }

    /// Add a 2D Sérsic profile to the cutout (for generating test images).
    #[allow(clippy::too_many_arguments)]
    pub fn add_sersic(
        &mut self,
        center_row: f64,
        center_col: f64,
        i_eff: f64,
        r_eff_pix: f64,
        n: f64,
        axis_ratio: f64,
        pa_rad: f64,
    ) {
        let b_n = 1.9992 * n - 0.3271;
        let (sin_pa, cos_pa) = pa_rad.sin_cos();
        for r in 0..self.height {
            for c in 0..self.width {
                let dy = r as f64 - center_row;
                let dx = c as f64 - center_col;
                let x_rot = dx * cos_pa + dy * sin_pa;
                let y_rot = -dx * sin_pa + dy * cos_pa;
                let r_ell =
                    (x_rot * x_rot + (y_rot / axis_ratio).powi(2)).sqrt();
                if r_ell < 1e-10 {
                    self.data[r * self.width + c] += i_eff * b_n.exp();
                } else {
                    let ratio = r_ell / r_eff_pix;
                    self.data[r * self.width + c] +=
                        i_eff * (-b_n * (ratio.powf(1.0 / n) - 1.0)).exp();
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_pixel_to_sky_center() {
        let cutout = Cutout::zeros(10, 10, 0.262, 180.0, 45.0);
        let (ra, dec) = cutout.pixel_to_sky(5.0, 5.0);
        assert!((ra - 180.0).abs() < 1e-6);
        assert!((dec - 45.0).abs() < 1e-6);
    }

    #[test]
    fn test_sky_to_pixel_roundtrip() {
        let cutout = Cutout::zeros(10, 10, 0.262, 180.0, 45.0);
        let (ra, dec) = cutout.pixel_to_sky(3.0, 7.0);
        let (row, col) = cutout.sky_to_pixel(ra, dec);
        assert!((row - 3.0).abs() < 1e-6);
        assert!((col - 7.0).abs() < 1e-6);
    }

    #[test]
    fn test_sub_image() {
        let mut cutout = Cutout::zeros(20, 20, 0.262, 180.0, 45.0);
        cutout.set(10, 10, 99.0);
        let sub = cutout.sub_image(10, 10, 2).unwrap();
        assert_eq!(sub.width, 5);
        assert_eq!(sub.height, 5);
        // Center of sub-image should have the value we set
        assert_relative_eq!(sub.get(2, 2).unwrap(), 99.0);
    }

    #[test]
    fn test_add_gaussian() {
        let mut cutout = Cutout::zeros(51, 51, 0.262, 180.0, 45.0);
        cutout.add_gaussian(25.0, 25.0, 100.0, 5.0, 3.0, 0.0);
        // Center should be brightest
        let center = cutout.get(25, 25).unwrap();
        let edge = cutout.get(25, 0).unwrap();
        assert!(center > 99.0);
        assert!(edge < 1.0);
    }

    #[test]
    fn test_add_sersic() {
        let mut cutout = Cutout::zeros(51, 51, 0.262, 180.0, 45.0);
        cutout.add_sersic(25.0, 25.0, 100.0, 10.0, 1.0, 1.0, 0.0);
        let center = cutout.get(25, 25).unwrap();
        let r_eff = cutout.get(25, 35).unwrap(); // 10 pixels from center
        // At r_eff, I should be I_e = 100
        assert_relative_eq!(r_eff, 100.0, epsilon = 1.0);
        assert!(center > r_eff);
    }

    #[test]
    fn test_size_arcsec() {
        let cutout = Cutout::zeros(100, 80, 0.262, 0.0, 0.0);
        let (w, h) = cutout.size_arcsec();
        assert_relative_eq!(w, 26.2, epsilon = 0.01);
        assert_relative_eq!(h, 20.96, epsilon = 0.01);
    }
}
