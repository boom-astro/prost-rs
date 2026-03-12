"""Tests for the prost Python bindings."""

import numpy as np
import pytest

from prost import (
    AssociationConfig,
    AssociationResult,
    Cutout,
    DetectedSource,
    ExtractionConfig,
    FitConfig,
    GalaxyCandidate,
    HostCandidate,
    SersicFit,
    Transient,
    associate_host,
    compute_dlr,
    elliptical_radius,
    estimate_background,
    extract_sources,
    fit_sersic,
    offset_likelihood,
    redshift_likelihood,
    sersic_profile,
)


# ---------------------------------------------------------------------------
# Transient
# ---------------------------------------------------------------------------

class TestTransient:
    def test_basic(self):
        t = Transient(ra=180.0, dec=45.0)
        assert t.ra == 180.0
        assert t.dec == 45.0
        assert t.redshift is None

    def test_with_redshift(self):
        t = Transient(ra=0.0, dec=0.0, redshift=0.05, redshift_err=0.001)
        assert t.redshift == pytest.approx(0.05)
        assert t.redshift_err == pytest.approx(0.001)

    def test_repr(self):
        t = Transient(ra=180.0, dec=45.0, name="SN2020test")
        assert "180" in repr(t)


# ---------------------------------------------------------------------------
# GalaxyCandidate
# ---------------------------------------------------------------------------

class TestGalaxyCandidate:
    def test_basic(self):
        g = GalaxyCandidate(
            ra=180.0, dec=45.0,
            a_arcsec=5.0, b_arcsec=3.0, pa_deg=30.0,
        )
        assert g.a_arcsec == pytest.approx(5.0)
        assert g.b_arcsec == pytest.approx(3.0)
        assert g.pa_deg == pytest.approx(30.0)

    def test_from_tractor(self):
        g = GalaxyCandidate.from_tractor(
            ra=180.0, dec=45.0,
            shape_r=2.5, shape_e1=0.3, shape_e2=0.1,
        )
        assert g.a_arcsec > 0
        assert g.b_arcsec > 0
        assert g.a_arcsec >= g.b_arcsec

    def test_from_tractor_invalid(self):
        with pytest.raises(ValueError):
            GalaxyCandidate.from_tractor(
                ra=0.0, dec=0.0,
                shape_r=-1.0, shape_e1=0.0, shape_e2=0.0,
            )


# ---------------------------------------------------------------------------
# Association
# ---------------------------------------------------------------------------

class TestAssociation:
    def test_single_galaxy(self):
        t = Transient(ra=180.0, dec=45.0)
        g = GalaxyCandidate(
            ra=180.0 + 2.0 / 3600, dec=45.0,
            a_arcsec=5.0, b_arcsec=3.0, pa_deg=0.0,
        )
        result = associate_host(t, [g])
        assert len(result) == 1
        assert result.best_host is not None
        assert result.p_none >= 0
        total = sum(c.posterior for c in result.candidates) + result.p_none
        assert total == pytest.approx(1.0, abs=1e-6)

    def test_posteriors_sum_to_one(self):
        t = Transient(ra=180.0, dec=45.0)
        galaxies = [
            GalaxyCandidate(
                ra=180.0 + i / 3600, dec=45.0,
                a_arcsec=3.0, b_arcsec=2.0, pa_deg=0.0,
            )
            for i in range(1, 4)
        ]
        result = associate_host(t, galaxies)
        total = sum(c.posterior for c in result.candidates) + result.p_none
        assert total == pytest.approx(1.0, abs=1e-6)

    def test_redshift_breaks_degeneracy(self):
        t = Transient(ra=180.0, dec=45.0, redshift=0.05, redshift_err=0.001)
        g_match = GalaxyCandidate(
            ra=180.0 + 3.0 / 3600, dec=45.0,
            a_arcsec=5.0, b_arcsec=3.0, pa_deg=0.0,
            redshift=0.05, redshift_err=0.002,
            objname="match",
        )
        g_wrong = GalaxyCandidate(
            ra=180.0 + 2.5 / 3600, dec=45.0,
            a_arcsec=5.0, b_arcsec=3.0, pa_deg=0.0,
            redshift=0.5, redshift_err=0.002,
            objname="wrong",
        )
        result = associate_host(t, [g_match, g_wrong])
        assert result.best_host.galaxy.objname == "match"

    def test_empty_candidates(self):
        t = Transient(ra=180.0, dec=45.0)
        with pytest.raises(ValueError):
            associate_host(t, [])

    def test_config(self):
        cfg = AssociationConfig(max_fractional_offset=5.0, max_candidates=3)
        assert cfg.max_fractional_offset == pytest.approx(5.0)
        assert cfg.max_candidates == 3


# ---------------------------------------------------------------------------
# DLR / Likelihood
# ---------------------------------------------------------------------------

class TestDLR:
    def test_compute_dlr(self):
        d = compute_dlr(
            transient_ra=180.0, transient_dec=45.0,
            galaxy_ra=180.001, galaxy_dec=45.0,
            a_arcsec=5.0, b_arcsec=3.0, pa_deg=0.0,
        )
        assert "separation_arcsec" in d
        assert "directional_radius" in d
        assert "fractional_offset" in d
        assert d["separation_arcsec"] > 0
        assert d["fractional_offset"] > 0

    def test_offset_likelihood_scalar(self):
        l = offset_likelihood(1.0)
        assert l > 0

    def test_offset_likelihood_array(self):
        fo = np.array([0.1, 0.5, 1.0, 2.0, 5.0])
        result = offset_likelihood(fo)
        assert isinstance(result, np.ndarray)
        assert len(result) == 5
        assert all(r > 0 for r in result)

    def test_redshift_likelihood_matching(self):
        l = redshift_likelihood(0.05, 0.001, 0.05, 0.001)
        assert l == pytest.approx(1.0, abs=0.01)

    def test_redshift_likelihood_none(self):
        l = redshift_likelihood(None, None, 0.05, 0.001)
        assert l == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Morphology utilities
# ---------------------------------------------------------------------------

class TestMorphology:
    def test_sersic_at_r_eff(self):
        val = sersic_profile(r=1.0, r_eff=1.0, n=1.0, i_eff=100.0)
        assert val == pytest.approx(100.0, rel=1e-4)

    def test_sersic_decreasing(self):
        i1 = sersic_profile(0.5, 2.0, 2.0, 100.0)
        i2 = sersic_profile(2.0, 2.0, 2.0, 100.0)
        i3 = sersic_profile(5.0, 2.0, 2.0, 100.0)
        assert i1 > i2 > i3

    def test_elliptical_radius_circular(self):
        r = elliptical_radius(3.0, 4.0, 1.0, 0.0)
        assert r == pytest.approx(5.0, rel=1e-6)


# ---------------------------------------------------------------------------
# Cutout
# ---------------------------------------------------------------------------

class TestCutout:
    def test_zeros(self):
        c = Cutout.zeros(64, 64, 0.262, 180.0, 45.0)
        assert c.width == 64
        assert c.height == 64
        assert c.pixel_scale == pytest.approx(0.262)

    def test_from_array(self):
        data = np.ones((32, 32), dtype=np.float64) * 100.0
        c = Cutout(data=data, pixel_scale=0.262,
                   center_ra=180.0, center_dec=45.0)
        assert c.width == 32
        assert c.height == 32
        assert c.get(0, 0) == pytest.approx(100.0)

    def test_get_set(self):
        c = Cutout.zeros(10, 10, 0.262, 0.0, 0.0)
        assert c.get(5, 5) == pytest.approx(0.0)
        c.set(5, 5, 42.0)
        assert c.get(5, 5) == pytest.approx(42.0)

    def test_pixel_sky_roundtrip(self):
        c = Cutout.zeros(64, 64, 0.262, 180.0, 45.0)
        ra, dec = c.pixel_to_sky(32.0, 32.0)
        row, col = c.sky_to_pixel(ra, dec)
        assert row == pytest.approx(32.0, abs=0.01)
        assert col == pytest.approx(32.0, abs=0.01)

    def test_size_arcsec(self):
        c = Cutout.zeros(100, 80, 0.262, 0.0, 0.0)
        w, h = c.size_arcsec()
        assert w == pytest.approx(26.2, abs=0.01)
        assert h == pytest.approx(20.96, abs=0.01)

    def test_data_property(self):
        data = np.arange(16, dtype=np.float64).reshape(4, 4)
        c = Cutout(data=data, pixel_scale=1.0, center_ra=0.0, center_dec=0.0)
        out = c.data
        assert out.shape == (4, 4)
        assert out[0, 0] == pytest.approx(0.0)
        assert out[3, 3] == pytest.approx(15.0)

    def test_sub_image(self):
        c = Cutout.zeros(20, 20, 0.262, 0.0, 0.0)
        c.set(10, 10, 99.0)
        sub = c.sub_image(10, 10, 2)
        assert sub is not None
        assert sub.width == 5
        assert sub.height == 5

    def test_repr(self):
        c = Cutout.zeros(64, 64, 0.262, 180.0, 45.0)
        assert "64x64" in repr(c)


# ---------------------------------------------------------------------------
# Source extraction
# ---------------------------------------------------------------------------

class TestExtraction:
    @staticmethod
    def _make_cutout_with_source():
        c = Cutout.zeros(101, 101, 0.262, 180.0, 45.0)
        # Flat background
        data = np.ones((101, 101), dtype=np.float64) * 10.0
        c = Cutout(data=data, pixel_scale=0.262,
                   center_ra=180.0, center_dec=45.0)
        c.add_gaussian(50.0, 50.0, 500.0, 8.0, 4.0, 0.3)
        return c

    def test_estimate_background(self):
        c = self._make_cutout_with_source()
        bg = estimate_background(c)
        assert bg.global_mean > 0
        assert bg.global_rms >= 0
        assert bg.width == 101
        assert bg.height == 101

    def test_background_arrays(self):
        c = self._make_cutout_with_source()
        bg = estimate_background(c)
        assert bg.level.shape == (101, 101)
        assert bg.rms.shape == (101, 101)

    def test_extract_sources(self):
        c = self._make_cutout_with_source()
        cfg = ExtractionConfig(thresh_sigma=3.0, min_pixels=5)
        sources = extract_sources(c, cfg)
        assert len(sources) >= 1
        s = sources[0]
        assert abs(s.x - 50.0) < 3.0
        assert abs(s.y - 50.0) < 3.0
        assert s.flux > 0
        assert s.snr > 0

    def test_extract_no_sources(self):
        data = np.ones((64, 64), dtype=np.float64) * 100.0
        c = Cutout(data=data, pixel_scale=0.262,
                   center_ra=0.0, center_dec=0.0)
        sources = extract_sources(c)
        assert len(sources) == 0

    def test_detected_source_methods(self):
        c = self._make_cutout_with_source()
        sources = extract_sources(c, ExtractionConfig(thresh_sigma=3.0))
        s = sources[0]
        assert s.a_arcsec(0.262) == pytest.approx(s.a_pix * 0.262)
        assert s.b_arcsec(0.262) == pytest.approx(s.b_pix * 0.262)
        assert 0.0 <= s.axis_ratio() <= 1.0
        assert 0.0 <= s.ellipticity() <= 1.0

    def test_extraction_config_defaults(self):
        cfg = ExtractionConfig()
        assert cfg.thresh_sigma == pytest.approx(1.5)
        assert cfg.min_pixels == 5


# ---------------------------------------------------------------------------
# Sersic fitting
# ---------------------------------------------------------------------------

class TestSersicFit:
    def test_fit_exponential(self):
        data = np.ones((101, 101), dtype=np.float64) * 50.0
        c = Cutout(data=data, pixel_scale=0.262,
                   center_ra=180.0, center_dec=45.0)
        c.add_sersic(50.0, 50.0, 200.0, 8.0, 1.0, 0.6, 0.5)

        sources = extract_sources(c, ExtractionConfig(thresh_sigma=2.0))
        assert len(sources) >= 1

        fit = fit_sersic(c, sources[0])
        assert abs(fit.x - 50.0) < 3.0
        assert abs(fit.y - 50.0) < 3.0
        assert fit.r_eff_pix > 0
        assert fit.r_eff > 0
        assert 0.3 <= fit.n <= 8.0
        assert 0.0 < fit.axis_ratio <= 1.0
        assert fit.ndata > 0
        assert fit.chi2_reduced >= 0

    def test_fit_config(self):
        cfg = FitConfig(max_sersic_n=6.0, min_sersic_n=0.5)
        assert cfg.max_sersic_n == pytest.approx(6.0)
        assert cfg.min_sersic_n == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# End-to-end: image -> association
# ---------------------------------------------------------------------------

class TestEndToEnd:
    def test_image_to_association(self):
        """Full pipeline: synthetic image -> extract -> fit -> associate."""
        data = np.ones((101, 101), dtype=np.float64) * 20.0
        c = Cutout(data=data, pixel_scale=0.262,
                   center_ra=180.0, center_dec=45.0)
        # Galaxy at center
        c.add_sersic(50.0, 50.0, 300.0, 10.0, 1.0, 0.7, 0.0)

        sources = extract_sources(c, ExtractionConfig(thresh_sigma=2.0))
        assert len(sources) >= 1

        fit = fit_sersic(c, sources[0])

        a = fit.r_eff / fit.axis_ratio
        b = fit.r_eff
        if a < b:
            a, b = b, a

        galaxy = GalaxyCandidate(
            ra=fit.ra, dec=fit.dec,
            a_arcsec=a, b_arcsec=b,
            pa_deg=fit.pa_deg,
            shape_from_image=True,
        )

        # Transient offset from galaxy center
        transient = Transient(ra=180.0 + 1.0 / 3600, dec=45.0)
        result = associate_host(transient, [galaxy])

        assert len(result) == 1
        assert result.best_host is not None
        assert result.best_host.posterior > 0
        total = sum(c.posterior for c in result.candidates) + result.p_none
        assert total == pytest.approx(1.0, abs=1e-6)
