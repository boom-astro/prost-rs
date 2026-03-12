# Changelog

## 0.1.0 (2026-03-13)

Initial release.

### Catalog-based association
- DLR-based geometric ranking with ellipse-aware separation
- Bayesian posterior scoring with Gamma(0.75) offset likelihood
- Gaussian redshift likelihood for breaking positional degeneracies
- Null hypothesis handling (hostless, unobserved, outside-radius)
- Legacy Survey Tractor shape parameter conversion (`from_tractor`)

### Image-based analysis
- `Cutout` struct for 2D image data with TAN-projection WCS
- Background estimation via sigma-clipped statistics on a mesh grid
- Source extraction with connected-component labeling (8-connectivity)
- Flux-weighted second-moment shape measurements
- 2D Sersic profile fitting via Levenberg-Marquardt optimizer
- Synthetic image generation (`add_gaussian`, `add_sersic`)

### Python bindings
- Full PyO3 bindings for all classes and functions
- NumPy integration for image arrays and vectorized likelihood
- `prost` Python package with `pip install` support via maturin

### Integration tests
- Known host validation (AT2017gfo, SN 2014J, SN 2011fe)
- BOOM API integration tests (ZTF transients + LS_DR10)
- Legacy Survey / Pan-STARRS cutout download tests
