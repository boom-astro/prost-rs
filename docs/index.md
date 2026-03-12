# prost

**Probabilistic host galaxy association for astronomical transients.**

prost identifies the most likely host galaxy for a transient event
(supernova, kilonova, GRB, etc.) using directional light radius (DLR)
ranking and Bayesian posterior scoring. It combines the fast geometric
ranking of GHOST3 with the probabilistic framework of Prost.

## Features

- **DLR-based ranking** — Ellipse-aware separation metric that accounts for galaxy morphology
- **Bayesian posteriors** — Gamma(0.75) offset likelihood with redshift and absolute magnitude components
- **Tractor shape support** — Direct ingestion of Legacy Survey shape parameters
- **Null hypothesis handling** — Proper p(none) for hostless, unobserved, and outside-radius scenarios
- **Image-based source extraction** — Detect galaxies directly from FITS cutouts when catalog coverage is missing
- **Sersic profile fitting** — Levenberg-Marquardt fitting of 2D Sersic models to derive morphology from images
- **Rust performance** — Core algorithms in Rust with Python bindings via PyO3

## Two Workflows

### Catalog-based association

When you have galaxy catalogs (e.g., Legacy Survey via BOOM), provide
positions and shape parameters directly:

```python
from prost import Transient, GalaxyCandidate, associate_host

transient = Transient(ra=197.4504, dec=-23.3815, redshift=0.00978)

galaxy = GalaxyCandidate(
    ra=197.4487, dec=-23.3839,
    a_arcsec=22.0, b_arcsec=16.0, pa_deg=75.0,
    redshift=0.00978, objname="NGC 4993",
)

result = associate_host(transient, [galaxy])
print(f"P(host) = {result.best_host.posterior:.3f}")
```

### Image-based association

When no catalog is available, extract galaxies from a FITS cutout and
fit their morphology:

```python
from prost import Cutout, extract_sources, fit_sersic, associate_host, Transient, GalaxyCandidate

cutout = Cutout(data=image_array, pixel_scale=0.262,
                center_ra=185.48, center_dec=4.47)

sources = extract_sources(cutout)
fit = fit_sersic(cutout, sources[0])

galaxy = GalaxyCandidate(
    ra=fit.ra, dec=fit.dec,
    a_arcsec=fit.r_eff / fit.axis_ratio,
    b_arcsec=fit.r_eff,
    pa_deg=fit.pa_deg, shape_from_image=True,
)

result = associate_host(Transient(ra=185.48, dec=4.47), [galaxy])
```

See the [Image-Based Association](examples/image-based-association.md)
example for the full pipeline with a real FITS cutout.

## Architecture

```
prost-rs/
├── src/                   # Rust core library
│   ├── associate.rs       # Bayesian host association
│   ├── cutout.rs          # Image cutout + WCS transforms
│   ├── dlr.rs             # Directional light radius
│   ├── ellipse.rs         # Ellipse geometry + Tractor conversion
│   ├── morphology.rs      # Sérsic fitting (Levenberg-Marquardt)
│   ├── source.rs          # Source extraction + background estimation
│   ├── likelihood.rs      # Offset + redshift likelihoods
│   ├── prior.rs           # Prior distributions
│   └── types.rs           # Transient, GalaxyCandidate, HostCandidate
├── python/                # PyO3 bindings (maturin)
│   └── src/lib.rs         # prost_extension module
├── prost/                 # Python package
│   ├── __init__.py        # Public API
│   └── _prost.py          # Extension wrapper
└── tests/                 # Rust integration tests
```

## References

- Gagliano et al. (2021) — [PROST](https://ui.adsabs.harvard.edu/abs/2021ApJ...908..170G): Probabilistic association framework
- Gagliano et al. (2021) — [GHOST](https://ui.adsabs.harvard.edu/abs/2021ApJ...911L..14G): Galaxy-transient association catalog
- Sullivan et al. (2006) — [DLR host matching](https://doi.org/10.1086/506137): Directional light radius methodology
