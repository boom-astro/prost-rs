# prost

**Probabilistic host galaxy association for astronomical transients.**

prost identifies the most likely host galaxy for a transient event
(supernova, kilonova, GRB, etc.) using directional light radius (DLR)
ranking and Bayesian posterior scoring.

## Features

- **DLR-based ranking** — Ellipse-aware separation accounting for galaxy morphology
- **Bayesian posteriors** — Gamma(0.75) offset likelihood with redshift scoring
- **Tractor shape support** — Direct ingestion of Legacy Survey parameters
- **Image-based extraction** — Detect galaxies from FITS cutouts when catalogs are missing
- **Sersic profile fitting** — LM fitting of 2D Sersic models for morphology
- **Rust performance** — Core algorithms in Rust with Python bindings via PyO3

## Quick Start

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

## Installation

```bash
# Requires Rust toolchain: https://rustup.rs/
git clone https://github.com/boom-astro/prost-rs.git
cd prost-rs
pip install ./python
pip install -e .
```

## Documentation

Full documentation: `mkdocs serve` (local) or see `docs/`.

- [Quick Start](docs/getting-started/quickstart.md)
- [Physics](docs/physics/association.md)
- [API Reference](docs/api/index.md)
- [Image-Based Association](docs/examples/image-based-association.md)
- [BOOM Integration](docs/examples/boom-integration.md)

## References

- Gagliano et al. (2021) — [PROST](https://doi.org/10.3847/1538-4357/abd02b)
- Gagliano et al. (2021) — [GHOST](https://doi.org/10.3847/2041-8213/abf1e4)

## License

MIT
