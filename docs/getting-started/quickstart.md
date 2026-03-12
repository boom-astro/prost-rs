# Quick Start

## Catalog-Based Association

### Basic example

```python
from prost import Transient, GalaxyCandidate, associate_host

# Define a transient
transient = Transient(ra=180.0, dec=45.0)

# Define candidate galaxies
galaxies = [
    GalaxyCandidate(
        ra=180.0 + 2.0/3600, dec=45.0,
        a_arcsec=5.0, b_arcsec=3.0, pa_deg=30.0,
        objname="Galaxy A",
    ),
    GalaxyCandidate(
        ra=180.0, dec=45.0 + 8.0/3600,
        a_arcsec=3.0, b_arcsec=2.0, pa_deg=0.0,
        objname="Galaxy B",
    ),
]

# Run association
result = associate_host(transient, galaxies)

# Inspect results
print(f"Best host: {result.best_host.galaxy.objname}")
print(f"Posterior: {result.best_host.posterior:.3f}")
print(f"Separation: {result.best_host.separation_arcsec:.1f}\"")
print(f"Fractional offset: {result.best_host.fractional_offset:.2f}")
print(f"P(no host): {result.p_none:.4f}")

for c in result.candidates:
    print(f"  rank={c.dlr_rank} {c.galaxy.objname}: "
          f"FO={c.fractional_offset:.2f}, P={c.posterior:.3f}")
```

### Using Tractor shape parameters

When working with Legacy Survey data, galaxies are parameterized with
Tractor shape parameters (`shape_r`, `shape_e1`, `shape_e2`) instead of
explicit semi-axes. Use `GalaxyCandidate.from_tractor()`:

```python
galaxy = GalaxyCandidate.from_tractor(
    ra=180.0, dec=45.0,
    shape_r=2.5, shape_e1=0.3, shape_e2=0.1,
    objname="LS source 12345",
)
print(f"a={galaxy.a_arcsec:.2f}\", b={galaxy.b_arcsec:.2f}\", PA={galaxy.pa_deg:.1f}°")
```

### With redshift information

When both transient and galaxy have redshifts, the association uses a
Gaussian redshift likelihood to break degeneracies:

```python
transient = Transient(ra=180.0, dec=45.0, redshift=0.05, redshift_err=0.001)

g1 = GalaxyCandidate(
    ra=180.0 + 3.0/3600, dec=45.0,
    a_arcsec=5.0, b_arcsec=3.0, pa_deg=0.0,
    redshift=0.05, redshift_err=0.002,  # matching redshift
    objname="match_z",
)
g2 = GalaxyCandidate(
    ra=180.0 + 2.5/3600, dec=45.0,  # closer, but wrong z
    a_arcsec=5.0, b_arcsec=3.0, pa_deg=0.0,
    redshift=0.5, redshift_err=0.002,
    objname="wrong_z",
)

result = associate_host(transient, [g1, g2])
# The redshift-matching galaxy wins even though it's slightly farther
print(result.best_host.galaxy.objname)  # "match_z"
```

### Configuration

```python
from prost import AssociationConfig

config = AssociationConfig(
    max_fractional_offset=10.0,  # Maximum FO to consider
    min_b_arcsec=0.05,           # Floor for tiny galaxies
    max_candidates=10,           # Max results returned
    use_redshift=True,           # Use redshift likelihood
)

result = associate_host(transient, galaxies, config)
```

### Low-level functions

```python
from prost import compute_dlr, offset_likelihood

# Compute DLR for a single pair
dlr = compute_dlr(
    transient_ra=180.0, transient_dec=45.0,
    galaxy_ra=180.001, galaxy_dec=45.001,
    a_arcsec=5.0, b_arcsec=3.0, pa_deg=30.0,
)
print(f"Separation: {dlr['separation_arcsec']:.1f}\"")
print(f"Directional radius: {dlr['directional_radius']:.2f}\"")
print(f"Fractional offset: {dlr['fractional_offset']:.2f}")

# Evaluate the offset likelihood
import numpy as np
fo = np.linspace(0.01, 10, 100)
likelihood = offset_likelihood(fo)  # vectorized
```

---

## Image-Based Association

When no catalog is available — or catalog shapes are unreliable — you can
extract galaxies directly from a FITS cutout.

### Load a cutout

```python
import numpy as np
from prost import Cutout

# From a FITS file (via astropy)
from astropy.io import fits
hdu = fits.open("cutout.fits")
data = hdu[0].data.astype(np.float64)

cutout = Cutout(
    data=data,
    pixel_scale=0.262,     # arcsec/pixel (Legacy Survey)
    center_ra=185.4788,
    center_dec=4.4700,
    survey="legacy",
)
print(cutout)
# Cutout(256x256, scale=0.262"/px, center=(185.4788, 4.4700))
```

### Estimate background

```python
from prost import estimate_background, ExtractionConfig

config = ExtractionConfig(back_size=64)
bg = estimate_background(cutout, config)
print(f"Background: mean={bg.global_mean:.4f}, rms={bg.global_rms:.4f}")
# bg.level and bg.rms are 2D numpy arrays
```

### Extract sources

```python
from prost import extract_sources

sources = extract_sources(cutout, config)
print(f"Detected {len(sources)} sources")

for s in sources[:3]:
    print(f"  ({s.ra:.4f}, {s.dec:.4f}) flux={s.flux:.0f} "
          f"SNR={s.snr:.1f} a={s.a_pix:.1f}px")
```

### Fit Sersic profiles

```python
from prost import fit_sersic, FitConfig

fit_cfg = FitConfig(max_sersic_n=8.0)

for src in sources:
    if src.a_pix < 2.0:
        continue  # skip unresolved (likely stars)

    try:
        fit = fit_sersic(cutout, src, fit_cfg)
    except RuntimeError:
        continue

    if fit.converged:
        print(f"  n={fit.n:.2f}, r_eff={fit.r_eff:.2f}\", "
              f"q={fit.axis_ratio:.2f}, chi2r={fit.chi2_reduced:.2f}")
```

### Build candidates and associate

```python
from prost import GalaxyCandidate, Transient, associate_host

galaxy_candidates = []
for src in sources:
    if src.a_pix < 2.0:
        continue

    try:
        fit = fit_sersic(cutout, src, fit_cfg)
    except RuntimeError:
        continue

    if not fit.converged:
        continue

    a = fit.r_eff / fit.axis_ratio
    b = fit.r_eff
    if a < b:
        a, b = b, a

    galaxy_candidates.append(GalaxyCandidate(
        ra=fit.ra, dec=fit.dec,
        a_arcsec=a, b_arcsec=b, pa_deg=fit.pa_deg,
        shape_from_image=True,
    ))

transient = Transient(ra=185.4788, dec=4.4700, name="SN 2020jfo")
result = associate_host(transient, galaxy_candidates)

if result.best_host:
    print(f"Best host: sep={result.best_host.separation_arcsec:.1f}\", "
          f"FO={result.best_host.fractional_offset:.2f}, "
          f"P={result.best_host.posterior:.3f}")
```

For a complete walkthrough with a real Legacy Survey FITS download, see
[Image-Based Association](../examples/image-based-association.md).

### Synthetic test image

You can also build test images programmatically:

```python
cutout = Cutout.zeros(101, 101, 0.262, 180.0, 45.0)
cutout.add_sersic(
    center_row=50.0, center_col=50.0,
    i_eff=200.0, r_eff_pix=10.0,
    n=1.0,           # exponential disk
    axis_ratio=0.6,
    pa_rad=0.5,
)

sources = extract_sources(cutout)
fit = fit_sersic(cutout, sources[0])
print(f"Recovered: n={fit.n:.2f}, r_eff={fit.r_eff_pix:.1f}px, q={fit.axis_ratio:.2f}")
```
