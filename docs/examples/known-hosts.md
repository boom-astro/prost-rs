# Known Transient–Host Pairs

Validation against well-studied transient events with unambiguous host
galaxy associations.

## AT2017gfo in NGC 4993

The electromagnetic counterpart to gravitational wave event GW170817,
located ~10" from the center of NGC 4993 (an E/S0 galaxy at z=0.00978).

```python
from prost import Transient, GalaxyCandidate, associate_host

transient = Transient(
    ra=197.4504, dec=-23.3815,
    redshift=0.009787, redshift_err=0.00015,
    name="AT2017gfo",
)

ngc4993 = GalaxyCandidate(
    ra=197.4487, dec=-23.3839,
    a_arcsec=22.0, b_arcsec=16.0, pa_deg=75.0,
    redshift=0.009787, redshift_err=0.00015,
    mag=13.05, objname="NGC 4993",
)

# Add a decoy galaxy at similar separation but wrong redshift
decoy = GalaxyCandidate(
    ra=197.46, dec=-23.37,
    a_arcsec=3.0, b_arcsec=2.0, pa_deg=30.0,
    redshift=0.08, objname="field galaxy",
)

result = associate_host(transient, [ngc4993, decoy])

print(f"Best host: {result.best_host.galaxy.objname}")
# Best host: NGC 4993
print(f"P(NGC 4993) = {result.best_host.posterior:.3f}")
print(f"P(none) = {result.p_none:.4f}")
```

## SN 2014J in M82

One of the nearest Type Ia supernovae, discovered in the starburst
galaxy M82 at a distance of ~3.5 Mpc.

```python
transient = Transient(
    ra=148.9256, dec=69.6739,
    redshift=0.000677, name="SN 2014J",
)

m82 = GalaxyCandidate(
    ra=148.9685, dec=69.6797,
    a_arcsec=335.0, b_arcsec=155.0, pa_deg=65.0,
    redshift=0.000677, mag=8.41, objname="M82",
)

result = associate_host(transient, [m82])
print(f"FO = {result.best_host.fractional_offset:.3f}")
# SN 2014J is well within M82's disk (FO < 1)
```

## SN 2011fe in M101

A prototypical Type Ia supernova in the nearby face-on spiral M101.

```python
transient = Transient(
    ra=210.7742, dec=54.2739,
    redshift=0.000804, name="SN 2011fe",
)

m101 = GalaxyCandidate(
    ra=210.8024, dec=54.3488,
    a_arcsec=840.0, b_arcsec=840.0, pa_deg=0.0,
    redshift=0.000804, mag=7.86, objname="M101",
)

result = associate_host(transient, [m101])
print(f"Separation: {result.best_host.separation_arcsec:.0f}\"")
print(f"FO = {result.best_host.fractional_offset:.3f}")
# ~280" from center, well within the ~14' disk
```
