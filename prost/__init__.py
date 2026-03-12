"""
prost — Probabilistic host galaxy association for astronomical transients.

A Rust-accelerated Python package for identifying host galaxies of transient
events using directional light radius (DLR) ranking and Bayesian posterior
scoring.

Quick start::

    from prost import Transient, GalaxyCandidate, associate_host

    transient = Transient(ra=197.4504, dec=-23.3815, redshift=0.00978)
    galaxy = GalaxyCandidate(
        ra=197.4487, dec=-23.3839,
        a_arcsec=22.0, b_arcsec=16.0, pa_deg=75.0,
        redshift=0.00978, objname="NGC 4993",
    )
    result = associate_host(transient, [galaxy])
    print(result.best_host)
"""

from ._prost import (
    Transient,
    GalaxyCandidate,
    HostCandidate,
    AssociationConfig,
    AssociationResult,
    associate_host,
    compute_dlr,
    offset_likelihood,
    redshift_likelihood,
    sersic_profile,
    elliptical_radius,
    # Image-based analysis
    Cutout,
    ExtractionConfig,
    DetectedSource,
    Background,
    FitConfig,
    SersicFit,
    extract_sources,
    estimate_background,
    fit_sersic,
)

__all__ = [
    # Catalog-based association
    "Transient",
    "GalaxyCandidate",
    "HostCandidate",
    "AssociationConfig",
    "AssociationResult",
    "associate_host",
    "compute_dlr",
    "offset_likelihood",
    "redshift_likelihood",
    "sersic_profile",
    "elliptical_radius",
    # Image-based analysis
    "Cutout",
    "ExtractionConfig",
    "DetectedSource",
    "Background",
    "FitConfig",
    "SersicFit",
    "extract_sources",
    "estimate_background",
    "fit_sersic",
]
