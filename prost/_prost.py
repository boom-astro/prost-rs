"""
Wrapper around the Rust extension module.

Imports the compiled prost_extension (.so) and re-exports its classes
and functions with optional Python-level convenience additions.
"""

try:
    from . import prost_extension as _ext
except ImportError:
    import prost_extension as _ext

# Re-export all classes and functions from the extension
Transient = _ext.Transient
GalaxyCandidate = _ext.GalaxyCandidate
HostCandidate = _ext.HostCandidate
AssociationConfig = _ext.AssociationConfig
AssociationResult = _ext.AssociationResult

associate_host = _ext.associate_host
compute_dlr = _ext.compute_dlr
offset_likelihood = _ext.offset_likelihood
redshift_likelihood = _ext.redshift_likelihood
sersic_profile = _ext.sersic_profile
elliptical_radius = _ext.elliptical_radius

# Image-based analysis
Cutout = _ext.Cutout
ExtractionConfig = _ext.ExtractionConfig
DetectedSource = _ext.DetectedSource
Background = _ext.Background
FitConfig = _ext.FitConfig
SersicFit = _ext.SersicFit

extract_sources = _ext.extract_sources
estimate_background = _ext.estimate_background
fit_sersic = _ext.fit_sersic

# BOOM client
BoomConfig = _ext.BoomConfig
BoomClient = _ext.BoomClient
