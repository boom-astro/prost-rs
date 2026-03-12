// Catalog-based association (offline equivalent of BOOM's src/utils/host.rs)
pub mod associate;
pub mod dlr;
pub mod ellipse;
pub mod errors;
pub mod likelihood;
pub mod prior;
pub mod types;

// Image-based analysis (not in BOOM — cutouts, source extraction, morphology)
pub mod cutout;
pub mod morphology;
pub mod source;

pub use associate::{associate_host, AssociationConfig, AssociationResult};
pub use types::{GalaxyCandidate, HostCandidate, Transient};
