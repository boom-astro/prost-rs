# Installation

## Python Package

### From source (recommended during development)

```bash
# Clone the repository
git clone https://github.com/mcoughli/prost-rs.git
cd prost-rs

# Build and install the Rust extension
pip install ./python

# Install the Python package
pip install -e .
```

### Requirements

- Python >= 3.9
- Rust toolchain (install via [rustup](https://rustup.rs/))
- NumPy >= 1.20

### Optional dependencies

For working with FITS cutouts (image-based workflow):

```bash
pip install astropy requests
```

For development and testing:

```bash
pip install pytest astropy requests
```

## Rust Library Only

If you only need the Rust crate (e.g., for integration into another Rust project):

```toml
# In your Cargo.toml
[dependencies]
prost-rs = { git = "https://github.com/mcoughli/prost-rs.git" }
```

To enable HTTP-based catalog queries (reqwest + tokio):

```toml
prost-rs = { git = "https://github.com/mcoughli/prost-rs.git", features = ["catalogs"] }
```

## HPC Clusters

On shared HPC systems without root access:

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# Build in a virtual environment
python -m venv ~/.venvs/prost
source ~/.venvs/prost/bin/activate
pip install maturin numpy
cd prost-rs
maturin develop --release -m python/Cargo.toml
pip install -e .
```

## Building the Documentation

```bash
pip install mkdocs-material pymdownx-extensions
mkdocs serve     # local preview at http://127.0.0.1:8000
mkdocs build     # static site in site/
```

## Verifying the Installation

```python
import prost

# Catalog-based
t = prost.Transient(ra=180.0, dec=45.0)
print(t)  # Transient(ra=180.000000, dec=45.000000, name=None)

# Image-based
c = prost.Cutout.zeros(64, 64, 0.262, 180.0, 45.0)
print(c)  # Cutout(64x64, scale=0.262"/px, center=(180.0000, 45.0000))
```

## Troubleshooting

### `ImportError: cannot import name 'prost_extension'`

The Rust extension was not compiled or is not on the Python path.
Ensure you ran both `pip install ./python` and `pip install -e .`,
or `maturin develop -m python/Cargo.toml` on HPC.

### `error[E0463]: can't find crate for 'prost_rs'`

The Python extension depends on the main crate via a relative path.
Make sure you're building from the repository root, not from `python/`.

### Slow build on first compile

The first build compiles all Rust dependencies (~2-3 minutes).
Subsequent builds only recompile changed files and are much faster.

### `maturin` not finding Python

On HPC systems, ensure the virtual environment is activated before
running `maturin develop`. Maturin uses the active Python interpreter.
