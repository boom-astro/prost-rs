# Contributing

## Development Setup

```bash
git clone https://github.com/boom-astro/prost-rs.git
cd prost-rs

# Rust tests
cargo test

# Python extension
pip install maturin numpy pytest astropy
maturin develop -m python/Cargo.toml
pip install -e .
pytest
```

## Running Tests

```bash
# Unit tests (fast, no network)
cargo test

# Integration tests (require network / credentials)
cargo test --features catalogs -- --ignored

# Python tests
pytest
```

## Code Style

- Rust: `cargo fmt` and `cargo clippy`
- Python: standard PEP 8

## Adding a New Module

1. Create `src/your_module.rs`
2. Add `pub mod your_module;` to `src/lib.rs`
3. Add Python bindings in `python/src/lib.rs`
4. Re-export in `prost/_prost.py` and `prost/__init__.py`
5. Add API docs in `docs/api/index.md`
6. Add to `mkdocs.yml` nav if it warrants its own page

## Submitting Changes

1. Create a branch from `main`
2. Ensure `cargo test` and `cargo clippy` pass
3. Update docs and CHANGELOG.md
4. Open a pull request
