# steering-functions

This repository now contains:

- Legacy C++ implementation under `src/steering_functions`
- A new Rust workspace and crate under `rust/steering_functions`

## Rust quick start

From the repository root:

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --all-targets
```

## CI

GitHub Actions CI is configured in `.github/workflows/ci.yml` and runs Rust format, lint, and tests on every push and pull request.
