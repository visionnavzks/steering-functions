# steering-functions

This repository now contains:

- Legacy C++ implementation under `/home/runner/work/steering-functions/steering-functions/src/steering_functions`
- A new Rust workspace and crate under `/home/runner/work/steering-functions/steering-functions/rust/steering_functions`

## Rust quick start

From `/home/runner/work/steering-functions/steering-functions`:

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --all-targets
```

## CI

GitHub Actions CI is configured in `/home/runner/work/steering-functions/steering-functions/.github/workflows/ci.yml` and runs Rust format, lint, and tests on every push and pull request.
