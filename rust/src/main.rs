fn main() {
    println!(
        "Rust visualizer is exposed through the Python package.\n\
         From rust/:\n\
         1. uv run --with maturin maturin develop\n\
         2. uv run --with matplotlib python -m steering_functions_rust.visualize"
    );
}
