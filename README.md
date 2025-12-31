# LigeSIS-PCS

A Rust implementation of **LigeSIS** (Ligero + SIS) Polynomial Commitment Scheme, based on the [HyperPlonk](https://github.com/EspressoSystems/hyperplonk) framework.

## Project Structure

```
├── ligesis-pcs/          # Main PCS library
│   ├── src/
│   │   ├── ligesis/      # LigeSIS PCS implementation
│   │   ├── deepfold/     # DeepFold PCS (used by LigeSIS)
│   │   ├── ligero/       # Ligero PCS
│   │   ├── sumcheck/     # SumCheck protocol
│   │   ├── hash/         # Merkle tree utilities
│   │   └── lib.rs
│   └── dTests/           # Distributed testing
├── arithmetic/           # Field arithmetic utilities
├── transcript/           # Fiat-Shamir transcript
├── deNetwork/            # Distributed networking
└── util/                 # General utilities
```

## Building

Requires **Rust Nightly**.

```bash
cargo build --release
```

For optimal performance:
```bash
RUSTFLAGS='-C target-cpu=native -C target-feature=+bmi2,+adx' cargo build --release
```

## Testing

Run all tests:
```bash
cargo test
```

Run tests with timing output:
```bash
cargo test -p ligesis-pcs --features print-trace -- --nocapture
```

Run specific PCS tests:
```bash
cargo test -p ligesis-pcs test_ligesis_pcs
cargo test -p ligesis-pcs test_deepfold_pcs
cargo test -p ligesis-pcs test_ligero_pcs
```

## Benchmarking

Run the LigeSIS PCS benchmark:
```bash
cargo bench --package ligesis-pcs --bench ligesis_bench --features print-trace
```

The benchmark parameters can be adjusted in `ligesis-pcs/benches/ligesis_bench.rs`:
- `MU`: Number of polynomial variables (default: 18)
- `ITERATIONS`: Number of iterations per operation (default: 1)

For optimal benchmark performance:
```bash
RUSTFLAGS='-C target-cpu=native -C target-feature=+bmi2,+adx' cargo bench --package ligesis-pcs --bench ligesis_bench
```

## Distributed Testing

Run distributed LigeSIS test with 4 parties:
```bash
cd ligesis-pcs
./dTests/run.sh dLigesis
```

The distributed test configuration is in `ligesis-pcs/dTests/data/4`.

## Features

- `print-trace`: Enable timing output for performance profiling

## License

MIT License
