# LigeSIS-PCS

A Rust implementation of **Polynomial Commitment Schemes (PCS)** for multilinear polynomials, based on the [HyperPlonk](https://github.com/EspressoSystems/hyperplonk) framework.

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

### LigeSIS Benchmark

```bash
cargo bench --package ligesis-pcs --bench ligesis_bench --features print-trace
```

Options:
- `-m, --mu <MU>`: Number of polynomial variables (default: 24)
- `-i, --iterations <N>`: Number of iterations per operation (default: 1)

```bash
cargo bench --package ligesis-pcs --bench ligesis_bench -- --mu 20 --iterations 3
```

### DeepFold Benchmark

Single polynomial benchmark:
```bash
cargo bench --package ligesis-pcs --bench deepfold_bench --features print-trace -- --mu 20
```

Batch open/verify benchmark:
```bash
cargo bench --package ligesis-pcs --bench deepfold_bench -- --test-batch --num-polys 5 --mu 18
```

Options:
- `-m, --mu <MU>`: Number of polynomial variables (default: 20)
- `-i, --iterations <N>`: Number of iterations per operation (default: 1)
- `--test-batch`: Run batch open/verify benchmark instead of single
- `-n, --num-polys <N>`: Number of polynomials for batch benchmark (default: 3)

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
