# LigeSIS-PCS

A Rust implementation of **Polynomial Commitment Schemes (PCS)** for multilinear polynomials, based on the [HyperPlonk](https://github.com/EspressoSystems/hyperplonk) framework.

## Testing

Run all tests:
```bash
cargo test
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

### Local Mode

Run distributed LigeSIS test locally:
```bash
cd ligesis-pcs/dTests
python3 run.py dLigesis              # 4 parties (default)
python3 run.py dLigesis -n 8         # 8 parties
python3 run.py dLigesis -m 24        # mu=24
python3 run.py dLigesis --trace      # Enable internal timing
```

Options:
- `-n, --num-parties <N>`: Number of parties (default: 4, must be power of 2)
- `-m, --mu <MU>`: Number of polynomial variables (default: 20)
- `--trace`: Enable internal timing output
- `--port <PORT>`: Base port (default: 18000)

### Remote Mode (Multi-Server)

1. Create a `servers.json` config file:
```json
{
    "servers": [
        {"host": "10.128.0.2", "ssh_host": "35.202.139.171"},
        {"host": "10.128.0.3", "ssh_host": "104.197.202.243"},
        {"host": "10.128.0.4", "ssh_host": "34.72.91.60"},
        {"host": "10.128.0.5", "ssh_host": "34.69.184.100"}
    ],
    "user": "ubuntu",
    "ssh_key": "~/.ssh/id_ed25519",
    "remote_dir": "~/ligesis-pcs",
    "network_port": 18000
}
```

Config fields:
- `host`: Internal IP for inter-node communication
- `ssh_host`: (Optional) Public IP for SSH access (defaults to `host`)
- `user`: SSH username (global, can be overridden per-server)
- `ssh_key`: (Optional) Path to SSH private key
- `remote_dir`: Code location on remote servers
- `network_port`: Port for distributed protocol communication

2. Run the test:
```bash
# Sync code + build + run
python3 run.py dLigesis --servers servers.json --sync --build -m 24

# Sync and run (if already built)
python3 run.py dLigesis --servers servers.json --sync -m 24

# Run only (if code already synced and built)
python3 run.py dLigesis --servers servers.json -m 24

# With internal timing output
python3 run.py dLigesis --servers servers.json -m 24 --trace
```

Remote mode options:
- `--sync`: Sync local code to all remote servers
- `--build`: Build on remote servers before running
- `--trace`: Enable internal timing output

Requirements:
- SSH key-based authentication to all servers
- Rust toolchain installed on all servers
- Network connectivity between servers on the specified port

## Features

- `print-trace`: Enable timing output for performance profiling

## License

MIT License
