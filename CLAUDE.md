# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

LigeSIS-PCS is a Rust implementation of Polynomial Commitment Schemes (PCS) for multilinear polynomials, based on the HyperPlonk framework. It includes three PCS implementations: LigeSIS, DeepFold, and Ligero.

## Build Commands

```bash
# Build all packages
cargo build --release

# Build with timing output (for performance profiling)
cargo build --release --features print-trace
```

## Testing

```bash
# Run all tests
cargo test

# Run specific PCS tests
cargo test -p ligesis-pcs test_ligesis_pcs
cargo test -p ligesis-pcs test_deepfold_pcs
cargo test -p ligesis-pcs test_ligero_pcs
```

## Benchmarking

```bash
# LigeSIS benchmark
cargo bench --package ligesis-pcs --bench ligesis_bench --features print-trace -- --mu 20

# DeepFold single polynomial benchmark
cargo bench --package ligesis-pcs --bench deepfold_bench --features print-trace -- --mu 20

# DeepFold batch benchmark
cargo bench --package ligesis-pcs --bench deepfold_bench -- --test-batch --num-polys 5 --mu 18
```

## Distributed Testing

```bash
# Local distributed testing (4 parties default)
cd ligesis-pcs/dTests
python3 run.py dLigesis -m 20                      # Basic test
python3 run.py dLigesis -n 8 -m 24 --trace         # 8 parties, mu=24, with timing

# Remote distributed testing (multi-server)
python3 run.py dLigesis --servers servers.json --sync --build -m 28  # Full deploy
python3 run.py dLigesis --servers servers.json -m 28                  # Run only
```

Available distributed test examples: `dLigesis`, `dDeepFold`, `dDeepFoldBatch`, `dMerkle`, `dChunkedBatch`, `dMultiChunkedBatchBench`, `dMultiChunkedBatchProfile`

## Architecture

### Workspace Structure

- `ligesis-pcs/` - Main PCS implementations (LigeSIS, DeepFold, Ligero)
- `arithmetic/` - Polynomial arithmetic (multilinear extensions, virtual polynomials, univariate)
- `deNetwork/` - Distributed network layer for MPC-style protocols
- `transcript/` - Fiat-Shamir transcript (Merlin-based)
- `util/` - Timing utilities
- `third_party/` - Patched arkworks algebra libraries

### Core Modules (ligesis-pcs/src/)

**PCS Implementations:**
- `ligesis/` - LigeSIS PCS using SIS hash + RS encoding + DeepFold
- `deepfold/` - DeepFold PCS (FRI-style folding with Merkle commitments)
- `ligero/` - Ligero PCS

**Key DeepFold submodules:**
- `commit.rs` - Polynomial commitment (FFT-based with Merkle trees)
- `open.rs` - Single/batch opening proofs
- `ext.rs` - Extension field support (128-bit soundness)
- `multi.rs` - Multi-polynomial chunked commits
- `chunked_batch.rs` - Distributed batch commit/open with chunking

**Supporting modules:**
- `sumcheck/` - Interactive sumcheck protocol with distributed support (`d_prove`)
- `ext_sumcheck.rs` - Extension field sumcheck
- `hash/` - SHA256-based Merkle trees

### Distributed Protocol Pattern

The codebase uses `deNetwork::DeMultiNet` for MPC-style distribution:

```rust
use deNetwork::{DeMultiNet as Net, DeNet};

// Check if running distributed
if Net::is_init() {
    // Master-worker pattern
    if Net::am_master() {
        let gathered = Net::send_to_master(&local_data);  // Master gathers
        let result = compute(gathered.unwrap());
        Net::recv_from_master_uniform(Some(result));       // Master broadcasts
    } else {
        Net::send_to_master(&local_data);                  // Workers send
        let result = Net::recv_from_master_uniform(None);  // Workers receive
    }
}
```

### Key Types

- `DenseMultilinearExtension<F>` - Multilinear polynomial stored as evaluations over hypercube
- `VirtualPolynomial<F>` - Sum of products of MLEs (for sumcheck)
- `IOPTranscript<F>` - Fiat-Shamir transcript for challenge generation
- `DeepFoldCommitment` - Merkle root commitment
- `DeepFoldProverCommitmentAdvice` - Prover state (FFT results, Merkle tree)

### Field Types

The codebase uses Goldilocks field (`FGoldilocks`) with quadratic extension for 128-bit soundness:
- Base field: 64-bit Goldilocks prime
- Extension field: Quadratic extension for extension-field sumcheck/opening

## Features

- `print-trace` - Enable ark-std timing macros for performance profiling
- `debug` - Additional debug assertions
