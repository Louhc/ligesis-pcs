# Code Review - Potential Issues Found

## 1. Unused `rate` field in DeepFoldSRS (Minor)

**Location:** `ligesis-pcs/src/deepfold/mod.rs:66-68`

```rust
/// Code rate multiplier (e.g., 4 for 1/4 rate, 8 for 1/8 rate)
/// len_l0 = (1 << max_mu) * rate
pub rate: usize,
```

**Issue:** The `rate` field is stored but never read anywhere in the codebase. The rate is implicitly encoded in `l0.size()` (which equals `(1 << max_mu) * rate`), so all operations use `l0.size()` directly.

**Impact:** None - this is just an unused field, not a bug.

**Recommendation:** Either:
- Remove the field if not needed
- Or add a getter method if explicit rate access is desired in the future


## 2. Inconsistent `max_mu` default in `gen_srs_for_testing` (Potential Issue)

**Location:** `ligesis-pcs/src/ligesis/mod.rs:452-455`

```rust
let deepfold_srs = DeepFoldPCS::<F>::gen_srs_for_testing(
    rng,
    log_m + 7,  // <-- Should be log_m + 9 for consistency
)?;
```

**Context:**
- `mat_a` size = `c * eta * 2^log_m` = `8 * 64 * 2^log_m` = `2^(3+6+log_m)` = `2^(log_m + 9)`
- `gen_srs_for_testing` uses `log_m + 7` for max_mu
- `gen_with_params` (distributed) defaults to `log_m + 9`

**Impact:**
- In `gen_srs_for_testing`, max_mu = log_m + 7, meaning:
  - max_mu handles polynomials up to 2^(log_m + 7) = 128 * 2^log_m
  - But mat_a has 512 * 2^log_m elements
  - This requires chunking (4 chunks)
- This should still work because the chunked batch commit handles polynomials larger than max_mu via chunking.
- However, it's inconsistent with the distributed code which uses log_m + 9.

**Recommendation:** Consider updating to `log_m + 9` for consistency:
```rust
let deepfold_srs = DeepFoldPCS::<F>::gen_srs_for_testing(
    rng,
    log_m + 9,  // Match mat_a size
)?;
```


## 3. Double broadcast in d_multi_chunked_batch_open_at_ext_point (Minor redundancy)

**Location:** `ligesis-pcs/src/ligesis/open.rs:618-638` and `chunked_batch.rs:4003-4022`

The `points_ext`, `point_to_commit`, and `point_to_poly` are broadcast twice:
1. First in `ligesis_d_open` (lines 618-638)
2. Again in `d_multi_chunked_batch_open_at_ext_point` (lines 4003-4022)

**Impact:** Minor performance overhead from redundant network operations.

**Recommendation:** Remove the first broadcast since the function already does it internally.


## 4. No issues found in core logic

The following components were reviewed and found to be correct:
- Distributed commit (`ligesis_d_commit`)
- Distributed open (`ligesis_d_open`)
- Extension field sumcheck (`run_ext_sumcheck`, `verify_ext_sumcheck_with_subclaim`)
- Verification logic (`ligesis_verify`)
- Chunked batch commit/open (`d_chunked_batch_commit`, `d_multi_chunked_batch_open_at_ext_point`)
- mat_a resize fix (correctly uses `log_m + 9` in distributed setup)


## Summary

| Issue | Severity | Status |
|-------|----------|--------|
| Unused `rate` field | Minor | Can be cleaned up |
| Inconsistent max_mu default | Low | Works via chunking, but should be consistent |
| Double broadcast | Minor | Redundant but harmless |
| Core protocol logic | N/A | Correct |
