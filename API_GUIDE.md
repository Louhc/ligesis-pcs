# ligesis-pcs 接口与流程

本文档仅覆盖 `ligesis-pcs/` 内 DeepFold / LigeSIS / Ligero 及 dLigesis（分布式基准）相关的主要公开接口与协议步骤。工具性辅助函数（例如 `utils.rs`、Merkle 细节封装等）不单独展开。

## DeepFold

### 核心类型（`ligesis-pcs/src/deepfold/mod.rs`）
- `DeepFoldPCS`: PCS 实现体，提供 `setup/commit/open/verify` 等接口。
- `DeepFoldSRS`: SRS 参数（`max_mu`、`l0`、`rate` 等）。
- `DeepFoldProverParam` / `DeepFoldVerifierParam`: 由 SRS 派生的 prover/verifier 参数。
- `DeepFoldCommitment`: Merkle 根（`rt0`）+ 多项式大小 `mu`。
- `DeepFoldProverCommitmentAdvice`: prover 侧缓存（`f0/v0`、Merkle 树、分布式 upper_tree 等）。
- `DeepFoldProof`: 线性多项式+Merkle 证明（FRI folding 相关）。
- `DeepFoldBatchedProof`: 带 SumCheck 的 batch 证明（多项式联合开证明）。
- `DeepFoldExtProof` / `DeepFoldExtBatchedProof`: 扩域点（128-bit soundness）版本证明结构。
- `DeepFoldBatchCommitment` / `DeepFoldBatchProverAdvice` / `DeepFoldBatchProof`: 列哈希批量承诺/开证明结构。
- `DeepFoldBatchMultiCommitment` / `DeepFoldBatchMultiProverAdvice` / `DeepFoldBatchMultiProof`: 分块批量承诺/证明结构。
- `MultiChunkedBatchProof` / `MultiChunkedBatchExtProof`: 多承诺多点分块开证明结构（含扩域版本）。

### 基础接口（单多项式）
- `DeepFoldPCS::gen_srs_with_rate` (`deepfold/mod.rs`):
  - 原理：按 `max_mu` 和 `rate` 生成可配置冗余度的 SRS。
- `DeepFoldPCS::compute_value_from_proof` / `compute_value_from_proof_distributed`:
  - 从 `DeepFoldProof` 中提取被证明的值（分布式版本用于一致取值）。
- `deepfold_commit` (`deepfold/commit.rs`):
  - 原理：对多项式系数做 FFT 得到 RS 码字 `v0`，对 `v0` 建 Merkle 树承诺。
  - 关键步骤：evals->coeffs -> FFT -> Merkle root -> 返回 `DeepFoldCommitment` + advice。
- `deepfold_open` (`deepfold/open.rs`):
  - 原理：FRI folding 逐轮线性化，Merkle 证明一致性。
  - 关键步骤：生成 `alpha/r` 挑战 -> 逐轮折叠 -> 生成线性多项式和 Merkle proof。
- `deepfold_verify` (`deepfold/verify.rs`):
  - 原理：重放挑战，校验线性多项式一致性和 Merkle proof；做共线性检查。
  - 关键步骤：验证 `rt0`、线性关系、Merkle proof、共线性。

### 扩域点接口（128-bit）
- `deepfold_open_at_ext_point`, `deepfold_verify_at_ext_point` (`deepfold/ext.rs`):
  - 原理：挑战使用扩域点，折叠值在扩域中计算；Merkle 仍在基域。
  - 关键步骤：扩域折叠 -> 基域 Merkle 证明 -> 扩域线性关系校验。
- `deepfold_batch_open_at_ext_point`, `deepfold_batch_verify_at_ext_point`, `deepfold_d_batch_open_at_ext_point`:
  - 原理：批量 + 扩域点结合，用 SumCheck 聚合后再做扩域 folding。

### Batch（SumCheck 聚合）
- `deepfold_batch_open` / `deepfold_batch_verify` (`deepfold/open.rs`, `deepfold/verify.rs`):
  - 原理：多多项式用 SumCheck 合并到单点，再做一次 DeepFold folding。
  - 关键步骤：SumCheck 得随机点 -> 线性组合折叠 -> Merkle proof 校验。

### Batch（列哈希）
- `batch_commit` (`deepfold/chunked_batch.rs`):
  - 原理：将多个多项式 FFT 结果按列哈希，列哈希做 Merkle 树。
  - 关键步骤：对每个多项式 FFT -> 按列 hash -> Merkle root。
- `batch_open` / `batch_verify`:
  - 原理：多个多项式共享折叠挑战，FRI 一致性共用随机点。
  - 关键步骤：折叠 + Merkle proof；verify 重放挑战并校验。

### Chunked Batch（重点）
- `chunked_batch_commit` (`deepfold/chunked_batch.rs`):
  - 原理：多项式按 `base_mu` 切成 chunk，所有 chunk 统一 batch commit。
  - 关键步骤：split -> `batch_commit` -> 记录 chunk 元信息（每个多项式的 chunks 数）。
- `chunked_batch_open`:
  - 原理：对“多多项式、不同点”做 SumCheck 聚合到随机点 `r`，再对所有 chunks 在 `r_low` 开证明。
  - 关键步骤：
    1) 计算各多项式 claimed value；
    2) SumCheck: `Σ_i γ^i f_i(x) eq(point_i, x)`；
    3) 取 `r_low`（前 `base_mu` 维）评估所有 chunk；
    4) `batch_open` 生成证明。
- `chunked_batch_verify`:
  - 原理：验证 SumCheck，按 `f(x_low,x_high)=Σ_b f_b(x_low)*eq(b,x_high)` 合并 chunk。
  - 关键步骤：
    1) 验 SumCheck 得到 `r`；
    2) 用 chunk 值合成 `f_i(r)` 并校验期望值；
    3) `batch_verify` 验证所有 chunk 的开证明。
- `compute_claimed_values_from_proof`: 从证明中提取 claimed values（轻量辅助）。

### 分布式（distributed）
- `deepfold_d_commit` (`deepfold/commit.rs`):
  - 原理：各方持有本地 evals，master 汇总并构建 Merkle 上层树。
  - 关键步骤：gather evals -> master FFT -> 分发叶子 hash -> 分布式 Merkle。
- `deepfold_d_commit_v2`, `deepfold_d_commit_full_poly_v2`, `deepfold_batch_d_commit_v2`:
  - 原理：按列分发 RS 码字，列哈希作为 leaf，分布式 Merkle 构建。
  - 关键步骤：master FFT -> 按列切分 -> 各方 hash 列 -> 合成全局 Merkle。
- `deepfold_d_open`, `deepfold_d_batch_open` (`deepfold/open.rs`):
  - 原理：master 生成挑战并广播；分布式 Merkle 证明通过 `d_prove` 聚合上下层树。
  - 关键步骤：广播挑战 -> 分布式 Merkle proof -> master 汇总证明。
- `d_chunked_batch_commit`, `d_chunked_batch_open` (`deepfold/chunked_batch.rs`):
  - 原理：分块 + 分布式列哈希。小多项式先 gather 到 master，大多项式本地 split。
  - 关键步骤：
    1) split/gather；
    2) 本地 FFT -> master 汇总并重排；
    3) master 分发列 -> 各方 hash -> 分布式 Merkle；
    4) open 阶段做分布式 SumCheck + 批量开证明。
- `multi_chunked_batch_open`, `multi_chunked_batch_verify`:
  - 原理：对多个 `DeepFoldBatchMultiCommitment` 的多点开证明，用一次 SumCheck 聚合。
  - 关键步骤：重构多项式 -> SumCheck 得 `r` -> 合并 chunk -> FRI folding -> mt0 一致性检查。
- `d_multi_chunked_batch_open`:
  - 原理：多承诺 + 分布式 SumCheck；各方协同完成 Merkle 证明。
- `multi_chunked_batch_open_at_ext_point`, `d_multi_chunked_batch_open_at_ext_point`, `multi_chunked_batch_verify_at_ext_point`:
  - 原理：与上面相同，但 evaluation 点在扩域；用于 128-bit soundness。

## LigeSIS

### 核心类型（`ligesis-pcs/src/ligesis/*.rs`）
- `LigeSISPCS`: LigeSIS PCS 实现体。
- `LigeSISSRS`: 参数生成（含 `mat_a`、DeepFold SRS 等）。
- `LigeSISProverParam` / `LigeSISVerifierParam`: prover/verifier 参数。
- `LigeSISCommitment`: 对 `mat_h` 的 chunked batch 承诺。
- `LigeSISProverCommitmentAdvice`: `mat_f_prime`、`mat_h` 与其 DeepFold advice。
- `ExtSumCheckWithReductionProof`: 扩域 SumCheck 证明封装。
- `LigeSISProof`: 包含多个扩域 SumCheck 证明 + DeepFold multi-chunked batch 证明。

### Setup 与参数
- `LigeSISSRS::gen_with_params` (`ligesis/mod.rs`):
  - 原理：生成 SIS 矩阵 `mat_a` + DeepFold SRS（支持 `base_mu` 与 `code_rate`）。
- `LigeSISPCS::d_setup`:
  - 原理：分布式环境下预提交 `mat_a`（`d_chunked_batch_commit`）。
  - 关键步骤：DeepFold setup -> `mat_a` padding -> 分布式提交 -> master 生成 verifier 参数。

### Commit / Open / Verify
- `ligesis_commit` (`ligesis/commit.rs`):
  - 原理：对多项式做 RS 编码并计算 SIS Hash `H`，用 DeepFold chunked batch 承诺 `H`。
  - 关键步骤：pad -> RS encode -> SIS hash -> `chunked_batch_commit`。
- `ligesis_d_commit`:
  - 原理：各方对本地块计算 SIS hash；master 汇总 `H` 并分发；DeepFold 分布式提交。
  - 关键步骤：本地 RS+hash -> gather -> master 合并 -> `d_chunked_batch_commit`。
- `ligesis_open` (`ligesis/open.rs`):
  - 原理：扩域 SumCheck 保证编码/一致性，随后用 DeepFold 多点扩域开证明。
  - 关键步骤：
    1) 构造 `a`, `bI`, `rs_a` 并合并承诺；
    2) 多个扩域 SumCheck（`bI`、`rs_a`、`mat_g` 等）；
    3) `multi_chunked_batch_open_at_ext_point` 打开所需点。
- `ligesis_d_open`:
  - 原理：分布式 SumCheck + 分布式扩域多点开证明。
  - 关键步骤：master 广播挑战 -> 分布式 SumCheck -> `d_multi_chunked_batch_open_at_ext_point`。
- `ligesis_verify` (`ligesis/verify.rs`):
  - 原理：验证多个扩域 SumCheck 子声明，并复核 DeepFold 批量开证明。
  - 关键步骤：重放挑战 -> 校验 SumCheck 子声明 -> `multi_chunked_batch_verify_at_ext_point`。
- `LigeSISPCS::compute_value_from_proof`:
  - 从 `LigeSISProof` 中提取 claimed value（使用扩域证明里的值）。

## Ligero

### 核心类型（`ligesis-pcs/src/ligero/mod.rs`）
- `LigeroPCS`: Ligero PCS 实现体。
- `LigeroCommitment`: Merkle 根 + `num_vars`。
- `LigeroProof`: `f0/f1` 向量、被抽样列与 Merkle 证明。

### Commit / Open / Verify
- `commit`:
  - 原理：将多项式 reshape 为矩阵并 RS 编码，每列 hash 后构建 Merkle 树承诺。
  - 关键步骤：pad -> RS encode -> 列哈希 -> Merkle root。
- `open`:
  - 原理：采样随机向量 `r`，生成 `f0/f1`，抽样列并提供 Merkle proof。
  - 关键步骤：生成 `r` -> 计算 `f0/f1` -> 抽样列 -> Merkle proof。
- `verify`:
  - 原理：检查 `Enc(f)` 与抽样列一致性，并验证 Merkle proof。
  - 关键步骤：重放挑战 -> 校验 RS 关系 -> Merkle 验证。

## dLigesis（分布式基准）

`ligesis-pcs/dTests/dLigesis.rs`

- `test_multi`:
  - 原理：分布式端到端基准（`d_commit` / `d_open` / `verify`）。
  - 关键步骤：
    1) master 生成 `LigeSISSRS` 并广播；
    2) 所有节点 `LigeSISPCS::d_setup`；
    3) 每轮生成随机多项式与点 -> `d_commit` -> `d_open`；
    4) master 侧 `verify` 并统计时间/通信/证明大小。
- `main`: 通过 `common::network_run` 解析参数并启动 `test_multi`。
