# LigeSIS Distributed Benchmark (mu=28, 4 parties)

## 总览

| 阶段 | P0 (Master) | P1 (Worker) | 通信量 (Master) |
|------|-------------|-------------|-----------------|
| Setup | 2.11s | 1.46s | sent=128MB, recv=216MB |
| Commit | 21.16s | 20.51s | sent=288MB, recv=432MB |
| Open | 16.92s | 17.57s | sent=766MB, recv=599MB |
| Verify | 597ms | - | - |
| **Total** | **41.06s** | **40.46s** | **2.43 GB** |

## Commit 阶段详情

| 操作 | P0 | P1 | 说明 |
|------|-----|-----|------|
| RSEncode (256x262144) | 11.21s | 13.48s | P1 慢 20% (云性能波动) |
| SISHash (8x256x524288) | 2.44s | 2.97s | |
| GatherMatH (33MB) | 3.99s | 13ms | P0 等待其他 party |
| AssembleAndDistribute | 91ms | 1.24s | P1 等待 master 分发 |
| **DeepFoldChunkedBatch** | **3.27s** | **2.58s** | |
| - FFT | 805ms | 971ms | |
| - Gather | 932ms | 38ms | |
| - DistColData | 200ms | 476ms | |
| - ColHash | 268ms | 241ms | |
| - Merkle | 622ms | 152ms | |

## Open 阶段详情

| 操作 | P0 | P1 | 说明 |
|------|-----|-----|------|
| ComputeA | 579ms | 1.23s | |
| ComputeBI | 792ms | 847ms | |
| ComputeRSA | 60ms | 6ms | |
| **CombinedCommit** | **5.68s** | **4.91s** | a, bI, rs_a 合并提交 |
| - FFT | 2.51s | 1.84s | bI 的 FFT (2^21 元素) |
| - Gather | 909ms | 52ms | |
| - DistColData | 424ms | 980ms | |
| - ColHash | 495ms | 354ms | |
| - Merkle | 893ms | 150ms | |
| ExtSumchecks | 2.29s | - | Master only |
| **DeepFold** | **7.32s** | **7.32s** | |
| - ComputeClaims | 986ms | 1.05s | |
| - Sumcheck | 2.06s | 2.06s | 分布式 sumcheck |
| - EvalAndCombine | 2.30s | - | Master only |
| - ComputeCombinedV0 | 1.12s | 1.27s | |
| - Folding | 760ms | 622ms | |

## 瓶颈分析

| 瓶颈 | 耗时 | 原因 |
|------|------|------|
| RSEncode | ~12s | 计算密集，服务器性能差异大 |
| GatherMatH | ~4s | 同步等待最慢 party |
| FFT (bI) | ~2.5s | bI 多项式 FFT (2^21 元素) |
| ExtSumchecks | ~2.3s | Master only 计算 |
| DeepFold.Sumcheck | ~2s | 分布式 sumcheck 通信 |
| DeepFold.EvalAndCombine | ~2.3s | Master only 计算 |

## 优化记录

- **mat_a commit 移到 Setup**: Open 阶段减少 ~1.3s，通信减少 ~344MB
