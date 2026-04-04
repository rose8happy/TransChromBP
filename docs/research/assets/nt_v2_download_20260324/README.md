# NT v2 下载辅助资产

本目录保存 `nucleotide-transformer-v2-500m-multi-species` 下载与最小 smoke 验证过程中沉淀下来的小型可复用资产：

- `nt_v2_500m_resolved_urls.tsv`：当时解析出的直链清单
- `start_6000_nucleotide_transformer_download.sh`：镜像/HF 下载脚本
- `start_6000_nucleotide_transformer_direct_urls.sh`：直链补下载脚本
- `smoke_test_nucleotide_transformer_hf.py`：本地文件 smoke test

不纳入本仓库的内容：

- 大模型权重、压缩包、解包样例目录、缓存锁文件
- 远端下载过程产生的原始日志全文
