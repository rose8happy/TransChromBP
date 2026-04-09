# 6000 / 6002 远端前置资源盘点（2026-04-07）

## 目的

在当前 `multiscale decoder probe` 双机并行进行时，回答一个更务实的问题：

> 现在两台机器还有没有必须提前下载的数据集，或必须提前安装的环境/工具，才能不耽误下一阶段工作？

本报告只回答“对接下来 1-2 个阶段最相关的前置条件”，不把所有潜在大模型方向都提前铺开。

## 结论先行

### 1. 对当前主线实验，没有必须立刻补装的大包

当前 `6000 msdec_v1 short10` 与 `6002 msdec_v1_s2 short10` 都已经实际起跑，说明：

- 6000 当前主线训练环境足够
- 6002 当前单卡结构实验环境足够
- tutorial / GM12878 / K562 这批现有训练数据也足够

所以从“马上继续跑结构实验”这个角度看，**没有必须立刻下载的新数据集，也没有必须立刻安装的新训练环境**。

### 2. 真正的近端 blocker 只有两个

#### A. 6000 的 AlphaGenome 现在已经过了认证 smoke，阻塞已从“凭据缺失”变成“准备 matched slice”

实机检查结果：

- 6000 存在环境：`alphagenome`, `chrombpnet`, `nucleotide-transformer-py311`, `transchrombp`
- 6000 的 `alphagenome` 环境里 `alphagenome` 可 import，版本 `0.6.1`
- `2026-04-07` 新增认证检查结果：
  - `ALPHAGENOME_API_KEY = True`
  - `GOOGLE_API_KEY = False`
- 同轮真实 API smoke 已通过：
  - `dna_client.create(api_key=...)` 可成功建立 client
  - `output_metadata(Organism.HOMO_SAPIENS)` 返回正常
  - 已实际取回 `ATAC / DNase / RNA-seq` metadata

这意味着：

- **AlphaGenome SDK 已就绪**
- **当前 `ALPHAGENOME_API_KEY` 已足以支撑现阶段 SDK 调用**
- **`GOOGLE_API_KEY` 仍未配置，但对当前 `dna_client` 路径不是 blocker**

如果下一阶段要做 `AlphaGenome matched slice`，现在最该做的已经不是重装或继续折腾认证，而是开始挑 `matched loci`、确认目标 track/filter 口径，并设计第一轮窄切片对照。

#### B. 6002 缺的是“数据处理工具链”，不是训练环境

实机检查结果：

- 6002 当前只有一个环境：`/home/zhengwei/bylw_atac/.mamba/envs/transchrombp`
- 这个环境里训练基础包都在：
  - `torch=True`
  - `yaml=True`
  - `numpy=True`
  - `pandas=True`
- Python 侧 I/O 包部分可用：
  - `pyBigWig=True`
  - `pyfaidx=True`
  - `pysam=False`
- 但命令行工具链缺失：
  - `bedGraphToBigWig=None`
  - `samtools=None`
  - `bedtools=None`
  - `bgzip=None`
  - `tabix=None`

这意味着：

- 6002 **适合继续跑现成训练/评估**
- 6002 **不适合马上承担新的 raw dataset prep / bigWig 构建 / official preprocessing**

当前若要在 6002 上新做数据预处理，真正缺的是一套 `chrombpnet / bioinformatics tools` 侧环境，而不是 PyTorch 训练环境。

## 当前资源状态

### 6000

已验证：

- 环境：
  - `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet`
  - `/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp`
  - `/data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome`
  - `/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b`
- foundation 目录：
  - `alphagenome`: 存在
  - `genos`: 存在
  - `nucleotide_transformer`: 存在
  - `gengram`: 存在
- 缺失目录：
  - `segmentnt`
  - `grelu`
  - `hyenadna`

判断：

- 这些缺失目录**不是当前 `msdec / AlphaGenome matched slice` 的立即 blocker**
- 它们属于“如果后续真的转入第二层路线，再决定要不要补”的资源

### 6002

已验证：

- `transchrombp` 单卡训练环境可用
- tutorial / GM12878 / K562 数据目录存在
- 当前单卡 `msdec_v1_s2` 已实际跑起

缺口：

- 缺少 official/bioinformatics CLI 工具链
- 缺少 `pysam`

判断：

- 对当前结构实验无影响
- 对未来“6002 自己做新数据预处理”会构成硬门槛

## 建议优先级

### P0: 现在不建议动的事

1. 不建议在 6002 当前 `transchrombp` 环境上边跑训练边直接装新包
   - 6002 现在正跑 `msdec_v1_s2`
   - 直接改 active env 有不必要风险

2. 不建议现在就去大包下载 `SegmentNT / gReLU / HyenaDNA`
   - 当前主线还在 `multiscale decoder`
   - 这些下载既不解决当前 blocker，也会制造新的管理负担

### P1: 最值得提前准备的两件事

1. **准备 6000 的 AlphaGenome matched slice 任务定义**
   - 认证和 metadata smoke 已通过
   - 下一步应转向 `matched loci` 抽样、目标 track 筛选和最小 raw-track 对照脚本

2. **在 6002 准备独立的数据处理工具环境**
   - 不改当前正在用的 `transchrombp` env
   - 建议新建单独工具环境，至少包含：
     - `samtools`
     - `bedtools`
     - `htslib` / `bgzip` / `tabix`
     - `ucsc-bedgraphtobigwig`
     - `pysam`

### P2: 只有路线确认后才值得做的预下载

若后续真的要转向“参考外部多尺度 / segmentation / sequence-to-track 系统”，再考虑补：

1. `SegmentNT`
2. `gReLU / Borzoi / human-atac-catlas-model`
3. `HyenaDNA`

这些目前都不是马上需要。

## 我的执行判断

基于当前状态，**我不建议现在立刻启动任何大模型权重下载或在运行中的训练环境里装包**。

更合理的动作是：

1. 继续让 6000 / 6002 当前两条 `short10 smoke` 收口
2. 并行向你确认或准备 `AlphaGenome API key`
3. 等 6002 当前 run 结束后，再决定是否为它补一个独立 `tools` 环境

## 一句话总结

**当前主线实验不缺训练环境也不缺现成数据；真正该提前准备的是 6000 的 AlphaGenome 凭据，以及 6002 的独立数据处理工具环境，而不是继续下载新的大模型或数据集。**
