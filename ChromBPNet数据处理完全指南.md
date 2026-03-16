# 📊 ChromBPNet 数据处理完全指南

> 从原始文件到训练批次的完整流程

---

## 目录

1. [数据流程概览](#数据流程概览)
2. [输入文件格式](#输入文件格式)
3. [数据加载详解](#数据加载详解)
4. [批次生成器详解](#批次生成器详解)
5. [数据增强详解](#数据增强详解)
6. [完整代码注释](#完整代码注释)

---

## 数据流程概览

### 整体流程图

```
原始数据文件
├── peaks.bed          (峰区域)
├── nonpeaks.bed       (背景区域)
├── genome.fa          (参考基因组)
└── signal.bw          (ATAC-seq/DNase-seq信号)
         ↓
    ┌────────────────────────┐
    │  Step 1: 数据加载      │
    │  (data_utils.load_data)│
    └──────────┬─────────────┘
               ↓
    从文件中读取:
    - DNA序列 (序列矩阵)
    - 信号强度 (计数向量)
    - 坐标信息 (染色体+位置)
               ↓
    ┌────────────────────────┐
    │  Step 2: 批次生成      │
    │  (ChromBPNetBatchGen)  │
    └──────────┬─────────────┘
               ↓
    每个epoch:
    - 随机裁剪 (jittering)
    - 反向互补 (rev-comp)
    - 打乱顺序 (shuffling)
    - 负采样 (negative sampling)
               ↓
    ┌────────────────────────┐
    │  Step 3: 批次输出      │
    └────────────────────────┘
               ↓
    返回训练批次:
    - X: DNA序列 (batch, 2114, 4)
    - y: [Profile, Counts]
      - Profile: (batch, 1000)
      - Counts: (batch, 1)
```

---

## 输入文件格式

### 1. Peaks 文件 (narrowPeak格式)

**标准的10列BED文件：**

```
chr1  1000  1500  peak1  100  .  5.2  10.3  8.1  250
chr1  2000  2500  peak2  150  .  6.1  12.5  9.3  280
chr2  3000  3500  peak3  120  .  5.8  11.2  8.7  230
...
```

**列说明：**

| 列号 | 名称 | 说明 | 示例 |
|------|------|------|------|
| 1 | chr | 染色体 | chr1 |
| 2 | start | 起始位置 | 1000 |
| 3 | end | 结束位置 | 1500 |
| 4 | name | 峰名称 | peak1 |
| 5 | score | 分数 | 100 |
| 6 | strand | 链方向 | . |
| 7 | signalValue | 信号值 | 5.2 |
| 8 | pValue | p值 | 10.3 |
| 9 | qValue | q值 | 8.1 |
| 10 | **summit** | 🔑 峰值偏移 | 250 |

**关键：第10列 summit**

这是最重要的列！它指定了峰的精确中心：

```python
# 峰的实际中心位置
center = start + summit

# 例如
start = 1000
summit = 250
center = 1000 + 250 = 1250

# 提取序列范围（假设inputlen=2114）
seq_start = center - 2114//2 = 1250 - 1057 = 193
seq_end = center + 2114//2 = 1250 + 1057 = 2307
```

**可视化：**

```
染色体位置:  1000        1250          1500
            ↓           ↓             ↓
            |===========*=============|
            start     summit          end
                        ↑
                    峰的中心
                    
提取序列:            [=========2114bp=========]
                         以峰中心为中心
```

### 2. NonPeaks 文件

**格式与peaks相同，但代表背景区域（无TF结合）**

```
chr1  10000  10500  nonpeak1  0  .  0.5  1.2  0.8  250
chr1  20000  20500  nonpeak2  0  .  0.6  1.5  0.9  240
...
```

**用途：**
- 训练偏置模型（第一步）
- 作为负样本训练主模型（第二步）

### 3. 基因组文件 (FASTA)

**标准FASTA格式：**

```
>chr1
ATCGATCGATCGATCGATCGATCGATCG...
>chr2  
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
>chr3
...
```

**使用库：** `pyfaidx`

```python
genome = pyfaidx.Fasta("genome.fa")
# 提取序列
seq = genome['chr1'][1000:2000]  # 从chr1的1000到2000位置
```

### 4. BigWig 文件

**存储每个碱基的信号强度**

```
位置:    1000  1001  1002  1003  1004  ...
信号:     0.5   1.2   3.4   5.6   3.2  ...
```

**使用库：** `pyBigWig`

```python
bw = pyBigWig.open("signal.bw")
# 提取指定区域的信号
values = bw.values('chr1', 1000, 2000)  # 返回numpy数组
```

**注意：** 信号必须是 **+4/-4 shifted** 的

```
为什么需要shift？

Tn5转座酶切割DNA的位置：
   实际切割点
        ↓
5' -----*--------- 3'
3' ---------*----- 5'
        ↑
   测序得到的位置

需要shift才能得到真实的切割中心！
```

---

## 数据加载详解

### data_utils.py 核心函数

#### 函数1: get_seq()

```python
def get_seq(peaks_df, genome, width):
    """
    从基因组中提取DNA序列
    
    参数:
        peaks_df: DataFrame，包含chr, start, summit列
        genome: pyfaidx.Fasta对象
        width: 要提取的序列长度
    
    返回:
        one-hot编码的序列数组 (N, width, 4)
    
    例子:
        如果 width=2114，center=1250
        提取范围: [1250-1057, 1250+1057] = [193, 2307]
    """
    vals = []
    
    for i, r in peaks_df.iterrows():
        # 计算中心位置
        center = r['start'] + r['summit']
        
        # 计算提取范围
        start_pos = center - width//2
        end_pos = center + width//2
        
        # 提取序列
        sequence = str(genome[r['chr']][start_pos:end_pos])
        vals.append(sequence)
    
    # 转换为one-hot编码
    return one_hot.dna_to_one_hot(vals)
```

**One-hot编码详解：**

```python
# 输入: "ATCG"
# 输出:
# A: [1, 0, 0, 0]
# T: [0, 1, 0, 0]
# C: [0, 0, 1, 0]
# G: [0, 0, 0, 1]

# 最终形状: (序列长度, 4)
```

#### 函数2: get_cts()

```python
def get_cts(peaks_df, bw, width):
    """
    从BigWig文件中提取信号值
    
    参数:
        peaks_df: DataFrame
        bw: pyBigWig对象
        width: 提取宽度
    
    返回:
        信号数组 (N, width)
    
    处理:
        - 使用np.nan_to_num处理缺失值（转为0）
    """
    vals = []
    
    for i, r in peaks_df.iterrows():
        # 计算中心和范围
        center = r['start'] + r['summit']
        start_pos = center - width//2
        end_pos = center + width//2
        
        # 提取信号值
        signal = bw.values(r['chr'], start_pos, end_pos)
        
        # 处理NaN（基因组某些区域可能无数据）
        signal = np.nan_to_num(signal)
        vals.append(signal)
    
    return np.array(vals)
```

**信号示例：**

```python
# 假设width=10
signal = [0.5, 1.2, 3.4, 5.6, 8.2, 5.1, 2.3, 1.0, 0.6, 0.3]

# 这代表该区域每个碱基的ATAC-seq切割次数
# 峰值在位置4-5（值最大）
```

#### 函数3: get_coords()

```python
def get_coords(peaks_df, peaks_bool):
    """
    获取区域的坐标信息
    
    返回:
        坐标数组 [[chr, center, strand, is_peak], ...]
    
    用途:
        - 追踪每个样本的来源
        - 评估时使用
        - 可视化时使用
    """
    vals = []
    for i, r in peaks_df.iterrows():
        vals.append([
            r['chr'],              # 染色体
            r['start'] + r['summit'],  # 中心位置
            "f",                   # 正链（forward）
            peaks_bool             # 1=peak, 0=nonpeak
        ])
    return np.array(vals)
```

#### 主函数: load_data()

```python
def load_data(bed_regions, nonpeak_regions, genome_fasta, 
              cts_bw_file, inputlen, outputlen, max_jitter):
    """
    加载所有训练数据的主函数
    
    关键参数:
        inputlen: 输入序列长度（例如2114）
        outputlen: 输出预测长度（例如1000）
        max_jitter: 最大抖动范围（例如128）
    
    关键设计:
        - Peaks加载时额外加 2*max_jitter
        - 允许训练时随机裁剪（数据增强）
        
    返回:
        (peak_seqs, peak_cts, peak_coords,
         nonpeak_seqs, nonpeak_cts, nonpeak_coords)
    """
    
    # 打开文件
    cts_bw = pyBigWig.open(cts_bw_file)
    genome = pyfaidx.Fasta(genome_fasta)
    
    # 初始化
    train_peaks_seqs = None
    train_peaks_cts = None
    train_peaks_coords = None
    train_nonpeaks_seqs = None
    train_nonpeaks_cts = None
    train_nonpeaks_coords = None
    
    # 加载峰区域
    if bed_regions is not None:
        # 🔑 关键：为peaks额外加载更宽的区域
        # 这样训练时可以随机裁剪
        train_peaks_seqs, train_peaks_cts, train_peaks_coords = \
            get_seq_cts_coords(
                bed_regions,
                genome,
                cts_bw,
                inputlen + 2*max_jitter,   # ← 更宽！
                outputlen + 2*max_jitter,  # ← 更宽！
                peaks_bool=1
            )
    
    # 加载非峰区域
    if nonpeak_regions is not None:
        # 非峰区域不需要额外宽度
        # 因为信号本身就很低，不需要太多增强
        train_nonpeaks_seqs, train_nonpeaks_cts, train_nonpeaks_coords = \
            get_seq_cts_coords(
                nonpeak_regions,
                genome,
                cts_bw,
                inputlen,      # ← 精确长度
                outputlen,     # ← 精确长度
                peaks_bool=0
            )
    
    # 关闭文件
    cts_bw.close()
    genome.close()
    
    return (train_peaks_seqs, train_peaks_cts, train_peaks_coords,
            train_nonpeaks_seqs, train_nonpeaks_cts, train_nonpeaks_coords)
```

**为什么peaks需要额外宽度？**

```
假设: inputlen=2114, max_jitter=128

加载的实际宽度: 2114 + 2*128 = 2370

训练时随机裁剪:
[========2370========]
  ↓随机偏移↓
 [===2114===]  ← 可以向左偏移最多128
    [===2114===]  ← 中心
       [===2114===]  ← 可以向右偏移最多128

这样每个epoch看到的序列都略有不同！
→ 数据增强，防止过拟合
```

---

## 批次生成器详解

### ChromBPNetBatchGenerator 类

这是一个 **Keras Sequence** 类，负责：
1. 管理训练数据
2. 每个epoch进行数据增强
3. 按批次返回数据

#### 初始化

```python
class ChromBPNetBatchGenerator(keras.utils.Sequence):
    """
    批次生成器
    
    特点:
    - 继承keras.utils.Sequence（支持多进程）
    - 每个epoch重新裁剪和增强数据
    - 支持负采样比例调整
    """
    
    def __init__(self, peak_regions, nonpeak_regions, genome_fasta, 
                 batch_size, inputlen, outputlen, max_jitter, 
                 negative_sampling_ratio, cts_bw_file, 
                 add_revcomp, return_coords, shuffle_at_epoch_start):
        """
        参数说明:
            peak_regions: DataFrame，峰区域
            nonpeak_regions: DataFrame，非峰区域
            genome_fasta: 基因组文件路径
            batch_size: 批次大小（例如64）
            inputlen: 输入长度（2114）
            outputlen: 输出长度（1000）
            max_jitter: 最大抖动（128）
            negative_sampling_ratio: 负样本比例（0.0-1.0）
                - 1.0: 使用所有负样本
                - 0.1: 只用10%的负样本
            cts_bw_file: BigWig文件路径
            add_revcomp: 是否添加反向互补（True/False）
            return_coords: 是否返回坐标（True/False）
            shuffle_at_epoch_start: 是否在epoch开始时打乱（True/False）
        """
        
        # ========================================
        # Step 1: 加载所有数据
        # ========================================
        peak_seqs, peak_cts, peak_coords, \
        nonpeak_seqs, nonpeak_cts, nonpeak_coords = \
            data_utils.load_data(
                peak_regions, nonpeak_regions, genome_fasta, 
                cts_bw_file, inputlen, outputlen, max_jitter
            )
        
        # 保存到实例变量
        self.peak_seqs, self.nonpeak_seqs = peak_seqs, nonpeak_seqs
        self.peak_cts, self.nonpeak_cts = peak_cts, nonpeak_cts
        self.peak_coords, self.nonpeak_coords = peak_coords, nonpeak_coords
        
        # 保存参数
        self.negative_sampling_ratio = negative_sampling_ratio
        self.inputlen = inputlen
        self.outputlen = outputlen
        self.batch_size = batch_size
        self.add_revcomp = add_revcomp
        self.return_coords = return_coords
        self.shuffle_at_epoch_start = shuffle_at_epoch_start
        
        # ========================================
        # Step 2: 初始化时进行一次数据增强
        # ========================================
        self.crop_revcomp_data()
```

#### 负采样函数

```python
def subsample_nonpeak_data(nonpeak_seqs, nonpeak_cts, nonpeak_coords, 
                          peak_data_size, negative_sampling_ratio):
    """
    随机采样负样本
    
    为什么需要？
    - 通常nonpeaks数量 >> peaks数量
    - 使用所有nonpeaks会导致类别不平衡
    - 采样可以加快训练速度
    
    例子:
        peaks数量: 10000
        nonpeaks数量: 100000
        negative_sampling_ratio: 0.1
        
        采样后的nonpeaks: 10000 * 0.1 = 1000
        
        最终训练集: 10000 peaks + 1000 nonpeaks = 11000
    """
    # 计算要保留的负样本数量
    num_nonpeak_samples = int(negative_sampling_ratio * peak_data_size)
    
    # 随机选择索引
    nonpeak_indices_to_keep = np.random.choice(
        len(nonpeak_seqs), 
        size=num_nonpeak_samples, 
        replace=False  # 不重复采样
    )
    
    # 提取选中的样本
    nonpeak_seqs = nonpeak_seqs[nonpeak_indices_to_keep]
    nonpeak_cts = nonpeak_cts[nonpeak_indices_to_keep]
    nonpeak_coords = nonpeak_coords[nonpeak_indices_to_keep]
    
    return nonpeak_seqs, nonpeak_cts, nonpeak_coords
```

#### 数据裁剪和增强

```python
def crop_revcomp_data(self):
    """
    每个epoch调用一次
    
    执行:
    1. 随机裁剪peaks（jittering）
    2. 负采样nonpeaks
    3. 合并peaks和nonpeaks
    4. 反向互补增强
    5. 打乱顺序
    """
    
    # ========================================
    # 情况1: 同时有peaks和nonpeaks
    # ========================================
    if (self.peak_seqs is not None) and (self.nonpeak_seqs is not None):
        
        # Step 1: 随机裁剪峰数据
        # 从更宽的序列中随机裁剪出目标长度
        cropped_peaks, cropped_cnts, cropped_coords = \
            augment.random_crop(
                self.peak_seqs,      # (N, 2370, 4)
                self.peak_cts,       # (N, 1256)
                self.inputlen,       # 2114
                self.outputlen,      # 1000
                self.peak_coords
            )
        # 输出: (N, 2114, 4), (N, 1000), (N, 4)
        
        # Step 2: 负采样
        if self.negative_sampling_ratio < 1.0:
            sampled_nonpeak_seqs, sampled_nonpeak_cts, sampled_nonpeak_coords = \
                subsample_nonpeak_data(
                    self.nonpeak_seqs, 
                    self.nonpeak_cts, 
                    self.nonpeak_coords,
                    len(self.peak_seqs), 
                    self.negative_sampling_ratio
                )
            
            # Step 3: 合并数据
            self.seqs = np.vstack([cropped_peaks, sampled_nonpeak_seqs])
            self.cts = np.vstack([cropped_cnts, sampled_nonpeak_cts])
            self.coords = np.vstack([cropped_coords, sampled_nonpeak_coords])
        else:
            # 使用所有nonpeaks
            self.seqs = np.vstack([cropped_peaks, self.nonpeak_seqs])
            self.cts = np.vstack([cropped_cnts, self.nonpeak_cts])
            self.coords = np.vstack([cropped_coords, self.nonpeak_coords])
    
    # ========================================
    # 情况2: 只有peaks
    # ========================================
    elif self.peak_seqs is not None:
        cropped_peaks, cropped_cnts, cropped_coords = \
            augment.random_crop(
                self.peak_seqs, self.peak_cts, 
                self.inputlen, self.outputlen, self.peak_coords
            )
        
        self.seqs = cropped_peaks
        self.cts = cropped_cnts
        self.coords = cropped_coords
    
    # ========================================
    # 情况3: 只有nonpeaks
    # ========================================
    elif self.nonpeak_seqs is not None:
        self.seqs = self.nonpeak_seqs
        self.cts = self.nonpeak_cts
        self.coords = self.nonpeak_coords
    
    else:
        print("Both peak and non-peak arrays are empty")
    
    # ========================================
    # Step 4: 反向互补增强和打乱
    # ========================================
    self.cur_seqs, self.cur_cts, self.cur_coords = \
        augment.crop_revcomp_augment(
            self.seqs, self.cts, self.coords, 
            self.inputlen, self.outputlen, 
            self.add_revcomp,          # 是否添加反向互补
            shuffle=self.shuffle_at_epoch_start  # 是否打乱
        )
```

#### 获取批次

```python
def __getitem__(self, idx):
    """
    Keras会自动调用这个函数获取第idx个批次
    
    参数:
        idx: 批次索引 (0, 1, 2, ...)
    
    返回:
        (X, y) 或 (X, y, coords)
    """
    # 计算批次范围
    batch_start = idx * self.batch_size
    batch_end = (idx + 1) * self.batch_size
    
    # 提取批次数据
    batch_seq = self.cur_seqs[batch_start:batch_end]      # (B, 2114, 4)
    batch_cts = self.cur_cts[batch_start:batch_end]       # (B, 1000)
    batch_coords = self.cur_coords[batch_start:batch_end] # (B, 4)
    
    # 准备标签
    # y = [Profile, Counts]
    #   Profile: 原始计数 (B, 1000)
    #   Counts: log(1 + sum(计数)) (B, 1)
    
    if self.return_coords:
        return (
            batch_seq,  # X
            [
                batch_cts,  # y[0]: Profile标签
                np.log(1 + batch_cts.sum(-1, keepdims=True))  # y[1]: Counts标签
            ],
            batch_coords  # 坐标信息
        )
    else:
        return (
            batch_seq,  # X
            [
                batch_cts,  # y[0]: Profile标签
                np.log(1 + batch_cts.sum(-1, keepdims=True))  # y[1]: Counts标签
            ]
        )
```

**为什么Counts标签要用log(1 + sum)？**

```python
# 原始计数
profile = [0, 5, 10, 15, 10, 5, 0, ...]  # 长度1000
total_counts = sum(profile) = 45

# 模型预测的是log空间的counts
# 所以标签也要转换到log空间
log_counts = log(1 + 45) = log(46) = 3.83

# 为什么+1？
# 防止log(0)未定义
# 当总计数为0时，log(1+0) = log(1) = 0
```

#### Epoch结束时的操作

```python
def on_epoch_end(self):
    """
    Keras在每个epoch结束时自动调用
    
    作用:
    - 重新进行随机裁剪
    - 重新进行反向互补
    - 重新打乱数据
    
    结果:
    - 每个epoch训练的数据都略有不同
    - 增强泛化能力
    """
    self.crop_revcomp_data()
```

#### 其他方法

```python
def __len__(self):
    """
    返回每个epoch有多少个批次
    
    计算:
        总样本数 / 批次大小（向上取整）
    
    例子:
        样本数: 10050
        批次大小: 64
        批次数: ceil(10050/64) = 158
    """
    return math.ceil(self.seqs.shape[0] / self.batch_size)
```

---

## 数据增强详解

### 1. 随机裁剪 (Jittering)

**目的：** 让模型看到峰的不同"视角"

```python
def random_crop(seqs, cts, inputlen, outputlen, coords):
    """
    从更宽的序列中随机裁剪
    
    输入:
        seqs: (N, 2370, 4)    ← 宽序列
        cts: (N, 1256)        ← 宽信号
        inputlen: 2114        ← 目标序列长度
        outputlen: 1000       ← 目标信号长度
    
    输出:
        cropped_seqs: (N, 2114, 4)
        cropped_cts: (N, 1000)
    
    过程:
        对每个样本随机选择一个偏移量
        在[-128, +128]范围内
    """
```

**可视化：**

```
原始序列 (2370 bp):
[============================================]

可能的裁剪位置:
[===2114===]           ← 向左偏移128
    [===2114===]       ← 向左偏移64
        [===2114===]   ← 中心
            [===2114===]   ← 向右偏移64
                [===2114===]   ← 向右偏移128

每个epoch随机选择一个位置！
```

### 2. 反向互补 (Reverse Complement)

**目的：** DNA的两条链都包含信息

```python
def revcomp_augment(seqs, cts):
    """
    添加反向互补的序列
    
    DNA反向互补规则:
    A ↔ T
    C ↔ G
    方向相反
    
    例子:
    原序列: 5' ATCG 3'
           3' TAGC 5'
    
    反向互补: 5' CGAT 3'
             3' GCTA 5'
    """
    
    # 序列反向互补
    revcomp_seqs = seqs[:, ::-1, ::-1]  # 位置反向，碱基反向
    
    # 信号也要反向
    revcomp_cts = cts[:, ::-1]
    
    # 合并原始和反向互补
    aug_seqs = np.vstack([seqs, revcomp_seqs])
    aug_cts = np.vstack([cts, revcomp_cts])
    
    return aug_seqs, aug_cts
```

**为什么反向互补？**

```
DNA双链结构:
5' ATCGATCG 3'  ← 正链（我们测序的）
3' TAGCTAGC 5'  ← 反链

转录因子可以结合任一链！
添加反向互补 → 数据量翻倍 → 更好的泛化
```

### 3. 打乱 (Shuffling)

```python
def shuffle_data(seqs, cts, coords):
    """
    随机打乱样本顺序
    
    为什么？
    - 防止模型学习样本顺序
    - 确保每个批次都是混合的
    """
    indices = np.random.permutation(len(seqs))
    return seqs[indices], cts[indices], coords[indices]
```

---

## 完整数据流程示例

### 示例：训练ChromBPNet

```python
# ========================================
# 1. 准备输入文件
# ========================================
peaks_file = "peaks.bed"          # 10000个峰
nonpeaks_file = "nonpeaks.bed"    # 100000个背景区域
genome_file = "hg38.fa"
bigwig_file = "signal.bw"

# ========================================
# 2. 读取文件
# ========================================
peaks = pd.read_csv(peaks_file, sep='\t', 
                    names=['chr', 'start', 'end', 'name', 'score', 
                           'strand', 'signalValue', 'pValue', 'qValue', 'summit'])

nonpeaks = pd.read_csv(nonpeaks_file, sep='\t', names=...)

# ========================================
# 3. 创建生成器
# ========================================
train_generator = ChromBPNetBatchGenerator(
    peak_regions=peaks,
    nonpeak_regions=nonpeaks,
    genome_fasta=genome_file,
    batch_size=64,
    inputlen=2114,
    outputlen=1000,
    max_jitter=128,
    negative_sampling_ratio=0.1,  # 只用10%的nonpeaks
    cts_bw_file=bigwig_file,
    add_revcomp=True,              # 添加反向互补
    return_coords=False,
    shuffle_at_epoch_start=True    # 每个epoch打乱
)

# ========================================
# 4. 训练模型
# ========================================
model.fit(
    train_generator,
    epochs=50,
    verbose=1
)

# 生成器会自动:
# - 每个epoch重新裁剪数据
# - 添加反向互补
# - 打乱顺序
# - 按批次返回数据
```

### 数据流动追踪

```
Epoch 1:
  ├─ on_epoch_start
  │   ├─ 随机裁剪peaks: [2370,4] → [2114,4]
  │   ├─ 负采样nonpeaks: 100000 → 1000
  │   ├─ 合并: 10000 + 1000 = 11000
  │   ├─ 反向互补: 11000 → 22000
  │   └─ 打乱顺序
  │
  ├─ __getitem__(0): 返回批次0 (样本0-63)
  ├─ __getitem__(1): 返回批次1 (样本64-127)
  ├─ ...
  └─ __getitem__(343): 返回批次343 (最后一批)

Epoch 2:
  ├─ on_epoch_end (调用crop_revcomp_data)
  │   ├─ 重新随机裁剪 (不同的偏移！)
  │   ├─ 重新负采样 (不同的nonpeaks！)
  │   ├─ 重新反向互补
  │   └─ 重新打乱
  │
  └─ 重复训练...
```

---

## 总结

### 关键要点

1. **数据加载**
   - Peaks加载更宽的区域（允许jittering）
   - NonPeaks加载精确长度
   - 同时加载序列和信号

2. **批次生成**
   - 每个epoch重新增强数据
   - 支持负采样
   - 自动处理批次划分

3. **数据增强**
   - 随机裁剪（jittering）
   - 反向互补
   - 打乱顺序

4. **设计巧妙之处**
   - 使用Keras Sequence（支持多进程）
   - 每个epoch动态增强（无限数据变化）
   - 负采样解决类别不平衡

### 常见问题

**Q: 为什么需要jittering？**
A: 让模型学习峰的不同"视角"，提高泛化能力。

**Q: 为什么negative_sampling_ratio通常<1？**
A: 因为nonpeaks数量远大于peaks，全用会导致类别不平衡。

**Q: 反向互补会不会导致过拟合？**
A: 不会，因为DNA确实是双链的，这是真实的生物学特性。

**Q: 数据增强后数据量有多大？**
A: 
- 原始: N个peaks + M个nonpeaks
- 裁剪后: N + M*ratio
- 反向互补后: 2*(N + M*ratio)

---

## 下一步

现在您已经理解了数据处理流程，可以继续学习：
1. 模型训练流程
2. 模型评估方法
3. 如何解释模型预测

查看配套文件：
- `ChromBPNet完整项目深度解析.md` - 模型架构
- `ChromBPNet训练流程详解.md` - 训练过程
