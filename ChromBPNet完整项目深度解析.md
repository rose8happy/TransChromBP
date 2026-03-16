# 🔬 ChromBPNet 完整项目深度解析

> 从模型架构到数据处理的全方位讲解  
> 适合想要深入理解整个系统的学习者

---

## 📋 目录

1. [项目整体架构](#项目整体架构)
2. [两个模型的对比](#两个模型的对比)
3. [BPNet 模型详解](#bpnet-模型详解)
4. [ChromBPNet 带偏置模型详解](#chrombpnet-带偏置模型详解)
5. [数据处理流程](#数据处理流程)
6. [训练流程详解](#训练流程详解)
7. [预测流程详解](#预测流程详解)
8. [完整工作流程](#完整工作流程)

---

## 项目整体架构

### 核心思想：偏置因子分解 (Bias Factorization)

ChromBPNet 的核心创新是将 ATAC-seq/DNase-seq 信号**分解为两部分**：

```
总信号 = 酶偏好性信号 (Bias) + 真实生物学信号 (TF Binding)
```

**为什么需要这样做？**

ATAC-seq 和 DNase-seq 实验中使用的酶（Tn5 转座酶和 DNase I）本身对某些 DNA 序列有偏好性，这会**混淆**真实的转录因子结合信号。

```
观察到的信号
    ↓
┌────────────────────────┐
│  包含两种成分:         │
│  1. 酶的序列偏好      │  ← 噪声，需要去除
│  2. 真实的TF结合      │  ← 我们真正想要的
└────────────────────────┘
```

### 解决方案：两步训练

```
第一步: 训练偏置模型 (Bias Model)
  ↓
只使用背景区域（非peak）
学习酶的序列偏好性
  ↓
冻结这个模型
  
第二步: 训练ChromBPNet模型
  ↓
使用峰区域（peaks）
模型输出 = 无偏置预测 + 偏置预测
  ↓
学习真实的生物学信号
```

---

## 两个模型的对比

### 文件对比

| 特性 | bpnet_model.py | chrombpnet_with_bias_model.py |
|------|----------------|-------------------------------|
| **用途** | 基础模型 | 带偏置校正的完整模型 |
| **训练阶段** | 第一步：训练偏置模型 | 第二步：训练主模型 |
| **训练数据** | 背景区域（nonpeaks） | 峰区域（peaks） |
| **输入** | DNA序列 | DNA序列 |
| **输出** | Profile + Counts | Profile + Counts |
| **特殊处理** | 无 | 整合预训练的偏置模型 |
| **可训练** | 所有参数 | 主模型参数（偏置模型冻结） |

### 架构对比图

#### bpnet_model.py (简单模型)

```
输入: DNA序列
    ↓
┌─────────────────────┐
│   第一层卷积         │
│   膨胀卷积 × 9      │
│   残差连接          │
└─────────┬───────────┘
          │
    ┌─────┴─────┐
    ↓           ↓
┌─────────┐ ┌─────────┐
│Profile  │ │Counts   │
│预测     │ │预测     │
└─────────┘ └─────────┘
```

#### chrombpnet_with_bias_model.py (偏置因子分解模型)

```
输入: DNA序列
    │
    ├──────────────────────┐
    ↓                      ↓
┌──────────────┐    ┌──────────────┐
│ 偏置模型     │    │ 主BPNet模型  │
│ (冻结参数)   │    │ (可训练)     │
│              │    │              │
│ 学习:        │    │ 学习:        │
│ 酶序列偏好   │    │ TF结合信号   │
└──────┬───────┘    └──────┬───────┘
       │                   │
       ↓                   ↓
   Bias_Profile      TF_Profile
   Bias_Counts       TF_Counts
       │                   │
       └────────┬──────────┘
                ↓
         ┌──────────────┐
         │  组合输出    │
         ├──────────────┤
         │ Profile =    │
         │ TF + Bias    │
         │              │
         │ Counts =     │
         │ logsumexp(   │
         │  TF, Bias)   │
         └──────────────┘
```

---

## BPNet 模型详解

### 完整代码注释

```python
def getModelGivenModelOptionsAndWeightInits(args, model_params):
    """
    构建基础 BPNet 模型
    
    用途:
    1. 作为偏置模型（在背景区域训练）
    2. 作为独立模型（不考虑偏置校正）
    
    参数:
        args: 命令行参数
        model_params: 模型超参数字典
    
    返回:
        编译好的Keras模型
    """
    
    # ============================================================
    # 1. 参数设置
    # ============================================================
    
    # 固定参数
    conv1_kernel_size = 21      # 第一层卷积核大小
    profile_kernel_size = 75    # Profile头卷积核大小
    num_tasks = 1               # 单任务
    
    # 从配置读取的参数
    filters = int(model_params['filters'])                    # 通道数
    n_dil_layers = int(model_params['n_dil_layers'])          # 膨胀卷积层数
    counts_loss_weight = float(model_params['counts_loss_weight'])  # 损失权重
    sequence_len = int(model_params["inputlen"])              # 输入长度
    out_pred_len = int(model_params["outputlen"])             # 输出长度
    
    # ============================================================
    # 2. 设置随机种子（确保可重复）
    # ============================================================
    seed = args.seed
    np.random.seed(seed)    
    tf.random.set_seed(seed)
    rn.seed(seed)
    
    # ============================================================
    # 3. 构建模型
    # ============================================================
    
    # 输入层
    inp = Input(shape=(sequence_len, 4), name='sequence')
    
    # 第一层卷积（不使用膨胀）
    # 作用: 学习基本的DNA motifs
    x = Conv1D(filters,
               kernel_size=conv1_kernel_size,
               padding='valid', 
               activation='relu',
               name='bpnet_1st_conv')(inp)
    
    # 膨胀卷积残差模块
    for i in range(1, n_dil_layers + 1):
        # 膨胀卷积: dilation_rate = 2^i
        conv_x = Conv1D(filters, 
                        kernel_size=3, 
                        padding='valid',
                        activation='relu', 
                        dilation_rate=2**i,
                        name=f'bpnet_{i}conv')(x)
        
        # 对称裁剪以匹配维度
        x_len = int_shape(x)[1]
        conv_x_len = int_shape(conv_x)[1]
        x = Cropping1D((x_len - conv_x_len) // 2, 
                       name=f"bpnet_{i}crop")(x)
        
        # 残差连接
        x = add([conv_x, x])
    
    # 分支1: Profile预测
    prof_out_precrop = Conv1D(filters=num_tasks,
                              kernel_size=profile_kernel_size,
                              padding='valid',
                              name='prof_out_precrop')(x)
    
    cropsize = int(int_shape(prof_out_precrop)[1]/2) - int(out_pred_len/2)
    prof = Cropping1D(cropsize,
                      name='logits_profile_predictions_preflatten')(prof_out_precrop)
    
    profile_out = Flatten(name="logits_profile_predictions")(prof)
    
    # 分支2: Counts预测
    gap_combined_conv = GlobalAvgPool1D(name='gap')(x)
    count_out = Dense(num_tasks, name="logcount_predictions")(gap_combined_conv)
    
    # ============================================================
    # 4. 编译模型
    # ============================================================
    model = Model(inputs=[inp], outputs=[profile_out, count_out])
    
    model.compile(optimizer=Adam(learning_rate=args.learning_rate),
                  loss=[multinomial_nll, 'mse'],
                  loss_weights=[1, counts_loss_weight])
    
    return model
```

### 关键设计决策

#### 为什么这个模型可以既做偏置模型又做主模型？

**作为偏置模型时（第一步）：**
- 训练数据：背景区域（没有TF结合的区域）
- 学到的：酶的序列偏好性
- 例如：Tn5 喜欢在某些序列模式处切割

**作为独立模型时：**
- 训练数据：峰区域（有TF结合的区域）
- 学到的：TF结合信号 + 酶偏好性（混合）
- 问题：无法区分两者

---

## ChromBPNet 带偏置模型详解

### 核心创新：偏置因子分解

```python
def load_pretrained_bias(model_hdf5):
    """
    加载预训练的偏置模型并冻结参数
    
    为什么冻结？
    - 偏置模型已经在背景区域学好了酶偏好
    - 不希望在训练主模型时改变它
    - 确保偏置预测稳定
    """
    from tensorflow.keras.models import load_model
    from tensorflow.keras.utils import get_custom_objects
    
    # 注册自定义损失函数
    custom_objects = {"multinomial_nll": multinomial_nll, "tf": tf}    
    get_custom_objects().update(custom_objects)
    
    # 加载模型
    pretrained_bias_model = load_model(model_hdf5)
    
    # 🔒 冻结所有层（关键步骤！）
    num_layers = len(pretrained_bias_model.layers)
    for i in range(num_layers):
        pretrained_bias_model.layers[i].trainable = False
    
    return pretrained_bias_model
```

### 主模型构建

```python
def bpnet_model(filters, n_dil_layers, sequence_len, out_pred_len):
    """
    构建不带偏置的BPNet模型
    
    注意：层名称都加了 'wo_bias_bpnet_' 前缀
    这样可以区分主模型和偏置模型的层
    """
    
    conv1_kernel_size = 21
    profile_kernel_size = 75
    num_tasks = 1
    
    inp = Input(shape=(sequence_len, 4), name='sequence')
    
    # 第一层卷积
    x = Conv1D(filters,
               kernel_size=conv1_kernel_size,
               padding='valid', 
               activation='relu',
               name='wo_bias_bpnet_1st_conv')(inp)  # ← 注意前缀
    
    # 膨胀卷积层（与bpnet_model.py相同逻辑）
    for i in range(1, n_dil_layers + 1):
        conv_x = Conv1D(filters, 
                        kernel_size=3, 
                        padding='valid',
                        activation='relu', 
                        dilation_rate=2**i,
                        name=f'wo_bias_bpnet_{i}conv')(x)  # ← 注意前缀
        
        x_len = int_shape(x)[1]
        conv_x_len = int_shape(conv_x)[1]
        x = Cropping1D((x_len - conv_x_len) // 2, 
                       name=f"wo_bias_bpnet_{i}crop")(x)  # ← 注意前缀
        x = add([conv_x, x])
    
    # Profile分支
    prof_out_precrop = Conv1D(filters=num_tasks,
                              kernel_size=profile_kernel_size,
                              padding='valid',
                              name='wo_bias_bpnet_prof_out_precrop')(x)
    
    cropsize = int(int_shape(prof_out_precrop)[1]/2) - int(out_pred_len/2)
    prof = Cropping1D(cropsize,
                      name='wo_bias_bpnet_logitt_before_flatten')(prof_out_precrop)
    
    profile_out = Flatten(name="wo_bias_bpnet_logits_profile_predictions")(prof)
    
    # Counts分支
    gap_combined_conv = GlobalAvgPool1D(name='gap')(x)
    count_out = Dense(num_tasks, name="wo_bias_bpnet_logcount_predictions")(gap_combined_conv)
    
    # 🎯 关键：作为子模型
    model = Model(inputs=[inp], outputs=[profile_out, count_out], 
                  name="model_wo_bias")  # ← 命名很重要！
    
    return model
```

### 组合模型

```python
def getModelGivenModelOptionsAndWeightInits(args, model_params):
    """
    构建完整的ChromBPNet模型
    
    架构:
        输入
         ↓
    ┌────┴────┐
    ↓         ↓
  Bias    无偏置
  模型     模型
    ↓         ↓
  组合输出
    """
    
    # ============================================================
    # 1. 参数检查和设置
    # ============================================================
    assert "bias_model_path" in model_params.keys()  # 必须提供偏置模型
    
    filters = int(model_params['filters'])
    n_dil_layers = int(model_params['n_dil_layers'])
    counts_loss_weight = float(model_params['counts_loss_weight'])
    bias_model_path = model_params['bias_model_path']  # 偏置模型路径
    sequence_len = int(model_params['inputlen'])
    out_pred_len = int(model_params['outputlen'])
    
    # ============================================================
    # 2. 加载预训练模型
    # ============================================================
    
    # 加载并冻结偏置模型
    bias_model = load_pretrained_bias(bias_model_path)
    
    # 创建无偏置的主模型
    bpnet_model_wo_bias = bpnet_model(filters, n_dil_layers, 
                                      sequence_len, out_pred_len)
    
    # ============================================================
    # 3. 设置随机种子
    # ============================================================
    seed = args.seed
    np.random.seed(seed)    
    tf.random.set_seed(seed)
    rn.seed(seed)
    
    # ============================================================
    # 4. 组合两个模型
    # ============================================================
    
    # 创建新的输入层
    inp = Input(shape=(sequence_len, 4), name='sequence')
    
    # 获取两个模型的输出
    bias_output = bias_model(inp)              # [profile, counts]
    output_wo_bias = bpnet_model_wo_bias(inp)  # [profile, counts]
    
    # ============================================================
    # 5. 形状检查（防止出错）
    # ============================================================
    assert len(bias_output[1].shape) == 2       # counts: (None, 1)
    assert len(bias_output[0].shape) == 2       # profile: (None, out_pred_len)
    assert len(output_wo_bias[0].shape) == 2
    assert len(output_wo_bias[1].shape) == 2
    assert bias_output[1].shape[1] == 1
    assert bias_output[0].shape[1] == out_pred_len
    
    # ============================================================
    # 6. 组合输出 (关键创新!)
    # ============================================================
    
    # Profile: 简单相加
    # 原理: 信号可以线性叠加
    #       总Profile = TF引起的 + 酶偏好引起的
    profile_out = Add(name="logits_profile_predictions")(
        [output_wo_bias[0], bias_output[0]]
    )
    
    # Counts: logsumexp组合
    # 原理: counts在log空间，需要特殊组合
    #       总Counts = log(exp(TF_counts) + exp(bias_counts))
    concat_counts = Concatenate(axis=-1)(
        [output_wo_bias[1], bias_output[1]]
    )
    count_out = Lambda(
        lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True),
        name="logcount_predictions"
    )(concat_counts)
    
    # ============================================================
    # 7. 创建最终模型
    # ============================================================
    model = Model(inputs=[inp], outputs=[profile_out, count_out])
    
    model.compile(optimizer=Adam(learning_rate=args.learning_rate),
                  loss=[multinomial_nll, 'mse'],
                  loss_weights=[1, counts_loss_weight])
    
    return model
```

### 为什么这样组合？

#### Profile: 简单相加

```python
profile_out = Add()([output_wo_bias[0], bias_output[0]])
```

**数学原理：**
```
假设:
- TF在位置i产生信号 S_TF[i]
- 酶在位置i产生信号 S_bias[i]

观察到的总信号:
  S_total[i] = S_TF[i] + S_bias[i]

因此直接相加就可以！
```

**直观理解：**
```
位置:     1    2    3    4    5
TF信号:   0    5   10    5    0
酶信号:   2    2    2    2    2
总信号:   2    7   12    7    2  ← 简单相加
```

#### Counts: LogSumExp

```python
concat_counts = Concatenate(axis=-1)([output_wo_bias[1], bias_output[1]])
count_out = Lambda(lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True))(concat_counts)
```

**为什么不能简单相加？**

Counts 是在 **log 空间** 预测的：
```
模型输出: log(counts)
真实值:   counts
```

**数学原理：**
```
总计数 = TF产生的计数 + 酶产生的计数
       = exp(log_TF_counts) + exp(log_bias_counts)

取log:
log(总计数) = log(exp(log_TF) + exp(log_bias))
            = logsumexp([log_TF, log_bias])
```

**代码实现：**
```python
# 步骤1: 拼接两个log counts
concat = [log_TF_counts, log_bias_counts]  # shape: (batch, 2)

# 步骤2: logsumexp
result = log(exp(concat[:,0]) + exp(concat[:,1]))  # shape: (batch, 1)
```

**数值例子：**
```
log_TF_counts = 5.0      → TF计数 = exp(5.0) = 148
log_bias_counts = 3.0    → 酶计数 = exp(3.0) = 20

总计数 = 148 + 20 = 168
log(总计数) = log(168) = 5.12

logsumexp([5.0, 3.0]) = log(exp(5.0) + exp(3.0)) 
                       = log(148 + 20)
                       = log(168)
                       = 5.12  ✓
```

### 保存不带偏置的模型

```python
def save_model_without_bias(model, output_prefix):
    """
    保存只包含主模型（不含偏置）的版本
    
    用途:
    - 分析TF信号（不受酶偏好影响）
    - 迁移学习到其他数据集
    - 理解模型学到了什么
    """
    # 提取子模型的输出
    model_wo_bias = model.get_layer("model_wo_bias").output
    
    # 创建新模型
    model_without_bias = Model(
        inputs=model.get_layer("model_wo_bias").inputs,
        outputs=[model_wo_bias[0], model_wo_bias[1]]
    )
    
    print('save model without bias') 
    model_without_bias.save(output_prefix + "_nobias.h5")
```

---

## 数据处理流程

### 数据流程概览

```
原始数据
    ↓
┌──────────────────────────┐
│ 1. peaks.bed             │  峰区域（TF结合位点）
│ 2. nonpeaks.bed          │  背景区域（无TF结合）
│ 3. genome.fa             │  参考基因组
│ 4. signal.bw             │  信号文件(BigWig)
└────────────┬─────────────┘
             ↓
    ┌────────────────┐
    │ 数据加载器     │
    │ (data_utils)   │
    └────────┬───────┘
             ↓
    ┌────────────────┐
    │ Batch生成器    │
    │ (batchgen)     │
    └────────┬───────┘
             ↓
    ┌────────────────┐
    │ 数据增强       │
    │ - 随机裁剪     │
    │ - 反向互补     │
    └────────┬───────┘
             ↓
         训练模型
```

### 文件格式详解

#### 1. Peaks文件 (narrowPeak格式)

```
chr1  1000  1500  peak1  100  .  5.2  10.3  8.1  250
│     │     │     │      │    │  │    │     │    │
│     │     │     │      │    │  │    │     │    └─ 第10列: summit (峰值位置相对于start的偏移)
│     │     │     │      │    │  │    │     └────── 第9列: qValue
│     │     │     │      │    │  │    └──────────── 第8列: pValue
│     │     │     │      │    │  └───────────────── 第7列: signalValue
│     │     │     │      │    └──────────────────── 第6列: strand
│     │     │     │      └───────────────────────── 第5列: score
│     │     │     └──────────────────────────────── 第4列: name
│     │     └────────────────────────────────────── 第3列: end
│     └──────────────────────────────────────────── 第2列: start
└────────────────────────────────────────────────── 第1列: chromosome
```

**重要：**第10列的summit用于确定提取序列的中心：
```
中心位置 = start + summit
提取范围 = [中心 - inputlen/2, 中心 + inputlen/2]
```

#### 2. BigWig文件

存储每个碱基对的信号强度（ATAC-seq/DNase-seq切割位点）

```
位置:    1    2    3    4    5    6    7    8    ...
信号:   0.5  1.2  3.4  5.6  3.2  1.1  0.8  0.3  ...
```

### 数据加载器详解

<继续创建完整内容...>

