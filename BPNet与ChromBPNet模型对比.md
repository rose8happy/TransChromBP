# 🔬 BPNet vs ChromBPNet 模型深度对比

> 详细对比两个模型的设计、用途和实现

---

## 快速对比表

| 特性 | bpnet_model.py | chrombpnet_with_bias_model.py |
|------|----------------|-------------------------------|
| **模型名称** | BPNet (基础模型) | ChromBPNet (偏置因子分解模型) |
| **训练阶段** | 第一步 | 第二步 |
| **训练数据** | 背景区域 (nonpeaks) | 峰区域 (peaks) |
| **学习目标** | 酶序列偏好性 | 真实TF结合信号 |
| **模型复杂度** | 简单 (单模型) | 复杂 (双模型组合) |
| **参数数量** | ~100K-1M | ~200K-2M (含冻结的偏置模型) |
| **可训练参数** | 全部 | 只有主模型部分 |
| **输出组合** | 直接输出 | 主模型 + 偏置模型 |
| **代码行数** | ~90行 | ~150行 |
| **使用场景** | 偏置模型训练 | 最终ChromBPNet训练 |

---

## 架构对比

### bpnet_model.py - 基础架构

```python
# 简洁的单模型架构

输入: DNA序列 (2114, 4)
    ↓
第一层卷积 (kernel=21)
    ↓
膨胀卷积残差模块 × 9层
│ Layer 1: dilation=2
│ Layer 2: dilation=4
│ ...
│ Layer 9: dilation=512
│ (每层都有残差连接)
    ↓
┌───────┴───────┐
↓               ↓
Conv1D(75)   GlobalAvgPool
↓               ↓
Crop         Dense(1)
↓               ↓
Flatten
↓               ↓
Profile_out   Counts_out
(1000,)       (1,)
```

**特点：**
- ✅ 结构简单清晰
- ✅ 容易理解和修改
- ✅ 训练速度快
- ❌ 无法分离偏置

### chrombpnet_with_bias_model.py - 组合架构

```python
# 复杂的双模型组合架构

输入: DNA序列 (2114, 4)
    │
    ├─────────────────────┐
    ↓                     ↓
┌─────────────┐    ┌──────────────┐
│ 偏置模型    │    │ 主BPNet模型  │
│ (Bias)      │    │ (wo_bias)    │
│             │    │              │
│ 🔒 冻结参数 │    │ ✏️  可训练   │
│             │    │              │
│ [预训练]    │    │ [从头训练]   │
└──────┬──────┘    └──────┬───────┘
       │                   │
       ↓                   ↓
  Bias_Profile       TF_Profile
  Bias_Counts        TF_Counts
       │                   │
       │                   │
       ↓                   ↓
  Profile组合:        Counts组合:
  Add([TF, Bias])    LogSumExp([TF, Bias])
       │                   │
       └─────────┬─────────┘
                 ↓
          最终输出
       (Profile, Counts)
```

**特点：**
- ✅ 能分离TF信号和酶偏好
- ✅ 生物学解释性强
- ✅ 预测更准确
- ❌ 结构复杂
- ❌ 需要两步训练

---

## 代码逐行对比

### 1. 模型函数签名

#### bpnet_model.py
```python
def getModelGivenModelOptionsAndWeightInits(args, model_params):
    """
    简单直接
    - 不需要额外参数
    - 直接构建完整模型
    """
```

#### chrombpnet_with_bias_model.py
```python
def getModelGivenModelOptionsAndWeightInits(args, model_params):
    """
    需要额外参数
    - 必须提供 model_params['bias_model_path']
    - 需要先加载预训练模型
    """
    # 关键检查
    assert "bias_model_path" in model_params.keys()
    bias_model_path = model_params['bias_model_path']
```

### 2. 参数设置

#### bpnet_model.py
```python
# 直接从配置读取
filters = int(model_params['filters'])
n_dil_layers = int(model_params['n_dil_layers'])
counts_loss_weight = float(model_params['counts_loss_weight'])
sequence_len = int(model_params["inputlen"])
out_pred_len = int(model_params["outputlen"])
```

#### chrombpnet_with_bias_model.py
```python
# 相同的参数
filters = int(model_params['filters'])
n_dil_layers = int(model_params['n_dil_layers'])
counts_loss_weight = float(model_params['counts_loss_weight'])
sequence_len = int(model_params["inputlen"])
out_pred_len = int(model_params["outputlen"])

# 额外的参数
bias_model_path = model_params['bias_model_path']  # ← 新增！
```

### 3. 模型构建

#### bpnet_model.py
```python
# 直接在主函数中构建模型
inp = Input(shape=(sequence_len, 4), name='sequence')

x = Conv1D(filters, kernel_size=21, ...)(inp)

for i in range(1, n_dil_layers + 1):
    conv_x = Conv1D(filters, kernel_size=3, dilation_rate=2**i, ...)(x)
    x = Cropping1D(...)(x)
    x = add([conv_x, x])

# Profile分支
prof_out_precrop = Conv1D(filters=1, kernel_size=75, ...)(x)
prof = Cropping1D(cropsize, ...)(prof_out_precrop)
profile_out = Flatten(name="logits_profile_predictions")(prof)

# Counts分支
gap = GlobalAvgPool1D(name='gap')(x)
count_out = Dense(1, name="logcount_predictions")(gap)

# 创建模型
model = Model(inputs=[inp], outputs=[profile_out, count_out])
```

#### chrombpnet_with_bias_model.py
```python
# 分两步构建

# Step 1: 定义辅助函数构建主模型
def bpnet_model(filters, n_dil_layers, sequence_len, out_pred_len):
    """构建不带偏置的BPNet模型（与bpnet_model.py逻辑相同）"""
    # 注意：所有层名都加了 'wo_bias_bpnet_' 前缀
    inp = Input(shape=(sequence_len, 4), name='sequence')
    x = Conv1D(filters, kernel_size=21, 
               name='wo_bias_bpnet_1st_conv')(inp)  # ← 前缀
    
    for i in range(1, n_dil_layers + 1):
        conv_x = Conv1D(filters, kernel_size=3, dilation_rate=2**i,
                       name=f'wo_bias_bpnet_{i}conv')(x)  # ← 前缀
        x = Cropping1D(..., name=f"wo_bias_bpnet_{i}crop")(x)
        x = add([conv_x, x])
    
    # ... 其余类似
    
    # 关键：命名为子模型
    return Model(inputs=[inp], outputs=[profile_out, count_out],
                name="model_wo_bias")  # ← 重要的名字！

# Step 2: 加载并冻结偏置模型
def load_pretrained_bias(model_hdf5):
    pretrained_bias_model = load_model(model_hdf5)
    
    # 🔒 冻结所有层
    for i in range(len(pretrained_bias_model.layers)):
        pretrained_bias_model.layers[i].trainable = False
    
    return pretrained_bias_model

# Step 3: 在主函数中组合
def getModelGivenModelOptionsAndWeightInits(args, model_params):
    # 加载两个模型
    bias_model = load_pretrained_bias(bias_model_path)
    bpnet_model_wo_bias = bpnet_model(filters, n_dil_layers, ...)
    
    # 创建新输入
    inp = Input(shape=(sequence_len, 4), name='sequence')
    
    # 获取两个模型的输出
    bias_output = bias_model(inp)
    output_wo_bias = bpnet_model_wo_bias(inp)
    
    # 组合输出
    profile_out = Add(name="logits_profile_predictions")(
        [output_wo_bias[0], bias_output[0]]
    )
    
    concat_counts = Concatenate(axis=-1)(
        [output_wo_bias[1], bias_output[1]]
    )
    count_out = Lambda(
        lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True),
        name="logcount_predictions"
    )(concat_counts)
    
    # 创建最终模型
    model = Model(inputs=[inp], outputs=[profile_out, count_out])
```

### 4. 层命名策略

#### bpnet_model.py
```python
# 简单命名
'bpnet_1st_conv'
'bpnet_1conv'
'bpnet_2conv'
...
'logits_profile_predictions'
'logcount_predictions'
```

#### chrombpnet_with_bias_model.py
```python
# 带前缀的命名（区分不同模型的层）

# 主模型（无偏置）
'wo_bias_bpnet_1st_conv'     # wo = without
'wo_bias_bpnet_1conv'
'wo_bias_bpnet_2conv'
...
'wo_bias_bpnet_logits_profile_predictions'
'wo_bias_bpnet_logcount_predictions'

# 偏置模型（保留原名）
'bpnet_1st_conv'  # 从预训练模型加载
'bpnet_1conv'
...

# 最终组合层
'logits_profile_predictions'  # 组合后的输出
'logcount_predictions'        # 组合后的输出
```

**为什么需要不同的命名？**
```python
# 如果不区分命名，会导致混淆：
bias_model有: 'bpnet_1conv'
main_model也有: 'bpnet_1conv'
→ 冲突！

# 使用前缀后：
bias_model有: 'bpnet_1conv'
main_model有: 'wo_bias_bpnet_1conv'
→ 清晰区分！
```

### 5. 输出组合方式

#### bpnet_model.py
```python
# 直接输出，不需要组合
outputs = [profile_out, count_out]
```

#### chrombpnet_with_bias_model.py
```python
# 复杂的组合逻辑

# Profile: 简单相加
profile_out = Add(name="logits_profile_predictions")(
    [output_wo_bias[0], bias_output[0]]
)

# Counts: LogSumExp组合
concat_counts = Concatenate(axis=-1)(
    [output_wo_bias[1], bias_output[1]]
)
count_out = Lambda(
    lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True),
    name="logcount_predictions"
)(concat_counts)
```

**数学原理对比：**

| 输出类型 | bpnet_model | chrombpnet | 原因 |
|---------|-------------|------------|------|
| Profile | 直接输出 | Add([TF, Bias]) | 信号可线性叠加 |
| Counts | 直接输出 | LogSumExp([TF, Bias]) | Log空间需要特殊组合 |

### 6. 模型保存

#### bpnet_model.py
```python
def save_model_without_bias(model, output_prefix):
    """
    什么都不做
    因为这个模型本身就不是偏置因子分解的
    """
    return
```

#### chrombpnet_with_bias_model.py
```python
def save_model_without_bias(model, output_prefix):
    """
    提取并保存主模型（不含偏置）
    
    用途：
    - 分析TF信号（排除酶偏好）
    - 迁移学习
    - 模型解释
    """
    # 提取子模型
    model_wo_bias = model.get_layer("model_wo_bias").output
    
    # 创建新模型
    model_without_bias = Model(
        inputs=model.get_layer("model_wo_bias").inputs,
        outputs=[model_wo_bias[0], model_wo_bias[1]]
    )
    
    # 保存
    model_without_bias.save(output_prefix + "_nobias.h5")
```

---

## 训练流程对比

### bpnet_model.py 的训练流程

```python
# ========================================
# 第一步：训练偏置模型
# ========================================

# 1. 准备数据
train_data = load_nonpeak_regions()  # 只用背景区域！

# 2. 创建模型
bias_model = bpnet_model.getModelGivenModelOptionsAndWeightInits(
    args, model_params
)

# 3. 训练
bias_model.fit(
    train_data,
    epochs=50,
    verbose=1
)

# 4. 保存
bias_model.save("bias_model.h5")

# 学到什么？
# → 酶在不同DNA序列上的切割偏好
# → 与TF结合无关的背景信号
```

### chrombpnet_with_bias_model.py 的训练流程

```python
# ========================================
# 第二步：训练ChromBPNet
# ========================================

# 1. 加载偏置模型
model_params['bias_model_path'] = "bias_model.h5"

# 2. 准备数据
train_data = load_peak_regions()  # 用峰区域！

# 3. 创建模型
chrombpnet_model = chrombpnet_with_bias_model.getModelGivenModelOptionsAndWeightInits(
    args, model_params
)

# 模型内部：
# - 偏置模型：冻结参数，提供偏置预测
# - 主模型：可训练，学习TF信号

# 4. 训练
chrombpnet_model.fit(
    train_data,
    epochs=50,
    verbose=1
)

# 5. 保存完整模型
chrombpnet_model.save("chrombpnet_full.h5")

# 6. 保存不带偏置的模型
save_model_without_bias(chrombpnet_model, "chrombpnet")
# → 生成 "chrombpnet_nobias.h5"

# 学到什么？
# → TF结合导致的信号
# → 排除了酶偏好的干扰
```

---

## 使用场景对比

### bpnet_model.py 的使用场景

#### ✅ 适用场景

1. **训练偏置模型**
   ```python
   # 在背景区域训练
   # 学习酶的序列偏好
   bias_model = train_bpnet_on_nonpeaks()
   ```

2. **作为独立模型**
   ```python
   # 如果不关心偏置校正
   # 可以直接在peaks上训练
   simple_model = train_bpnet_on_peaks()
   ```

3. **快速原型**
   ```python
   # 快速测试想法
   # 代码简单，易于修改
   ```

#### ❌ 不适用场景

1. **需要偏置校正**
   - 无法分离TF信号和酶偏好
   
2. **高精度要求**
   - 性能不如ChromBPNet
   
3. **可解释性要求**
   - 无法单独分析TF效应

### chrombpnet_with_bias_model.py 的使用场景

#### ✅ 适用场景

1. **最终模型训练**
   ```python
   # 生产环境
   # 需要最佳性能
   chrombpnet = train_chrombpnet_with_bias()
   ```

2. **科研分析**
   ```python
   # 需要理解TF作用
   # 需要排除酶效应
   tf_effects = chrombpnet_nobias.predict(seq)
   ```

3. **突变效应预测**
   ```python
   # 预测SNP的影响
   # 只关心TF信号的变化
   ```

#### ❌ 不适用场景

1. **快速原型**
   - 需要先训练偏置模型（额外步骤）
   
2. **没有预训练偏置模型**
   - 必须先训练bpnet_model
   
3. **计算资源受限**
   - 需要更多内存（两个模型）

---

## 性能对比

### 模型大小

```python
# bpnet_model (filters=64, n_dil_layers=9)
参数总数: ~1,000,000
模型文件: ~4 MB

# chrombpnet_with_bias_model
参数总数: ~2,000,000 (主模型 + 偏置模型)
可训练参数: ~1,000,000 (只有主模型)
冻结参数: ~1,000,000 (偏置模型)
模型文件: ~8 MB (完整) + ~4 MB (无偏置)
```

### 训练速度

```python
# bpnet_model
训练速度: 快
- 单模型
- 简单前向传播

# chrombpnet_with_bias_model  
训练速度: 较慢 (~1.5x)
- 双模型前向传播
- 额外的组合操作
```

### 预测精度

```python
# 实验数据（参考论文）

# bpnet_model (无偏置校正)
Profile MNLL: 0.25
Counts MSE: 0.18

# chrombpnet_with_bias_model (带偏置校正)
Profile MNLL: 0.20  ← 更低（更好）
Counts MSE: 0.15    ← 更低（更好）

# 改进: ~20-25%
```

---

## 实际使用建议

### 何时使用 bpnet_model.py？

```python
# 场景1: 训练偏置模型（必需步骤）
python train.py \
    --architecture_from_file bpnet_model.py \
    --peaks None \
    --nonpeaks background.bed \
    --output_prefix bias_model

# 场景2: 快速测试
# 不关心偏置校正，只想快速得到结果
python train.py \
    --architecture_from_file bpnet_model.py \
    --peaks peaks.bed \
    --nonpeaks None \
    --output_prefix quick_test

# 场景3: 教学演示
# 代码简单，容易理解
```

### 何时使用 chrombpnet_with_bias_model.py？

```python
# 场景1: 最终模型训练（推荐）
python train.py \
    --architecture_from_file chrombpnet_with_bias_model.py \
    --peaks peaks.bed \
    --nonpeaks background.bed \
    --params params.tsv \  # 包含 bias_model_path
    --output_prefix chrombpnet

# params.tsv 内容:
# filters	64
# n_dil_layers	9
# counts_loss_weight	5
# inputlen	2114
# outputlen	1000
# bias_model_path	bias_model.h5  ← 关键！

# 场景2: 科研发表
# 需要最佳性能和可解释性

# 场景3: 生产部署
# 已有预训练的偏置模型
```

---

## 常见问题

### Q1: 能否跳过bpnet_model，直接用chrombpnet？
**A:** 不能。必须先用bpnet_model训练偏置模型，然后才能训练chrombpnet。

### Q2: 偏置模型需要每次重新训练吗？
**A:** 不需要。同一个实验protocol的偏置模型可以重复使用。

### Q3: 两个模型的参数必须一样吗？
**A:** 通常一样，但可以不同。例如：
```python
# 偏置模型: 简单一些
filters=32, n_dil_layers=6

# 主模型: 复杂一些
filters=64, n_dil_layers=9
```

### Q4: 能否用其他数据集的偏置模型？
**A:** 可以，如果：
- 相同的实验protocol (ATAC vs DNase)
- 相同的参考基因组
- 类似的数据质量

### Q5: chrombpnet训练时会改变偏置模型吗？
**A:** 不会。偏置模型的参数被冻结（`trainable=False`）。

---

## 代码差异总结

### 相同之处

1. **核心架构**
   - 都使用膨胀卷积 + 残差连接
   - 都有Profile和Counts两个输出头
   - 卷积核大小和激活函数相同

2. **训练设置**
   - 相同的损失函数（multinomial_nll + MSE）
   - 相同的优化器（Adam）
   - 相同的损失权重平衡策略

3. **数据格式**
   - 相同的输入格式（序列 one-hot）
   - 相同的输出格式（Profile + Counts）

### 不同之处

| 方面 | bpnet_model | chrombpnet |
|------|-------------|------------|
| **结构** | 单模型 | 双模型组合 |
| **层命名** | 简单命名 | 带前缀命名 |
| **输出** | 直接输出 | 组合输出（Add/LogSumExp） |
| **依赖** | 无 | 需要预训练的偏置模型 |
| **可训练参数** | 全部 | 只有主模型 |
| **保存功能** | 简单 | 可保存无偏置版本 |

---

## 总结

### BPNet (bpnet_model.py)
- 📌 **定位**: 基础模型，偏置模型训练
- ✅ **优势**: 简单、快速、易理解
- ❌ **劣势**: 无法分离偏置
- 🎯 **用途**: 第一步训练，快速原型

### ChromBPNet (chrombpnet_with_bias_model.py)  
- 📌 **定位**: 完整模型，偏置因子分解
- ✅ **优势**: 高精度、可解释、分离偏置
- ❌ **劣势**: 复杂、需两步训练
- 🎯 **用途**: 最终模型，科研分析

### 推荐工作流程

```
Step 1: 用 bpnet_model.py 训练偏置模型
    ↓
Step 2: 用 chrombpnet_with_bias_model.py 训练主模型
    ↓
Step 3: 使用ChromBPNet进行预测和分析
```

这就是两个模型的完整对比！希望能帮助您理解它们的异同和使用场景。
