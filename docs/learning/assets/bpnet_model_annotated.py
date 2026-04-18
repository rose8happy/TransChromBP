"""
BPNet 模型 - 超详细注释版（小白友好） 🎓

这个文件是 bpnet_model.py 的完全注释版本
每一行都有详细的中文解释，帮助初学者理解
"""

import numpy as np
from tensorflow.keras.backend import int_shape
from tensorflow.keras.layers import Input, Cropping1D, add, Conv1D, GlobalAvgPool1D, Dense, Flatten
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model
from chrombpnet.training.utils.losses import multinomial_nll
import tensorflow as tf
import random as rn
import os

# 设置随机种子的哈希值，确保实验可重复
os.environ['PYTHONHASHSEED'] = '0'

def getModelGivenModelOptionsAndWeightInits(args, model_params):
    """
    构建 BPNet 模型的主函数

    参数:
        args: 命令行参数对象，包含学习率、随机种子等
        model_params: 字典，包含模型的超参数

    返回:
        model: 编译好的 Keras 模型
    """

    # ============================================================
    # 第一部分：设置模型参数
    # ============================================================

    # 默认参数（固定值，不能被覆盖）
    conv1_kernel_size = 21      # 第一层卷积的卷积核大小 (21 个碱基)
                                # 为什么是 21？这个大小适合捕获常见的 DNA motif

    profile_kernel_size = 75    # Profile 预测分支的卷积核大小 (75 个碱基)
                                # 为什么是 75？足够覆盖大多数转录因子结合位点

    num_tasks = 1               # 任务数量（1 = 单任务学习）
                                # 如果有多个细胞系，可以设置 >1 进行多任务学习

    # 从 model_params 字典中读取可配置的参数
    filters = int(model_params['filters'])  # 卷积核数量（通道数）
                                             # 典型值: 64, 128, 256
                                             # 越大模型容量越大，但也更容易过拟合

    n_dil_layers = int(model_params['n_dil_layers'])  # 膨胀卷积层数
                                                       # 典型值: 6-11
                                                       # 决定了感受野的大小

    counts_loss_weight = float(model_params['counts_loss_weight'])  # Counts 损失的权重
                                                                     # 典型值: 1-10
                                                                     # 用于平衡两个任务

    sequence_len = int(model_params["inputlen"])    # 输入序列长度（bp）
                                                     # 典型值: 2114

    out_pred_len = int(model_params["outputlen"])   # 输出预测长度（bp）
                                                     # 典型值: 1000
                                                     # 通常小于输入长度

    # 打印参数，方便调试和记录
    print("params:")
    print("filters:" + str(filters))
    print("n_dil_layers:" + str(n_dil_layers))
    print("conv1_kernel_size:" + str(conv1_kernel_size))
    print("profile_kernel_size:" + str(profile_kernel_size))
    print("counts_loss_weight:" + str(counts_loss_weight))

    # ============================================================
    # 第二部分：设置随机种子（确保实验可重复）
    # ============================================================

    seed = args.seed  # 从命令行参数获取随机种子
    np.random.seed(seed)        # NumPy 的随机种子
    tf.random.set_seed(seed)    # TensorFlow 的随机种子
    rn.seed(seed)               # Python 内置 random 模块的种子

    # 为什么需要设置随机种子？
    # - 神经网络的初始权重是随机的
    # - 训练过程中的数据打乱是随机的
    # - 设置种子后，每次运行结果都一样，便于调试和比较

    # ============================================================
    # 第三部分：定义输入层
    # ============================================================

    # 输入是 DNA 序列的 one-hot 编码
    # 形状: (batch_size, sequence_len, 4)
    # - sequence_len: 序列长度（例如 2114）
    # - 4: 代表 A, T, C, G 四种碱基
    #
    # 举例：序列 "ATCG" 会被编码为:
    # A = [1, 0, 0, 0]
    # T = [0, 1, 0, 0]
    # C = [0, 0, 1, 0]
    # G = [0, 0, 0, 1]
    inp = Input(shape=(sequence_len, 4), name='sequence')

    # ============================================================
    # 第四部分：第一层卷积（不使用膨胀）
    # ============================================================

    # 这一层的作用：学习识别基本的 DNA motifs（序列模式）
    #
    # 参数解释:
    # - filters: 学习多少个不同的模式（每个卷积核学习一个模式）
    # - kernel_size=21: 每次看 21 个连续的碱基
    # - padding='valid': 不填充，所以输出会变短
    # - activation='relu': ReLU 激活函数，引入非线性
    #
    # 输入形状: (batch_size, 2114, 4)
    # 输出形状: (batch_size, 2094, filters)
    # 长度变化: 2114 - 21 + 1 = 2094
    x = Conv1D(filters,
               kernel_size=conv1_kernel_size,
               padding='valid',
               activation='relu',
               name='bpnet_1st_conv')(inp)

    # ============================================================
    # 第五部分：膨胀卷积残差模块（核心！）⭐
    # ============================================================

    # 这个循环构建多层膨胀卷积，每层的膨胀率指数增长
    # 这是模型最重要的创新！

    # 准备层的名称（用于 Keras 中标识每一层）
    layer_names = [str(i) for i in range(1, n_dil_layers + 1)]

    # 循环构建 n_dil_layers 层膨胀卷积
    # 假设 n_dil_layers=9，那么会构建第 1,2,3,...,9 层
    for i in range(1, n_dil_layers + 1):

        # --------------------------------------------------
        # 步骤 1: 膨胀卷积
        # --------------------------------------------------

        # 膨胀率 = 2^i
        # i=1: dilation=2   → 感受野约 5 bp
        # i=2: dilation=4   → 感受野约 13 bp
        # i=3: dilation=8   → 感受野约 29 bp
        # ...
        # i=9: dilation=512 → 感受野约 2000 bp
        #
        # 什么是膨胀卷积？
        # 普通卷积: [*][*][*] 连续看 3 个位置
        # 膨胀=2:   [*]_[*]_[*] 跳着看，感受野更大
        # 膨胀=4:   [*]___[*]___[*] 跳更多，感受野更大

        conv_layer_name = 'bpnet_{}conv'.format(layer_names[i-1])

        conv_x = Conv1D(filters,
                        kernel_size=3,           # 卷积核大小固定为 3
                        padding='valid',         # 不填充
                        activation='relu',       # ReLU 激活
                        dilation_rate=2**i,      # 关键！膨胀率指数增长
                        name=conv_layer_name)(x)

        # 此时 conv_x 的长度会比 x 短（因为膨胀卷积）

        # --------------------------------------------------
        # 步骤 2: 对称裁剪（为了残差连接）
        # --------------------------------------------------

        # 问题: conv_x 和 x 的长度不同，无法直接相加
        # 解决: 把 x 裁剪到和 conv_x 一样的长度

        # 获取当前长度
        x_len = int_shape(x)[1]          # 例如: 1000
        conv_x_len = int_shape(conv_x)[1]  # 例如: 996

        # 确保长度差是偶数，这样才能对称裁剪
        assert((x_len - conv_x_len) % 2 == 0)

        # 计算需要裁剪多少
        # 例如: (1000 - 996) / 2 = 2
        # 意思是两边各裁剪 2 个位置
        crop_amount = (x_len - conv_x_len) // 2

        # Cropping1D(2) 会从两边各去掉 2 个位置
        # 裁剪前: [--][====================][--]
        # 裁剪后:     [====================]
        x = Cropping1D(crop_amount,
                       name="bpnet_{}crop".format(layer_names[i-1]))(x)

        # --------------------------------------------------
        # 步骤 3: 残差连接（跳跃连接）
        # --------------------------------------------------

        # 残差连接的核心思想：
        # 输出 = 原始输入 + 卷积后的特征
        #
        # 好处:
        # 1. 梯度可以直接流回去，避免梯度消失
        # 2. 保留了多尺度的特征信息
        # 3. 即使卷积层学不到东西，至少还有原始输入
        #
        # 图示:
        #     输入 x
        #      ├───────┐
        #      │  Conv │ → 学习新特征
        #      │   ↓   │
        #      │ conv_x│
        #      ├───────┘
        #      │
        #    相加 (残差连接)
        #      ↓
        #   新的 x

        x = add([conv_x, x])  # 逐元素相加

        # 现在 x 包含了:
        # - 之前所有层学到的特征
        # - 当前层学到的新特征
        # - 多尺度的感受野信息

    # 循环结束后，x 已经经过了多层膨胀卷积
    # 现在 x 包含了从局部到全局的多尺度特征

    # ============================================================
    # 第六部分：分支 1 - Profile 预测
    # ============================================================

    # Profile 预测的目标：预测信号在每个位置的强度分布
    # 输出是一个概率分布，表示每个位置的"开放程度"

    # --------------------------------------------------
    # 步骤 1.1: 大卷积核提取 profile 特征
    # --------------------------------------------------

    # 为什么用大卷积核 (75)？
    # - 转录因子结合位点通常长度为 6-20 bp
    # - 周围的序列上下文也很重要
    # - 75 bp 足够覆盖这些信息
    #
    # 这里把通道数降到 1（num_tasks=1）
    # 因为最终只需要预测一个任务的 profile
    prof_out_precrop = Conv1D(filters=num_tasks,
                              kernel_size=profile_kernel_size,
                              padding='valid',
                              name='prof_out_precrop')(x)

    # --------------------------------------------------
    # 步骤 1.2: 裁剪到目标输出长度
    # --------------------------------------------------

    # 目标: 把 prof_out_precrop 裁剪到 out_pred_len 的长度
    # 例如: 从 1200 裁剪到 1000

    # 计算裁剪大小
    # 当前长度的一半 - 目标长度的一半 = 需要裁剪的大小
    current_len = int_shape(prof_out_precrop)[1]
    cropsize = int(current_len / 2) - int(out_pred_len / 2)

    # 确保裁剪大小合理
    assert cropsize >= 0, "输出长度太大，无法裁剪！"
    assert (current_len % 2 == 0), "当前长度必须是偶数才能对称裁剪"

    # 对称裁剪
    # 例如: 从 1200 裁剪到 1000，两边各去掉 100
    prof = Cropping1D(cropsize,
                      name='logits_profile_predictions_preflatten')(prof_out_precrop)

    # 现在 prof 的形状是 (batch_size, out_pred_len, 1)

    # --------------------------------------------------
    # 步骤 1.3: 展平成一维向量
    # --------------------------------------------------

    # Flatten 把 (batch_size, out_pred_len, 1) 变成 (batch_size, out_pred_len)
    # 这个向量中的每个值代表对应位置的信号强度（logit 空间）
    #
    # 什么是 logit？
    # logit 是还没有归一化的原始预测值
    # 在损失函数中会转换成概率分布（通过 softmax）
    profile_out = Flatten(name="logits_profile_predictions")(prof)

    # profile_out 形状: (batch_size, out_pred_len)
    # 例如: (32, 1000) 表示批次大小 32，每个样本预测 1000 个位置

    # ============================================================
    # 第七部分：分支 2 - Counts 预测
    # ============================================================

    # Counts 预测的目标：预测总信号强度（总计数）
    # 这是一个标量值，表示整个区域的"开放程度"

    # --------------------------------------------------
    # 步骤 2.1: 全局平均池化
    # --------------------------------------------------

    # GlobalAvgPool1D 的作用：
    # 把空间维度"压缩"掉，只保留通道维度
    #
    # 输入 x 的形状: (batch_size, length, filters)
    # 输出 gap 的形状: (batch_size, filters)
    #
    # 计算方式：对每个通道，计算所有位置的平均值
    #
    # 举例:
    # x = [
    #   通道1: [0.1, 0.2, 0.3, ...] → 平均 = 0.2
    #   通道2: [0.5, 0.4, 0.6, ...] → 平均 = 0.5
    #   ...
    # ]
    # gap = [0.2, 0.5, ...]
    #
    # 为什么用全局平均池化？
    # - 它把整个序列的信息汇总成一个向量
    # - 对位置不敏感（无论 motif 在哪里，都会被捕获）
    # - 参数高效（不需要额外的参数）
    gap_combined_conv = GlobalAvgPool1D(name='gap')(x)

    # --------------------------------------------------
    # 步骤 2.2: 全连接层预测总计数
    # --------------------------------------------------

    # Dense 层学习如何组合全局特征来预测总计数
    # 输入: (batch_size, filters)
    # 输出: (batch_size, num_tasks)
    #
    # 计算方式: y = W·x + b
    # 其中 W 是权重矩阵，b 是偏置
    #
    # 输出是 log 空间的计数（log counts）
    # 为什么用 log？
    # - 计数可能跨越很大的范围（1 到 10000+）
    # - log 空间使得预测更稳定
    # - 在损失函数中会转回正常空间
    count_out = Dense(num_tasks, name="logcount_predictions")(gap_combined_conv)

    # count_out 形状: (batch_size, 1)
    # 每个样本一个标量值，表示预测的总计数

    # ============================================================
    # 第八部分：创建和编译模型
    # ============================================================

    # --------------------------------------------------
    # 创建 Keras Model 对象
    # --------------------------------------------------

    # 指定输入和输出
    # 输入: DNA 序列的 one-hot 编码
    # 输出: [profile 预测, counts 预测]
    model = Model(inputs=[inp], outputs=[profile_out, count_out])

    # --------------------------------------------------
    # 编译模型（指定优化器和损失函数）
    # --------------------------------------------------

    # 优化器: Adam
    # - 自适应学习率的优化算法
    # - 通常比 SGD 收敛更快
    # - learning_rate 从命令行参数获取（例如 0.001）

    # 损失函数:
    # [multinomial_nll, 'mse']
    #
    # 1. Profile 头: multinomial_nll (多项式负对数似然)
    #    - 适合预测概率分布
    #    - 把 logits 转换成概率，然后计算真实分布和预测分布的差异
    #    - 公式: -sum(y_true * log(softmax(y_pred)))
    #
    # 2. Counts 头: mse (均方误差)
    #    - 适合回归任务
    #    - 公式: mean((y_true - y_pred)^2)

    # 损失权重: [1, counts_loss_weight]
    # - Profile 损失的权重固定为 1
    # - Counts 损失的权重是可配置的（例如 5）
    # - 总损失 = 1 × profile_loss + 5 × counts_loss
    #
    # 为什么需要权重？
    # - 两个任务的损失值可能差异很大
    # - 权重用来平衡它们，避免某个任务主导训练

    model.compile(optimizer=Adam(learning_rate=args.learning_rate),
                  loss=[multinomial_nll, 'mse'],
                  loss_weights=[1, counts_loss_weight])

    # 返回编译好的模型，可以直接用于训练
    return model


def save_model_without_bias(model, output_prefix):
    """
    保存不带 bias 的模型

    注意: 这个函数在 bpnet_model.py 中实际上什么都不做
    因为这个架构本身就不是 bias-factorized 的

    这个函数存在的原因：
    - 所有模型架构文件都需要有这个函数
    - 这是为了代码的一致性和接口统一
    - 防止用户误用参数

    参数:
        model: Keras 模型
        output_prefix: 输出文件前缀
    """
    # 什么都不做，直接返回
    # 如果这是 chrombpnet_with_bias_model.py，这个函数会保存不带 bias 的模型
    return


# ============================================================
# 使用示例
# ============================================================

"""
如何使用这个模型？

1. 准备参数:
   model_params = {
       'filters': 64,
       'n_dil_layers': 9,
       'counts_loss_weight': 5.0,
       'inputlen': 2114,
       'outputlen': 1000
   }
   
   args = {
       'learning_rate': 0.001,
       'seed': 1234
   }

2. 创建模型:
   model = getModelGivenModelOptionsAndWeightInits(args, model_params)

3. 训练模型:
   model.fit(x_train, [y_profile_train, y_counts_train], 
             epochs=50, 
             batch_size=64)

4. 预测:
   profile_pred, counts_pred = model.predict(x_test)

5. 保存模型:
   model.save('my_bpnet_model.h5')
"""

# ============================================================
# 常见问题 FAQ
# ============================================================

"""
Q1: 为什么有两个输出头？
A1: 因为我们想同时学习两个相关的任务：
    - Profile: 信号的空间分布（形状）
    - Counts: 信号的总强度（大小）
    这两个任务互补，联合训练效果更好。

Q2: 为什么膨胀率是 2^i 而不是线性增长？
A2: 指数增长能更快地扩大感受野，用更少的层数看到更远的距离。
    例如: 9 层就能看到 2000+ bp，如果线性增长可能需要几十层。

Q3: 什么是残差连接？为什么需要它？
A3: 残差连接允许梯度直接传播回去，避免梯度消失。
    没有残差连接，深层网络很难训练。

Q4: 为什么 Profile 用 multinomial_nll，Counts 用 MSE？
A4: - Profile 是概率分布，multinomial_nll 是自然选择
    - Counts 是连续值，MSE 是标准的回归损失

Q5: 感受野到底是什么？
A5: 感受野是指输出中的一个位置能"看到"输入中多大的范围。
    感受野越大，模型能捕获的长距离依赖关系就越多。

Q6: 为什么输入长度 (2114) 大于输出长度 (1000)？
A6: 因为卷积会减小长度，而且我们只关心中心区域的预测。
    边缘部分的感受野不完整，所以会被裁剪掉。

Q7: 这个模型能用于其他任务吗？
A7: 可以！只要是序列到序列的任务，都可以借鉴这个架构：
    - 文本分类
    - 时间序列预测
    - 音频处理
    核心思想（膨胀卷积 + 残差连接）是通用的。

Q8: 训练需要多长时间？
A8: 取决于数据量和硬件：
    - 在 GPU 上，通常 1-2 天
    - 在 CPU 上，可能需要一周或更久
    - 建议使用 GPU 训练

Q9: 如何选择超参数？
A9: 常用的超参数范围：
    - filters: 64-256
    - n_dil_layers: 6-11
    - counts_loss_weight: 1-10
    - learning_rate: 0.0001-0.01
    建议从论文中的默认值开始，然后根据验证集表现调整。

Q10: 模型多大？需要多少内存？
A10: - 参数量: 约 100K-1M（取决于 filters 和 n_dil_layers）
     - 训练内存: 8-16 GB GPU 内存
     - 推理内存: 2-4 GB
"""
