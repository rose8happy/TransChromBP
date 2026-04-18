"""
BPNet 模型可视化工具 🎨

这个脚本帮助您可视化和理解 BPNet 模型的各个部分
运行这个脚本可以生成模型架构图和中间层输出
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import seaborn as sns

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def visualize_dna_encoding():
    """可视化 DNA 序列的 one-hot 编码"""

    # DNA 序列
    sequence = "ATCGATCG"
    bases = ['A', 'T', 'C', 'G']

    # One-hot 编码
    encoding = {
        'A': [1, 0, 0, 0],
        'T': [0, 1, 0, 0],
        'C': [0, 0, 1, 0],
        'G': [0, 0, 0, 1]
    }

    # 创建编码矩阵
    matrix = np.array([encoding[base] for base in sequence])

    # 可视化
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))

    # 上图：序列
    ax1.text(0.5, 0.5, sequence, ha='center', va='center',
             fontsize=24, family='monospace', weight='bold')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.axis('off')
    ax1.set_title('DNA 序列', fontsize=16, weight='bold', pad=20)

    # 下图：编码矩阵
    im = ax2.imshow(matrix.T, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
    ax2.set_yticks(range(4))
    ax2.set_yticklabels(bases, fontsize=12)
    ax2.set_xticks(range(len(sequence)))
    ax2.set_xticklabels(list(sequence), fontsize=12)
    ax2.set_xlabel('位置', fontsize=12, weight='bold')
    ax2.set_ylabel('碱基', fontsize=12, weight='bold')
    ax2.set_title('One-hot 编码 (形状: 8 × 4)', fontsize=16, weight='bold', pad=20)

    # 添加网格
    ax2.set_xticks(np.arange(len(sequence)) - 0.5, minor=True)
    ax2.set_yticks(np.arange(4) - 0.5, minor=True)
    ax2.grid(which='minor', color='gray', linewidth=0.5)

    # 添加数值
    for i in range(4):
        for j in range(len(sequence)):
            text = ax2.text(j, i, int(matrix[j, i]),
                          ha="center", va="center", color="black", fontsize=10)

    plt.colorbar(im, ax=ax2, label='值')
    plt.tight_layout()
    plt.savefig('1_dna_encoding.png', dpi=150, bbox_inches='tight')
    print("✅ 保存: 1_dna_encoding.png")
    plt.close()


def visualize_receptive_field():
    """可视化膨胀卷积的感受野增长"""

    n_layers = 9
    dilations = [2**i for i in range(1, n_layers + 1)]

    # 简化的感受野计算
    receptive_fields = []
    rf = 1
    for d in dilations:
        rf = rf + 2 * d  # 简化计算
        receptive_fields.append(rf)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # 左图：膨胀率增长
    layers = list(range(1, n_layers + 1))
    ax1.plot(layers, dilations, 'o-', linewidth=2, markersize=10, color='#2E86AB')
    ax1.fill_between(layers, dilations, alpha=0.3, color='#2E86AB')
    ax1.set_xlabel('层数', fontsize=14, weight='bold')
    ax1.set_ylabel('膨胀率 (Dilation Rate)', fontsize=14, weight='bold')
    ax1.set_title('膨胀率的指数增长', fontsize=16, weight='bold', pad=20)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_yscale('log')

    # 标注关键点
    for i, (layer, dilation) in enumerate(zip(layers, dilations)):
        if i % 2 == 0:
            ax1.annotate(f'{dilation}', xy=(layer, dilation),
                        xytext=(0, 10), textcoords='offset points',
                        ha='center', fontsize=10, weight='bold')

    # 右图：感受野增长
    ax2.plot(layers, receptive_fields, 's-', linewidth=2, markersize=10, color='#A23B72')
    ax2.fill_between(layers, receptive_fields, alpha=0.3, color='#A23B72')
    ax2.set_xlabel('层数', fontsize=14, weight='bold')
    ax2.set_ylabel('感受野大小 (bp)', fontsize=14, weight='bold')
    ax2.set_title('感受野的累积增长', fontsize=16, weight='bold', pad=20)
    ax2.grid(True, alpha=0.3, linestyle='--')

    # 标注关键层
    key_layers = [1, 5, 9]
    for layer in key_layers:
        idx = layer - 1
        ax2.annotate(f'{receptive_fields[idx]} bp',
                    xy=(layer, receptive_fields[idx]),
                    xytext=(10, 10), textcoords='offset points',
                    ha='left', fontsize=10, weight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    plt.tight_layout()
    plt.savefig('2_receptive_field_growth.png', dpi=150, bbox_inches='tight')
    print("✅ 保存: 2_receptive_field_growth.png")
    plt.close()


def visualize_convolution_types():
    """可视化不同类型的卷积"""

    fig, axes = plt.subplots(3, 1, figsize=(14, 10))

    sequence_len = 15

    # 准备数据
    convolution_types = [
        ("普通卷积 (dilation=1)", 1, '#2E86AB'),
        ("膨胀卷积 (dilation=2)", 2, '#A23B72'),
        ("膨胀卷积 (dilation=4)", 4, '#F18F01')
    ]

    for ax, (title, dilation, color) in zip(axes, convolution_types):
        # 画输入序列
        for i in range(sequence_len):
            rect = patches.Rectangle((i, 0), 0.9, 0.9,
                                     linewidth=1, edgecolor='black',
                                     facecolor='lightblue', alpha=0.5)
            ax.add_patch(rect)
            ax.text(i + 0.45, 0.45, str(i+1), ha='center', va='center', fontsize=8)

        # 画卷积核覆盖的位置 (kernel_size=3)
        kernel_positions = [0, dilation, 2*dilation]
        for pos in kernel_positions:
            if pos < sequence_len:
                rect = patches.Rectangle((pos, 0), 0.9, 0.9,
                                         linewidth=3, edgecolor=color,
                                         facecolor=color, alpha=0.7)
                ax.add_patch(rect)

        # 标注感受野
        if max(kernel_positions) < sequence_len:
            ax.annotate('', xy=(max(kernel_positions) + 0.9, 1.5),
                       xytext=(0, 1.5),
                       arrowprops=dict(arrowstyle='<->', lw=2, color=color))
            rf = max(kernel_positions) + 1
            ax.text((max(kernel_positions) + 0.9) / 2, 1.8,
                   f'感受野: {rf} 位置', ha='center', fontsize=11,
                   weight='bold', color=color)

        ax.set_xlim(-0.5, sequence_len)
        ax.set_ylim(-0.5, 2.5)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(title, fontsize=14, weight='bold', pad=20)

    plt.tight_layout()
    plt.savefig('3_convolution_types.png', dpi=150, bbox_inches='tight')
    print("✅ 保存: 3_convolution_types.png")
    plt.close()


def visualize_residual_connection():
    """可视化残差连接"""

    fig, ax = plt.subplots(figsize=(10, 8))

    # 定义位置
    positions = {
        'input': (2, 7),
        'conv': (5, 7),
        'crop': (2, 4),
        'add': (5, 4),
        'output': (5, 1)
    }

    colors = {
        'input': '#2E86AB',
        'conv': '#A23B72',
        'crop': '#F18F01',
        'add': '#C73E1D',
        'output': '#6A994E'
    }

    # 画框
    boxes = {
        'input': FancyBboxPatch((1, 6.5), 2, 1, boxstyle="round,pad=0.1",
                                edgecolor='black', facecolor=colors['input'],
                                linewidth=2, alpha=0.7),
        'conv': FancyBboxPatch((4, 6.5), 2, 1, boxstyle="round,pad=0.1",
                               edgecolor='black', facecolor=colors['conv'],
                               linewidth=2, alpha=0.7),
        'crop': FancyBboxPatch((1, 3.5), 2, 1, boxstyle="round,pad=0.1",
                               edgecolor='black', facecolor=colors['crop'],
                               linewidth=2, alpha=0.7),
        'add': FancyBboxPatch((4, 3.5), 2, 1, boxstyle="round,pad=0.1",
                              edgecolor='black', facecolor=colors['add'],
                              linewidth=2, alpha=0.7),
        'output': FancyBboxPatch((4, 0.5), 2, 1, boxstyle="round,pad=0.1",
                                 edgecolor='black', facecolor=colors['output'],
                                 linewidth=2, alpha=0.7)
    }

    for box in boxes.values():
        ax.add_patch(box)

    # 添加文字
    texts = {
        'input': '输入 x\n(长度: 1000)',
        'conv': 'Conv1D\n(膨胀卷积)',
        'crop': 'Cropping\n(对称裁剪)',
        'add': 'Add\n(残差连接)',
        'output': '输出 x\'\n(保留多尺度特征)'
    }

    for key, (x, y) in positions.items():
        ax.text(x, y, texts[key], ha='center', va='center',
               fontsize=11, weight='bold', color='white')

    # 画箭头
    arrows = [
        (positions['input'], positions['conv'], 'solid'),  # input -> conv
        (positions['input'], positions['crop'], 'solid'),  # input -> crop (跳跃)
        (positions['conv'], positions['add'], 'solid'),    # conv -> add
        (positions['crop'], positions['add'], 'solid'),    # crop -> add
        (positions['add'], positions['output'], 'solid')   # add -> output
    ]

    for start, end, style in arrows:
        arrow = FancyArrowPatch(start, end, arrowstyle='->',
                               mutation_scale=20, linewidth=2.5,
                               color='black', linestyle=style)
        ax.add_patch(arrow)

    # 标注跳跃连接
    ax.annotate('跳跃连接\n(残差)', xy=(1.5, 5.5), xytext=(0.2, 5.5),
               fontsize=10, weight='bold', color='red',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8),
               arrowprops=dict(arrowstyle='->', lw=2, color='red'))

    ax.set_xlim(0, 7)
    ax.set_ylim(0, 8.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('残差连接 (Residual Connection)',
                fontsize=16, weight='bold', pad=20)

    plt.tight_layout()
    plt.savefig('4_residual_connection.png', dpi=150, bbox_inches='tight')
    print("✅ 保存: 4_residual_connection.png")
    plt.close()


def visualize_profile_prediction():
    """可视化 Profile 预测示例"""

    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    # 生成模拟数据
    x = np.linspace(0, 1000, 1000)

    # 真实 profile (两个峰)
    true_profile = (30 * np.exp(-((x - 300)**2) / 5000) +
                   50 * np.exp(-((x - 700)**2) / 3000) +
                   5 * np.random.randn(1000))
    true_profile = np.maximum(true_profile, 0)

    # 预测 profile (略有偏差)
    pred_profile = (28 * np.exp(-((x - 305)**2) / 5200) +
                   48 * np.exp(-((x - 695)**2) / 3100) +
                   3 * np.random.randn(1000))
    pred_profile = np.maximum(pred_profile, 0)

    # 上图：真实 vs 预测
    axes[0].plot(x, true_profile, label='真实 Profile',
                linewidth=2, color='#2E86AB', alpha=0.7)
    axes[0].fill_between(x, true_profile, alpha=0.3, color='#2E86AB')
    axes[0].plot(x, pred_profile, label='预测 Profile',
                linewidth=2, color='#C73E1D', linestyle='--', alpha=0.7)
    axes[0].fill_between(x, pred_profile, alpha=0.2, color='#C73E1D')

    # 标注峰值
    peak1_idx = np.argmax(true_profile[200:400]) + 200
    peak2_idx = np.argmax(true_profile[600:800]) + 600

    axes[0].scatter([x[peak1_idx], x[peak2_idx]],
                   [true_profile[peak1_idx], true_profile[peak2_idx]],
                   s=200, color='red', marker='*', zorder=5,
                   label='转录因子结合位点')

    axes[0].set_xlabel('基因组位置 (bp)', fontsize=12, weight='bold')
    axes[0].set_ylabel('信号强度 (reads)', fontsize=12, weight='bold')
    axes[0].set_title('Profile 预测示例', fontsize=14, weight='bold', pad=15)
    axes[0].legend(fontsize=11, loc='upper right')
    axes[0].grid(True, alpha=0.3, linestyle='--')

    # 下图：误差分析
    error = pred_profile - true_profile
    axes[1].plot(x, error, linewidth=1.5, color='purple', alpha=0.7)
    axes[1].fill_between(x, error, alpha=0.3, color='purple')
    axes[1].axhline(y=0, color='black', linestyle='--', linewidth=1)

    # 标注统计信息
    mae = np.mean(np.abs(error))
    mse = np.mean(error**2)
    axes[1].text(50, np.max(error) * 0.8,
                f'MAE: {mae:.2f}\nMSE: {mse:.2f}',
                fontsize=11, weight='bold',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8))

    axes[1].set_xlabel('基因组位置 (bp)', fontsize=12, weight='bold')
    axes[1].set_ylabel('预测误差', fontsize=12, weight='bold')
    axes[1].set_title('预测误差分布', fontsize=14, weight='bold', pad=15)
    axes[1].grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()
    plt.savefig('5_profile_prediction.png', dpi=150, bbox_inches='tight')
    print("✅ 保存: 5_profile_prediction.png")
    plt.close()


def visualize_model_architecture():
    """可视化完整的模型架构"""

    fig, ax = plt.subplots(figsize=(12, 14))

    # 定义每个组件的位置和大小
    components = [
        # (x, y, width, height, label, color)
        (3, 13, 2, 0.8, 'Input\nDNA Sequence\n(2114, 4)', '#2E86AB'),
        (3, 11.5, 2, 0.8, 'Conv1D\nkernel=21\nfilters=64', '#A23B72'),
        (3, 10, 2, 0.8, 'Dilated Conv\nLayer 1\ndilation=2', '#F18F01'),
        (3, 8.8, 2, 0.8, 'Dilated Conv\nLayer 2-8\n...', '#C73E1D'),
        (3, 7.6, 2, 0.8, 'Dilated Conv\nLayer 9\ndilation=512', '#6A994E'),

        # 分支1: Profile
        (1, 5.5, 2, 0.8, 'Conv1D\nkernel=75\nfilters=1', '#2E86AB'),
        (1, 4, 2, 0.8, 'Cropping\nto 1000 bp', '#A23B72'),
        (1, 2.5, 2, 0.8, 'Flatten', '#F18F01'),
        (1, 1, 2, 0.8, 'Profile Output\n(1000,)', '#C73E1D'),

        # 分支2: Counts
        (5, 5.5, 2, 0.8, 'GlobalAvgPool', '#6A994E'),
        (5, 4, 2, 0.8, 'Dense\n1 unit', '#2E86AB'),
        (5, 1, 2, 0.8, 'Counts Output\n(1,)', '#A23B72'),
    ]

    # 画组件
    for x, y, w, h, label, color in components:
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                            edgecolor='black', facecolor=color,
                            linewidth=2, alpha=0.7)
        ax.add_patch(box)
        ax.text(x + w/2, y + h/2, label, ha='center', va='center',
               fontsize=9, weight='bold', color='white')

    # 画箭头
    arrows = [
        ((4, 13), (4, 12.3)),      # input -> conv1
        ((4, 11.5), (4, 10.8)),    # conv1 -> dil1
        ((4, 10), (4, 9.6)),       # dil1 -> dil2-8
        ((4, 8.8), (4, 8.4)),      # dil2-8 -> dil9
        ((4, 7.6), (2, 6.3)),      # dil9 -> profile branch
        ((4, 7.6), (6, 6.3)),      # dil9 -> counts branch
        ((2, 5.5), (2, 4.8)),      # conv75 -> crop
        ((2, 4), (2, 3.3)),        # crop -> flatten
        ((2, 2.5), (2, 1.8)),      # flatten -> output
        ((6, 5.5), (6, 4.8)),      # gap -> dense
        ((6, 4), (6, 1.8)),        # dense -> output
    ]

    for start, end in arrows:
        arrow = FancyArrowPatch(start, end, arrowstyle='->',
                               mutation_scale=15, linewidth=2,
                               color='black')
        ax.add_patch(arrow)

    # 标注分支
    ax.text(2, 6.8, 'Profile Branch', fontsize=12, weight='bold',
           color='red', ha='center',
           bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8))
    ax.text(6, 6.8, 'Counts Branch', fontsize=12, weight='bold',
           color='red', ha='center',
           bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8))

    # 标注感受野
    ax.text(5.5, 10, 'Receptive\nField:\n5 bp', fontsize=9,
           ha='left', style='italic')
    ax.text(5.5, 7.6, 'Receptive\nField:\n2000+ bp', fontsize=9,
           ha='left', style='italic')

    ax.set_xlim(0, 8)
    ax.set_ylim(0, 14.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('BPNet 完整架构', fontsize=18, weight='bold', pad=30)

    plt.tight_layout()
    plt.savefig('6_complete_architecture.png', dpi=150, bbox_inches='tight')
    print("✅ 保存: 6_complete_architecture.png")
    plt.close()


def visualize_training_curve():
    """可视化训练过程"""

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    epochs = np.arange(1, 51)

    # 模拟损失曲线
    train_loss = 5.0 * np.exp(-epochs / 10) + 0.3 + 0.1 * np.random.randn(50) * 0.1
    val_loss = 5.2 * np.exp(-epochs / 10) + 0.4 + 0.1 * np.random.randn(50) * 0.15

    # 左图：总损失
    axes[0].plot(epochs, train_loss, label='训练损失', linewidth=2, color='#2E86AB')
    axes[0].plot(epochs, val_loss, label='验证损失', linewidth=2, color='#C73E1D')
    axes[0].fill_between(epochs, train_loss, alpha=0.2, color='#2E86AB')
    axes[0].fill_between(epochs, val_loss, alpha=0.2, color='#C73E1D')

    axes[0].set_xlabel('Epoch', fontsize=12, weight='bold')
    axes[0].set_ylabel('Loss', fontsize=12, weight='bold')
    axes[0].set_title('训练过程 - 损失曲线', fontsize=14, weight='bold')
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3, linestyle='--')

    # 标注关键点
    axes[0].scatter([10, 30], [train_loss[9], train_loss[29]],
                   s=100, color='red', marker='o', zorder=5)
    axes[0].annotate('快速下降', xy=(10, train_loss[9]),
                    xytext=(15, 3), fontsize=10,
                    arrowprops=dict(arrowstyle='->', lw=1.5))
    axes[0].annotate('收敛', xy=(30, train_loss[29]),
                    xytext=(35, 1.5), fontsize=10,
                    arrowprops=dict(arrowstyle='->', lw=1.5))

    # 右图：分解损失
    profile_loss = 3.0 * np.exp(-epochs / 10) + 0.2
    counts_loss = 2.0 * np.exp(-epochs / 10) + 0.1

    axes[1].plot(epochs, profile_loss, label='Profile 损失',
                linewidth=2, color='#A23B72')
    axes[1].plot(epochs, counts_loss, label='Counts 损失',
                linewidth=2, color='#F18F01')
    axes[1].fill_between(epochs, profile_loss, alpha=0.2, color='#A23B72')
    axes[1].fill_between(epochs, counts_loss, alpha=0.2, color='#F18F01')

    axes[1].set_xlabel('Epoch', fontsize=12, weight='bold')
    axes[1].set_ylabel('Loss', fontsize=12, weight='bold')
    axes[1].set_title('训练过程 - 分解损失', fontsize=14, weight='bold')
    axes[1].legend(fontsize=11)
    axes[1].grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()
    plt.savefig('7_training_curve.png', dpi=150, bbox_inches='tight')
    print("✅ 保存: 7_training_curve.png")
    plt.close()


def main():
    """主函数：生成所有可视化"""

    print("=" * 60)
    print("🎨 BPNet 模型可视化工具")
    print("=" * 60)
    print()

    print("开始生成可视化图表...")
    print()

    # 生成所有图表
    visualize_dna_encoding()
    visualize_receptive_field()
    visualize_convolution_types()
    visualize_residual_connection()
    visualize_profile_prediction()
    visualize_model_architecture()
    visualize_training_curve()

    print()
    print("=" * 60)
    print("✅ 所有图表生成完成！")
    print("=" * 60)
    print()
    print("生成的图表:")
    print("  1️⃣  1_dna_encoding.png - DNA 序列编码")
    print("  2️⃣  2_receptive_field_growth.png - 感受野增长")
    print("  3️⃣  3_convolution_types.png - 卷积类型对比")
    print("  4️⃣  4_residual_connection.png - 残差连接")
    print("  5️⃣  5_profile_prediction.png - Profile 预测示例")
    print("  6️⃣  6_complete_architecture.png - 完整架构")
    print("  7️⃣  7_training_curve.png - 训练曲线")
    print()
    print("请查看生成的图片以更好地理解 BPNet 模型！")
    print()


if __name__ == "__main__":
    main()
