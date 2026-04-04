#!/usr/bin/env python3
"""
Genos-1.2B Feasibility Probe for TransChromBP Integration
==========================================================

Phase 0 验证脚本：回答三个核心问题
  1. Genos 双向 embedding 能否区分 peak vs nonpeak？
  2. 不同层的 embedding 哪些对 ATAC 区域最有区分度？
  3. 实际训练场景下的显存和吞吐是多少？

Usage:
    # 在 6000 服务器上运行（需 GPU）
    export CUDA_VISIBLE_DEVICES=1  # 用空闲卡，避免影响正在跑的 v2fix
    /data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b/bin/python scripts/genos_feasibility.py

    # 快速模式（仅 50 个区域，用于验证脚本能跑通）
    QUICK=1 /data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b/bin/python scripts/genos_feasibility.py
"""

import os
import json
import time
import argparse
import resource
import numpy as np

import torch
from transformers import AutoTokenizer, AutoModelForCausalLM


# ─── 配置 ─────────────────────────────────────────────────
GENOS_MODEL_PATH = "/data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos-1.2B"
GENOME_FASTA = "/data1/zhoujiazhen/bylw_atac/hg38.fa"
PEAKS_BED = "/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.peaks.bed"
NONPEAKS_BED = "/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.nonpeaks.bed"
FOLDS_JSON = "/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json"
INPUT_LEN = 2114
OUTPUT_DIR = "/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/genos_feasibility"

QUICK_MODE = os.environ.get("QUICK", "0") == "1"
N_SAMPLES = 50 if QUICK_MODE else 500  # 每类的采样数
BATCH_SIZE = 16
DEFAULT_TEST_SIZE = 0.2
DEFAULT_CV_FOLDS = 5

RC_MAP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
           "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}


def reverse_complement(seq: str) -> str:
    return "".join(RC_MAP.get(c, "N") for c in reversed(seq))


# ─── 数据加载 ──────────────────────────────────────────────
def load_folds(folds_json: str):
    with open(folds_json) as f:
        folds = json.load(f)
    # 取 validation 染色体用于验证（避免 train/test 偏差）
    val_chroms = set(folds.get("valid", folds.get("val", [])))
    if not val_chroms:
        val_chroms = set(folds.get("0", {}).get("valid", []))
    return val_chroms


def load_bed_regions(bed_path: str, val_chroms: set, n_samples: int, seed: int = 42):
    """从 BED 文件加载区域，只保留 validation 染色体，随机采样 n_samples 个。"""
    regions = []
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom = parts[0]
            if val_chroms and chrom not in val_chroms:
                continue
            start = int(parts[1])
            end = int(parts[2])
            center = (start + end) // 2
            # 以 center 为中心取 INPUT_LEN 窗口
            win_start = center - INPUT_LEN // 2
            win_end = win_start + INPUT_LEN
            if win_start >= 0:
                regions.append((chrom, win_start, win_end))
    rng = np.random.default_rng(seed)
    if len(regions) > n_samples:
        idx = rng.choice(len(regions), n_samples, replace=False)
        regions = [regions[i] for i in idx]
    return regions


def extract_sequences(regions, genome_fasta):
    """用 pysam 从 FASTA 提取序列。"""
    try:
        import pysam
    except ImportError:
        print("pysam not available, trying pyfaidx...")
        return extract_sequences_pyfaidx(regions, genome_fasta)

    fa = pysam.FastaFile(genome_fasta)
    seqs = []
    for chrom, start, end in regions:
        seq = fa.fetch(chrom, start, end).upper()
        if len(seq) < INPUT_LEN:
            seq = seq + "N" * (INPUT_LEN - len(seq))
        seqs.append(seq)
    fa.close()
    return seqs


def extract_sequences_pyfaidx(regions, genome_fasta):
    """pyfaidx 备选方案。"""
    from pyfaidx import Fasta
    fa = Fasta(genome_fasta)
    seqs = []
    for chrom, start, end in regions:
        seq = str(fa[chrom][start:end]).upper()
        if len(seq) < INPUT_LEN:
            seq = seq + "N" * (INPUT_LEN - len(seq))
        seqs.append(seq)
    return seqs


# ─── Genos 推理 ───────────────────────────────────────────
def load_genos(model_path: str, device: str):
    """加载 Genos-1.2B 模型和 tokenizer。"""
    print(f"Loading Genos from {model_path}...")
    tokenizer = AutoTokenizer.from_pretrained(model_path, local_files_only=True)
    model = AutoModelForCausalLM.from_pretrained(
        model_path,
        local_files_only=True,
        dtype=torch.bfloat16,
        attn_implementation="flash_attention_2",
    ).to(device).eval()
    print(f"Genos loaded. Device: {device}")
    return tokenizer, model


def tokenize_batch(tokenizer, sequences, device):
    """Tokenize DNA sequences without special tokens to preserve 1:1 base alignment."""
    tokens = tokenizer(
        sequences,
        return_tensors="pt",
        padding=True,
        add_special_tokens=False,
    )
    input_ids = tokens["input_ids"].to(device)
    attention_mask = tokens.get("attention_mask")
    if attention_mask is not None:
        attention_mask = attention_mask.to(device)
    if input_ids.shape[1] != INPUT_LEN:
        raise ValueError(
            f"Tokenized length mismatch: expected {INPUT_LEN}, got {input_ids.shape[1]}. "
            "This breaks base-wise alignment."
        )
    return {"input_ids": input_ids, "attention_mask": attention_mask}


@torch.no_grad()
def collect_layer_statistics(
    model,
    tokenizer,
    peak_sequences,
    nonpeak_sequences,
    device,
    batch_size=16,
    layers=None,
):
    """流式提取多层双向 embedding，并直接累计可行性分析所需统计量。"""
    if layers is None:
        layers = [-1]

    center_start = (INPUT_LEN - 1000) // 2
    center_end = center_start + 1000
    hidden_size = int(model.config.hidden_size)
    num_hidden_layers = int(model.config.num_hidden_layers)

    layer_stats = {}
    for l in layers:
        layer_name = f"layer_{l}" if l >= 0 else f"layer_{num_hidden_layers + l + 1}"
        layer_stats[l] = {
            "layer": layer_name,
            "peak_global": [],
            "peak_center": [],
            "nonpeak_global": [],
            "nonpeak_center": [],
            "peak_center_norm_sum": 0.0,
            "nonpeak_center_norm_sum": 0.0,
            "peak_count": 0,
            "nonpeak_count": 0,
            "var_count": 0,
            "var_mean": np.zeros((INPUT_LEN, hidden_size), dtype=np.float32),
            "var_m2": np.zeros((INPUT_LEN, hidden_size), dtype=np.float32),
        }

    def update_running_variance(state, batch_emb):
        batch_count = batch_emb.shape[0]
        batch_mean = batch_emb.mean(axis=0, dtype=np.float32)
        centered = batch_emb - batch_mean[None, :, :]
        np.square(centered, out=centered)
        batch_m2 = centered.sum(axis=0, dtype=np.float32)

        if state["var_count"] == 0:
            state["var_mean"][:] = batch_mean
            state["var_m2"][:] = batch_m2
            state["var_count"] = batch_count
            return

        prev_count = state["var_count"]
        total_count = prev_count + batch_count
        delta = batch_mean - state["var_mean"]
        state["var_mean"] += delta * (batch_count / total_count)
        state["var_m2"] += batch_m2 + np.square(delta) * (prev_count * batch_count / total_count)
        state["var_count"] = total_count

    def accumulate_batch_stats(state, batch_emb, class_name):
        global_mean = batch_emb.mean(axis=1, dtype=np.float32)
        center_mean = batch_emb[:, center_start:center_end, :].mean(axis=1, dtype=np.float32)

        state[f"{class_name}_global"].append(global_mean)
        state[f"{class_name}_center"].append(center_mean)
        state[f"{class_name}_center_norm_sum"] += float(np.linalg.norm(center_mean, axis=-1).sum())
        state[f"{class_name}_count"] += batch_emb.shape[0]
        update_running_variance(state, batch_emb)

    def process_sequences(sequences, class_name):
        n = len(sequences)
        for i in range(0, n, batch_size):
            batch_seqs = sequences[i : i + batch_size]
            batch_rc = [reverse_complement(s) for s in batch_seqs]

            fwd_tokens = tokenize_batch(tokenizer, batch_seqs, device)
            fwd_out = model(
                input_ids=fwd_tokens["input_ids"],
                attention_mask=fwd_tokens["attention_mask"],
                output_hidden_states=True,
                use_cache=False,
            )

            rev_tokens = tokenize_batch(tokenizer, batch_rc, device)
            rev_out = model(
                input_ids=rev_tokens["input_ids"],
                attention_mask=rev_tokens["attention_mask"],
                output_hidden_states=True,
                use_cache=False,
            )

            for l in layers:
                fwd_h = fwd_out.hidden_states[l][:, :INPUT_LEN, :]
                rev_h = rev_out.hidden_states[l][:, :INPUT_LEN, :].flip(1)
                bidir = ((fwd_h + rev_h) / 2.0).float().cpu().numpy()
                accumulate_batch_stats(layer_stats[l], bidir, class_name)

            del fwd_out, rev_out

            if (i // batch_size) % 10 == 0:
                print(f"  [{class_name}] Processed {min(i + batch_size, n)}/{n} sequences")

    process_sequences(peak_sequences, "peak")
    process_sequences(nonpeak_sequences, "nonpeak")
    return layer_stats


# ─── 分析 ─────────────────────────────────────────────────
def evaluate_linear_probe(X, y, seed=42, test_size=DEFAULT_TEST_SIZE, cv_folds=DEFAULT_CV_FOLDS):
    """Evaluate a linear probe with hold-out and stratified cross-validation."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, roc_auc_score
    from sklearn.model_selection import StratifiedKFold, train_test_split
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=test_size,
        random_state=seed,
        stratify=y,
    )

    holdout_model = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=1000, C=1.0, random_state=seed),
    )
    holdout_model.fit(X_train, y_train)
    holdout_pred = holdout_model.predict(X_test)
    holdout_prob = holdout_model.predict_proba(X_test)[:, 1]
    holdout_acc = accuracy_score(y_test, holdout_pred)
    holdout_auc = roc_auc_score(y_test, holdout_prob)

    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    n_splits = min(cv_folds, n_pos, n_neg)
    if n_splits < 2:
        raise ValueError(
            f"Not enough samples for stratified CV: positives={n_pos}, negatives={n_neg}, n_splits={n_splits}"
        )
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    cv_acc_scores = []
    cv_auc_scores = []
    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        fold_model = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=1000, C=1.0, random_state=seed + fold_idx),
        )
        fold_model.fit(X[train_idx], y[train_idx])
        fold_pred = fold_model.predict(X[test_idx])
        fold_prob = fold_model.predict_proba(X[test_idx])[:, 1]
        cv_acc_scores.append(accuracy_score(y[test_idx], fold_pred))
        cv_auc_scores.append(roc_auc_score(y[test_idx], fold_prob))

    return {
        "holdout_acc": float(holdout_acc),
        "holdout_auc": float(holdout_auc),
        "cv_acc_mean": float(np.mean(cv_acc_scores)),
        "cv_acc_std": float(np.std(cv_acc_scores)),
        "cv_auc_mean": float(np.mean(cv_auc_scores)),
        "cv_auc_std": float(np.std(cv_auc_scores)),
        "n_splits": int(n_splits),
        "test_size": float(test_size),
    }


def analyze_discriminability(layer_state, seed=42, test_size=DEFAULT_TEST_SIZE, cv_folds=DEFAULT_CV_FOLDS):
    """
    分析 peak vs nonpeak 的区分度。

    使用三个指标：
    1. Hold-out 线性探针（主证据）
    2. Stratified CV 线性探针（稳健性）
    3. Center window vs flanking 的信号比
    """
    layer_name = layer_state["layer"]
    peak_global = np.concatenate(layer_state["peak_global"], axis=0)
    nonpeak_global = np.concatenate(layer_state["nonpeak_global"], axis=0)
    peak_center = np.concatenate(layer_state["peak_center"], axis=0)
    nonpeak_center = np.concatenate(layer_state["nonpeak_center"], axis=0)

    print(f"\n{'='*60}")
    print(f"Layer: {layer_name}")
    print(f"  Peak pooled features: global={peak_global.shape}, center={peak_center.shape}")
    print(f"  Nonpeak pooled features: global={nonpeak_global.shape}, center={nonpeak_center.shape}")

    # 2. 线性探针
    X_global = np.vstack([peak_global, nonpeak_global])
    X_center = np.vstack([peak_center, nonpeak_center])
    y = np.array([1] * len(peak_global) + [0] * len(nonpeak_global))

    pooled_results = {}
    for name, X in [("global_mean", X_global), ("center_mean", X_center)]:
        pooled_results[name] = evaluate_linear_probe(
            X,
            y,
            seed=seed,
            test_size=test_size,
            cv_folds=cv_folds,
        )
        result = pooled_results[name]
        print(
            f"  {name} holdout: acc={result['holdout_acc']:.4f}  auc={result['holdout_auc']:.4f}"
        )
        print(
            f"  {name} {result['n_splits']}-fold CV: "
            f"acc={result['cv_acc_mean']:.4f}±{result['cv_acc_std']:.4f}  "
            f"auc={result['cv_auc_mean']:.4f}±{result['cv_auc_std']:.4f}"
        )

    # 3. Embedding norm 对比
    peak_norm = layer_state["peak_center_norm_sum"] / max(layer_state["peak_count"], 1)
    nonpeak_norm = layer_state["nonpeak_center_norm_sum"] / max(layer_state["nonpeak_count"], 1)
    print(f"  Center embedding L2 norm: peak={peak_norm:.4f}  nonpeak={nonpeak_norm:.4f}")

    # 4. 位置级方差（看哪些位置区分度高）
    center_start = (INPUT_LEN - 1000) // 2
    center_end = center_start + 1000
    per_pos_var = (layer_state["var_m2"] / max(layer_state["var_count"], 1)).mean(axis=-1)
    center_var = per_pos_var[center_start:center_start+1000].mean()
    flank_var = np.concatenate([
        per_pos_var[:center_start],
        per_pos_var[center_end:]
    ]).mean()
    center_var_ratio = center_var / flank_var if flank_var > 0 else float("nan")
    print(f"  Positional variance: center={center_var:.6f}  flank={flank_var:.6f}  ratio={center_var_ratio:.2f}")

    global_result = pooled_results["global_mean"]
    center_result = pooled_results["center_mean"]
    return {
        "layer": layer_name,
        "global_holdout_acc": global_result["holdout_acc"],
        "global_holdout_auc": global_result["holdout_auc"],
        "global_cv_acc_mean": global_result["cv_acc_mean"],
        "global_cv_acc_std": global_result["cv_acc_std"],
        "global_cv_auc_mean": global_result["cv_auc_mean"],
        "global_cv_auc_std": global_result["cv_auc_std"],
        "center_holdout_acc": center_result["holdout_acc"],
        "center_holdout_auc": center_result["holdout_auc"],
        "center_cv_acc_mean": center_result["cv_acc_mean"],
        "center_cv_acc_std": center_result["cv_acc_std"],
        "center_cv_auc_mean": center_result["cv_auc_mean"],
        "center_cv_auc_std": center_result["cv_auc_std"],
        "center_var_ratio": float(center_var_ratio),
        "test_size": float(test_size),
        "cv_folds": int(center_result["n_splits"]),
    }


def get_peak_rss_gb():
    """获取当前进程的峰值 RSS（Linux 上 ru_maxrss 以 KB 计）。"""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**2)


def benchmark_throughput(model, tokenizer, device, batch_sizes=(4, 8, 16)):
    """测量不同 batch size 下双向推理的吞吐和显存。"""
    print(f"\n{'='*60}")
    print("Throughput Benchmark (bidirectional, 2114bp)")
    print(f"{'batch':>6} {'fwd_time':>10} {'samples/s':>10} {'peak_mem':>10}")

    dummy_seq = "ACGT" * (INPUT_LEN // 4) + "AC"
    results = []

    for bs in batch_sizes:
        seqs = [dummy_seq] * bs
        rc_seqs = [reverse_complement(s) for s in seqs]

        try:
            if torch.cuda.is_available():
                torch.cuda.reset_peak_memory_stats(device)

            # Warmup
            warmup_tokens = tokenize_batch(tokenizer, seqs, device)
            _ = model(
                input_ids=warmup_tokens["input_ids"],
                attention_mask=warmup_tokens["attention_mask"],
                use_cache=False,
            )

            # Timed run (fwd + rev)
            t0 = time.time()
            for _ in range(3):
                fwd_tokens = tokenize_batch(tokenizer, seqs, device)
                model(
                    input_ids=fwd_tokens["input_ids"],
                    attention_mask=fwd_tokens["attention_mask"],
                    use_cache=False,
                )
                rev_tokens = tokenize_batch(tokenizer, rc_seqs, device)
                model(
                    input_ids=rev_tokens["input_ids"],
                    attention_mask=rev_tokens["attention_mask"],
                    use_cache=False,
                )
            if torch.cuda.is_available():
                torch.cuda.synchronize()
            elapsed = (time.time() - t0) / 3.0  # avg per iteration (fwd+rev)

            peak_mem = 0.0
            if torch.cuda.is_available():
                peak_mem = torch.cuda.max_memory_allocated(device) / 1024**3
            throughput = bs / elapsed

            print(f"{bs:>6} {elapsed:>10.3f}s {throughput:>10.1f} {peak_mem:>9.2f}GB")
            results.append({
                "batch_size": bs,
                "status": "ok",
                "bidir_time_s": round(elapsed, 3),
                "samples_per_s": round(throughput, 1),
                "peak_mem_gb": round(peak_mem, 2),
            })
        except RuntimeError as exc:
            if "out of memory" not in str(exc).lower():
                raise
            print(f"{bs:>6} {'OOM':>10} {'-':>10} {'-':>10}")
            results.append({
                "batch_size": bs,
                "status": "oom",
            })
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    return results


# ─── 主流程 ───────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Genos feasibility probe")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--n-samples", type=int, default=N_SAMPLES)
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    parser.add_argument("--test-size", type=float, default=DEFAULT_TEST_SIZE)
    parser.add_argument("--cv-folds", type=int, default=DEFAULT_CV_FOLDS)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", default=OUTPUT_DIR)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    device = args.device

    # 1. 加载 Genos
    tokenizer, model = load_genos(GENOS_MODEL_PATH, device)

    # 2. 加载数据
    print("\nLoading genomic regions...")
    val_chroms = load_folds(FOLDS_JSON)
    print(f"Validation chromosomes: {sorted(val_chroms)}")

    peak_regions = load_bed_regions(PEAKS_BED, val_chroms, args.n_samples, seed=args.seed)
    nonpeak_regions = load_bed_regions(NONPEAKS_BED, val_chroms, args.n_samples, seed=args.seed + 1)
    print(f"Peak regions: {len(peak_regions)}, Nonpeak regions: {len(nonpeak_regions)}")

    if len(peak_regions) == 0 or len(nonpeak_regions) == 0:
        print("ERROR: No regions found in validation chromosomes.")
        print("Trying all chromosomes as fallback...")
        peak_regions = load_bed_regions(PEAKS_BED, set(), args.n_samples, seed=args.seed)
        nonpeak_regions = load_bed_regions(NONPEAKS_BED, set(), args.n_samples, seed=args.seed + 1)
        print(f"Fallback peak regions: {len(peak_regions)}, nonpeak regions: {len(nonpeak_regions)}")

    print("Extracting sequences from genome...")
    peak_seqs = extract_sequences(peak_regions, GENOME_FASTA)
    nonpeak_seqs = extract_sequences(nonpeak_regions, GENOME_FASTA)
    print(f"Extracted: {len(peak_seqs)} peak, {len(nonpeak_seqs)} nonpeak sequences")

    # 3. 提取多层 embedding
    layers_to_check = [0, 3, 6, 9, -1]  # 第 0, 3, 6, 9, 12(last) 层
    print(f"\nExtracting bidirectional embeddings for layers {layers_to_check}...")

    t0 = time.time()
    layer_stats = collect_layer_statistics(
        model,
        tokenizer,
        peak_seqs,
        nonpeak_seqs,
        device,
        args.batch_size,
        layers_to_check,
    )
    extraction_time = time.time() - t0
    print(f"Extraction time: {extraction_time:.1f}s "
          f"({len(peak_seqs) + len(nonpeak_seqs)} sequences, "
          f"{2 * (len(peak_seqs) + len(nonpeak_seqs))} forward passes)")
    print(f"Peak host RSS so far: {get_peak_rss_gb():.2f} GB")

    # 4. 逐层分析
    all_results = []
    for l in layers_to_check:
        r = analyze_discriminability(
            layer_stats[l],
            seed=args.seed,
            test_size=args.test_size,
            cv_folds=args.cv_folds,
        )
        all_results.append(r)

    # 5. 吞吐 benchmark
    bench_results = benchmark_throughput(model, tokenizer, device)

    # 6. 保存结果
    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "n_peak": len(peak_seqs),
        "n_nonpeak": len(nonpeak_seqs),
        "input_len": INPUT_LEN,
        "extraction_time_s": round(extraction_time, 1),
        "quick_mode": QUICK_MODE,
        "batch_size": args.batch_size,
        "test_size": args.test_size,
        "cv_folds": args.cv_folds,
        "selection_metric": "center_cv_auc_mean",
        "analysis_mode": "streaming_pooled_stats",
        "peak_host_rss_gb": round(get_peak_rss_gb(), 2),
        "layer_analysis": all_results,
        "throughput_benchmark": bench_results,
    }

    report_path = os.path.join(args.output_dir, "feasibility_report.json")
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nReport saved to {report_path}")

    # 7. 打印摘要
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Samples: {len(peak_seqs)} peaks + {len(nonpeak_seqs)} nonpeaks")
    print(f"Extraction time: {extraction_time:.1f}s")
    print()
    print(
        f"{'Layer':<12} {'Global CV AUC':>14} {'Center CV AUC':>14} "
        f"{'Center Holdout':>16} {'Var Ratio':>12}"
    )
    for r in all_results:
        print(
            f"{r['layer']:<12} {r['global_cv_auc_mean']:>14.4f} "
            f"{r['center_cv_auc_mean']:>14.4f} "
            f"{r['center_holdout_auc']:>16.4f} {r['center_var_ratio']:>12.2f}"
        )
    print()
    print("Throughput (bidirectional):")
    for b in bench_results:
        if b.get("status") == "oom":
            print(f"  batch={b['batch_size']:>2}: OOM")
        else:
            print(
                f"  batch={b['batch_size']:>2}: {b['samples_per_s']:.1f} samples/s, "
                f"{b['peak_mem_gb']:.2f} GB"
            )

    # 8. 结论判断
    best_layer = max(all_results, key=lambda r: r["center_cv_auc_mean"])
    print(f"\nBest layer for peak/nonpeak discrimination: {best_layer['layer']} "
          f"(center CV AUC={best_layer['center_cv_auc_mean']:.4f})")

    if best_layer["center_cv_auc_mean"] > 0.7:
        print("✓ Genos embeddings show meaningful peak/nonpeak discriminability.")
        print("  → Phase 1 (dual-branch fusion) is worth pursuing.")
    elif best_layer["center_cv_auc_mean"] > 0.55:
        print("△ Genos embeddings show weak discriminability.")
        print("  → Proceed with caution; may need LoRA or different fusion strategy.")
    else:
        print("✗ Genos embeddings do not discriminate peak/nonpeak.")
        print("  → Foundation model features may not help for this task.")


if __name__ == "__main__":
    main()
