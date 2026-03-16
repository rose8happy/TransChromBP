# Bias Model Training Analysis (ChromBPNet Tutorial)

## Run Context
- Dataset root: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial`
- Output dir: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu_retry2.log`
- GPU: NVIDIA RTX A6000, cuDNN 8.1 (log shows device creation and cuDNN load)
- Model: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias/models/bias.h5`

## Data + Preprocessing
- Shift estimation: +0/+0
- BigWig non-zero entries (chr2): 21,159,500
- Nonpeaks input 477,258 -> after edge filter 434,778 used
- Peaks input 238,629 -> after edge filter 238,629 used
- Upper bound counts cutoff: 102.5
- Nonpeaks after cutoff/outlier removal: 221,738 -> 208,155
- Split sizes from log:
  - Train nonpeaks regions: 187,151
  - Valid nonpeaks regions: 21,004
- Training sampled 30,000 examples (log shows "Done N examples of 30000")

## Training Configuration
- Input len 2114, output len 1000
- Filters 128, dilation layers 4
- Batch size 64, learning rate 0.001
- Epochs 50, early_stop 5
- Counts loss weight 5.7

## Training Dynamics
- 20 epochs completed (early stopping triggered)
- Loss trend:
  - Train loss: 157.84 -> 154.06
  - Val loss: 160.61 -> 158.69 (best 158.64 at epoch 14)
  - Train logcount loss: 0.676 -> 0.419
  - Val logcount loss: 0.503 -> 0.314
  - Train profile loss: 153.99 -> 151.68
  - Val profile loss: 157.75 -> 156.90
- Training and validation losses decrease and then plateau, with no obvious overfitting.

## Evaluation Metrics (bias_metrics.json)
- Counts metrics (pearson):
  - peaks_and_nonpeaks: 0.296
  - nonpeaks: 0.532
  - peaks: -0.074
- Profile metrics (median JSD):
  - peaks_and_nonpeaks: 0.4726 (norm 0.2890)
  - nonpeaks: 0.5570 (norm 0.2246)
  - peaks: 0.3931 (norm 0.3651)
- Interpretation: the bias model captures background structure better than peaks, which is expected for a bias-only model.

## Key Artifacts
- Reports: `evaluation/overall_report.html`, `evaluation/overall_report.pdf`
- Curves: `evaluation/epoch_loss.png`
- Metrics: `evaluation/bias_metrics.json`
- Predictions: `evaluation/bias_predictions.h5`

## Local Artifacts (this repo)
- Training curve plot: `bias_training_curve.svg`
- Local log copy used to plot: `bias_training_log.csv`
