import os

import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer


MODEL_DIR = "/data1/zhoujiazhen/bylw_atac/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species"


def main() -> None:
	tokenizer = AutoTokenizer.from_pretrained(
		MODEL_DIR,
		trust_remote_code=True,
		local_files_only=True,
	)

	model = AutoModelForMaskedLM.from_pretrained(
		MODEL_DIR,
		trust_remote_code=True,
		local_files_only=True,
	)
	model.eval()

	sequences = [
		"ATTCCGATTCCGATTCCG",
		"ATTTCTCTCTCTCTCTGAGATCGATCGATCGAT",
	]

	max_length = min(getattr(tokenizer, "model_max_length", 256), 256)
	batch = tokenizer.batch_encode_plus(
		sequences,
		return_tensors="pt",
		padding="max_length",
		truncation=True,
		max_length=max_length,
	)

	attention_mask = batch["attention_mask"].bool()

	with torch.no_grad():
		outputs = model(
			batch["input_ids"],
			attention_mask=attention_mask,
			encoder_attention_mask=attention_mask,
			output_hidden_states=True,
		)

	hidden = outputs["hidden_states"][-1]
	mask = attention_mask.unsqueeze(-1)
	mean_emb = (hidden * mask).sum(dim=1) / mask.sum(dim=1)

	print("model_dir:", MODEL_DIR)
	print("torch device:", "cpu")
	print("input_ids shape:", tuple(batch["input_ids"].shape))
	print("attention_mask shape:", tuple(attention_mask.shape))
	print("last_hidden_state shape:", tuple(hidden.shape))
	print("mean_sequence_embedding shape:", tuple(mean_emb.shape))


if __name__ == "__main__":
	main()
