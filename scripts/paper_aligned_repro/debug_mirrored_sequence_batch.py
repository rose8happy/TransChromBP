#!/usr/bin/env python3

import argparse
import math

import numpy as np
import tensorflow as tf
from tensorflow import keras


class DummySequence(keras.utils.Sequence):
    def __init__(self, batch_size: int, seq_len: int, channels: int, output_dim: int, num_batches: int):
        self.batch_size = batch_size
        self.seq_len = seq_len
        self.channels = channels
        self.output_dim = output_dim
        self.num_batches = num_batches

    def __len__(self):
        return self.num_batches

    def __getitem__(self, idx):
        del idx
        x = np.zeros((self.batch_size, self.seq_len, self.channels), dtype=np.float32)
        y_profile = np.zeros((self.batch_size, self.output_dim), dtype=np.float32)
        y_count = np.zeros((self.batch_size, 1), dtype=np.float32)
        return x, [y_profile, y_count]


class ReplicaShapeProbe(keras.layers.Layer):
    def call(self, inputs):
        replica_ctx = tf.distribute.get_replica_context()
        replica_id = replica_ctx.replica_id_in_sync_group if replica_ctx is not None else tf.constant(-1, dtype=tf.int32)
        tf.print("[batch-debug] replica", replica_id, "input_shape", tf.shape(inputs))
        return inputs


def build_model(seq_len: int, channels: int, output_dim: int):
    inp = keras.Input(shape=(seq_len, channels), name="sequence")
    x = ReplicaShapeProbe(name="replica_shape_probe")(inp)
    x = keras.layers.Conv1D(8, kernel_size=3, padding="same", activation="relu")(x)
    x = keras.layers.GlobalAveragePooling1D()(x)
    profile = keras.layers.Dense(output_dim, name="profile")(x)
    count = keras.layers.Dense(1, name="count")(x)
    model = keras.Model(inputs=inp, outputs=[profile, count])
    model.compile(optimizer="adam", loss=["mse", "mse"])
    return model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--seq-len", type=int, default=2114)
    parser.add_argument("--channels", type=int, default=4)
    parser.add_argument("--output-dim", type=int, default=1000)
    parser.add_argument("--num-batches", type=int, default=1)
    args = parser.parse_args()

    gpus = tf.config.list_physical_devices("GPU")
    print("Visible GPUs:", len(gpus))
    for i, gpu in enumerate(gpus):
        print(f"  gpu[{i}] = {gpu}")

    strategy = None
    if len(gpus) > 1:
        strategy = tf.distribute.MirroredStrategy()
        print("Using MirroredStrategy with", strategy.num_replicas_in_sync, "replicas")

    seq = DummySequence(
        batch_size=args.batch_size,
        seq_len=args.seq_len,
        channels=args.channels,
        output_dim=args.output_dim,
        num_batches=args.num_batches,
    )
    print("Sequence batch_size =", args.batch_size)
    print("Sequence __len__ =", len(seq))

    if strategy is None:
        model = build_model(args.seq_len, args.channels, args.output_dim)
    else:
        with strategy.scope():
            model = build_model(args.seq_len, args.channels, args.output_dim)

    history = model.fit(seq, epochs=1, steps_per_epoch=min(1, len(seq)), verbose=1)
    print("fit history keys:", sorted(history.history.keys()))


if __name__ == "__main__":
    main()
