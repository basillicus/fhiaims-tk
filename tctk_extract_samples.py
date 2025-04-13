#!/usr/bin/env python3
"""
Extract training, test (and optionally validation) samples from a dataset.

This script reads all configurations from an input file (either in a NumPy or ASE-readable format)
and randomly splits them into training, test, and validation sets, writing each to a separate file.

Command-line arguments:
    - -i, --input: Input file containing all configurations (default "train_all.xyz").
    - --out_train: Output file for training configurations (default "train.xyz").
    - --out_test: Output file for test configurations (default "test.xyz").
    - --out_validation: Output file for validation configurations (default None).
    - --n_train: Number of training samples to extract (default 10).
    - --n_test: Number of test samples to extract (default 5).
    - --n_validation: Number of validation samples (default 0).

"""


import argparse
import numpy as np
from sklearn.model_selection import train_test_split

from ase.io import read, write


of_train = "train.xyz"
of_test = "test.xyz"
if_all = "train_all.xyz"

train_size = 10
test_size = 5

# Note on argparse: Ensure that the number of samples is consistent with your dataset size;
parser = argparse.ArgumentParser(
    prog="fhi_extract_samples.py",
    description="Extract n train samples and m test samples from a dataset",
)

parser.add_argument(
    "-i", "--input", default=if_all, help="Input file containing all configurations. Could be a structured numpy array or any valid format that ASE can read"
)
parser.add_argument(
    "--out_train",
    default=of_train,
    help="Output file containing training configurations. Format will be guessed from the extension filename",
)
parser.add_argument(
    "--out_test", default=of_test, help="Output file containing test configurations"
)
parser.add_argument(
    "--out_validation",
    default=None,
    help="Output file containing validation configurations",
)
parser.add_argument(
    "--n_test", default=test_size, type=int, help="Number of test configurations"
)
parser.add_argument(
    "--n_train", default=train_size, type=int, help="Number of train configurations"
)
parser.add_argument(
    "--n_validation", default=0, type=int, help="Number of validation configurations"
)

args = parser.parse_args()

input_file = args.input
of_train = args.out_train
of_test = args.out_test
of_valid = args.out_validation

train_size = args.n_train
test_size = args.n_test
valid_size = args.n_validation


file_extension = input_file.split(".")[-1]

if file_extension == "npy" or file_extension == "npz":
    all_structures = np.load(input_file)
else:  # Assume is a input format ASE can read
    all_structures = read(input_file, index=":")

if of_valid is not None and valid_size > 0:
    train, temp = train_test_split(
        all_structures, train_size=train_size, test_size=test_size + valid_size
    )
    test, valid = train_test_split(temp, train_size=test_size, test_size=valid_size)
    for geom in valid:
        write(of_valid, geom, append=True)
else:
    # validtion file is empty
    valid = []
    train, test = train_test_split(
        all_structures, train_size=train_size, test_size=test_size
    )

print("Train size:", len(train))
print("Test size:", len(test))
if len(valid) > 0:
    print("Valid size:", len(valid))

for geom in train:
    write(of_train, geom, append=True)
for geom in test:
    write(of_test, geom, append=True)
