#!/usr/bin/env python3

import argparse
import numpy as np
from sklearn.model_selection import train_test_split

from ase.io import read, write


of_train = 'train.xyz'
of_test = 'test.xyz'
if_all = 'train_all.xyz'

train_size = 100
test_size = 25

parser = argparse.ArgumentParser(
    prog='fhi_extract_samples.py',
    description='Extract n train samples and m test samples from a dataset',
        )

parser.add_argument('-i', '--input', default=if_all,
                    help='Input file containing all configurations')
parser.add_argument('--out_train', default=of_train,
                    help='Output file containing training configurations')
parser.add_argument('--out_test', default=of_test,
                    help='Output file containing test configurations')
parser.add_argument('--out_validation', default=None,
                    help='Output file containing validation configurations')
parser.add_argument('--n_test', default=test_size, type=int,
                    help='Number of test configurations')
parser.add_argument('--n_train', default=train_size,type=int,
                    help='Number of train configurations')
parser.add_argument('--n_validation', default=0, type=int,
                    help='Number of validation configurations')

args = parser.parse_args()

input_file = args.input
of_train = args.out_train
of_test = args.out_test
of_valid = args.out_validation

train_size = args.n_train
test_size = args.n_test
valid_size = args.n_validation


all_structures = read(input_file, index=':')

if of_valid is not None and valid_size > 0:
    train, temp = train_test_split(all_structures, train_size=train_size, test_size=test_size+valid_size)
    test, valid = train_test_split(temp, train_size=test_size, test_size=valid_size)
    for geom in valid:
        write(of_valid, geom, append=True)
else:
    # validtion file is empty
    valid=[]
    train, test = train_test_split(all_structures, train_size=train_size, test_size=test_size)

print('Train size:', len(train))
print('Test size:', len(test))
if len(valid) > 0:
    print('Valid size:', len(valid))

for geom in train:
    write(of_train, geom, append=True)
for geom in test:
    write(of_test, geom, append=True)
