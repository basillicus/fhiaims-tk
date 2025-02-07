import argparse
import numpy as np
from sklearn.model_selection import train_test_split

from ase.io import read, write


of_train = 'train.xyz'
of_test = 'test.xyz'
if_all = 'train_all.xyz'

all_structures = read('train_all.xyz', index=':')
train_size = 100
test_size = 25

parser = argparse.ArgumentParser(
    prog='fhi_extract_samples.py',
    description='Extract n train samples  and m test samples from a dataset',
)

# # TODO: Add parameters!
# parser.add_argument('-p', '--path', default=personal_path,
#                     help='set basis set path. Change the personal_path variable in the script so you do not have to set it everytime')
# parser.add_argument('-l', '--level', default='light', choices=['light', 'lightdense', 'light_spd', 'ligh_rm2', 'intermediate', 'intermediate_gw',  'tight', 'tight_gw', 'really-tight', 'really_tight_gw',  'minimal+s', ], help='set basis set level')
# train, test = train_test_split(all_structures, train_size=train_size, test_size=test_size)
prnt('Train size:', len(train))
print('Test size:', len(test))


for geom in train:
    write(of_train, geom, append=True)
for geom in test:
    write(of_test, geom, append=True)
