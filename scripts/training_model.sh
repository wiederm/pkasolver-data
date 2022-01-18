#!/bin/sh
# This script performs initial pre-training on the 
# ChEMBL dataset and then fine tuning on the experimental dataset

# path to scirpts
run=0
dir_path=$(dirname "$0")
data_path=${dir_path}/..

# start with pretraining on the CHEMBL data
python ${dir_path}/06_training.py --input ${data_path}/05_chembl_dataset_pyg.pkl --path ${data_path}/trained_models/training_run_${run} --epoch 1000
# transfer learning on the experimental data
python ${dir_path}/06_training.py --input ${data_path}/05_experimental_training_datasets_pyg.pkl --path ${data_path}/trained_models/training_run_${run} -r --epoch 1000 --reg ${data_path}/05_chembl_dataset_pyg.pkl
