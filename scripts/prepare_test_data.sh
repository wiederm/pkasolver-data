#!/bin/sh
# This script prepares the data for the test sets and fine tuning training set

# data path
data_path='../'

# split mols in protonated/deprotonated pairs with pka values
python 04_2_prepare_rest.py --input ${data_path}/Baltruschat/00_novartis_testdata.sdf --output ${data_path}/04_novartis_testdata_mols.pkl
# generate pyg input data
python 05_data_preprocess.py --input ${data_path}/04_novartis_testdata_mols.pkl --output ${data_path}/05_novartis_testdata_pyg_data.pkl 
# split mols in protonated/deprotonated pairs with pka values
python 04_2_prepare_rest.py --input ${data_path}/Baltruschat/00_AvLiLuMoVe_testdata.sdf --output ${data_path}/04_AvLiLuMoVe_testdata_mols.pkl
# generate pyg input data
python 05_data_preprocess.py --input ${data_path}/04_AvLiLuMoVe_testdata_mols.pkl --output ${data_path}/05_AvLiLuMoVe_testdata_pyg_data.pkl 
# split mols in protonated/deprotonated pairs with pka values
python 04_2_prepare_rest.py --input ${data_path}/Baltruschat/00_experimental_training_datasets.sdf --output ${data_path}/04_experimental_training_dataset.pkl
# generate pyg input data
python 05_data_preprocess.py --input ${data_path}/04_experimental_training_dataset.pkl --output ${data_path}/05_experimental_training_dataset_pyg.pkl 
