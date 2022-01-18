#!/bin/sh

# This script prepares the chembl data

# data path
dir_path=$(dirname "$0")
data_path=${dir_path}/..
# downloading mols from the ChEMBL
python ${dir_path}/00_download_mols_from_chembl.py --output ${data_path}/00_chembl_dataset.sdf.gz
# convert to mae input format
python ${dir_path}/01_convert_sdf_to_mae.py --input ${data_path}/00_chembl_dataset.sdf.gz --output ${data_path}/01_chembl_dataset.mae.gz
# predict microstate pKa values with EPIK
python ${dir_path}/02_predict_pka_with_epik.py --input ${data_path}/01_chembl_dataset.mae.gz --output ${data_path}/02_chembl_dataset.mae.gz
# convert to sdf file format
python ${dir_path}/03_convert_mae_to_sdf.py --input ${data_path}/02_chembl_dataset.mae.gz --output ${data_path}/03_chembl_dataset.sdf.gz
#filter mols that are present in test sets
python ${dir_path}/04_0_filter_testmols.py --input ${data_path}/03_chembl_dataset.sdf.gz  --output ${data_path}/04_chembl_dataset_filtered.sdf.gz --filter ${data_path}/00_AvLiLuMoVe_testdata.sdf,${data_path}/00_novartis_testdata.sdf
# split mols in protonated/deprotonated pairs with pka values
python ${dir_path}/04_1_split_epik_output.py --input ${data_path}/04_chembl_dataset_filtered.sdf.gz --output ${data_path}/04_chembl_dataset_pyg.pkl
# generate pyg input data
python ${dir_path}/05_data_preprocess.py --input ${data_path}/04_chembl_dataset_pyg.pkl --output ${data_path}/05_chembl_dataset_pyg.pkl