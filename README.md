Scripts and notebooks shown here are used to generate the dataset and plots for the work in [pkasolver]

# Data generation

The scripts in the `scripts` folder can be used to generate the full dataset.

1. `prepare_chembl_data.sh` --contains the pipeline to download molecular data from the CHEMBL-database, predict its pka-data with Schrödinger's Epik and finally convert it the pytorch geometric graph data, to be used in the initial training.

2. `prepare_test_data.sh` -- prepares the Novartis and Literature test data.

4. `training_model.sh` --script for training and fine tuning a model with the ChEMBL dataset and the experimental dataset.

# notebook

The `plotting.ipynb` Jupyter notebook can be used to generate all plots shown in `/plots`. 
To do so you first have to unzip `04_chembl_dataset_filtered.sdf.gz` and then run 
`python 04_1_split_epik_output.py --input 04_chembl_dataset_filtered.sdf.gz --output 04_chembl_dataset_pyg.pkl`
and
`python 05_data_preprocess.py --input 04_chembl_dataset_pyg.pkl --output 05_chembl_dataset_pyg.pkl`.
This generates all the data needed to rerun the notebook.

# python scripts

`00_download_mols_from_chembl.py`:
--input: None 
--output: path to output file (sdf.gz, sdf) 

--filters the molecules of the chembl database by the specified criteria (e.g. max number of rule of five violation = 1) and outputs them to a gzipped sdf file.

`01_convert_sdf_to_mae.py` 
--input: path to input file (sdf.gz, sdf)
--output: path to output file (mae.gz, mae)

--takes sdf file (can be gzipped) and converts it to Schrödinger maestro (mae) file. Schrödinger "ligprep" CLI-binary must be installed and path must be specified inside the script.  

`02_predict_pka_with_epik.py` 
--input: path to input file (mae.gz, mae)
--output: path to output file (mae.gz, mae)

--takes molecules from Schrödinger maestro (mae) file and returns new mae file containing Epik pka prediction data for each molecule. Schrödinger "Epik" CLI-binary must be installed and path must be specified inside the script.

`03_convert_mae_to_sdf.py` 
--input: path to input file (mae.gz, mae)
--output: path to output file (sdf.gz, sdf)

--takes Schröndiger maestro (mae) file (can be gzipped) and converts it to sdf file. Schrödinger "convert" CLI-binary must be installed and path must be specified inside the script.


`04_0_filter_testmols.py` 
--input: path to input file (mae.gz, mae)
--output: path to output file (sdf.gz, sdf)
--filter: path to output file (sdf, sdf.gz)
 

--takes sdf file of initial training molecules and sdf file of training molecules (both optionally gzipped) and returns only those initial training molecules not contained in the training molecules file as sdf file. 

`04_1_split_epik_output.py` 
--input: path to input file (sdf.gz, sdf)
--output: path to output file (pkl)

--takes sdf file with molecules containing Epik pka predictions in their properties and outputs a new sdf where those molecules containing more than one pka get duplicated so that every molecules only contains one pka value. The molecule associated with each pka is the protonated form of the respective pka reaction

`04_2_prepare_rest.py` 
--input: path to input file (sdf.gz, sdf)
--output: path to output file (pkl)

--takes sdf of molecule set containing pka data and returns it as a pkl file.

`05_data_preprocess.py` 
--input: path to input file (pkl)
--output: path to output file (pkl)

--takes pkl file of molecules containing pka data and returns pytorch geometric graph data containing protonated and deprotonated graphs for every pka

`06_training.py` 
--input: set of training molecules as pyg graphs (pkl)
--output: path to output file (pkl)
--model_name: name of model architecture used (string)

Optional parameters:
--model: path for saving model or containing model for retraining (pkl)
--val: set of validation molecules as pyg graphs (pkl)
--epochs: set number of training epochs (default == 1000)
--reg: regularization dataset
-r: flag for retraining model at path given by --model

--takes training set as pkl file and trains new model or retrains existing one. 