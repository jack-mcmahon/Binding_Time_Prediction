## Residence Time Prediction with Alpha-Synuclein

This project aims to explore the use of deep learning models for protein-ligand binding event residence time prediction using data from long-scale MD simulations. The data used in this project comes from the following DESRES paper: https://pubs.acs.org/doi/10.1021/jacs.1c07591

I use trajectories for Ligand 47 and Fasudil, though the project is intended as a general format that could be used with other ligands/trajectories with some data preprocessing.

### Prerequisites

This project requires the following software packages:
* numpy
* pandas
* scipy
* pytorch
* pytorch-tabnet
* sklearn
* meeko
* vina
* rdkit
* openbabel
* mdtraj
* ADFR: https://ccsb.scripps.edu/adfr/downloads/

## Usage

Basic Use: If you want the most straightforward introduction to the project, start with Run_pca_docking_comparison.ipynb to explore data preprocessing and train a TabNet model on PCA inputs. Data preprocessing assumes that molecular interactions data has previously been extracted from an MD simulation and is available in either numpy array or pandas dataframe format. The script is designed to process continuous binding events defined as the ligand being within 6A of the protein. Arguments for the dimensionality of the PCA input and TabNet construction are easily customizable and should be tuned for a specific dataset.

Docking Scores: After training a regression model, this project provides scripts to obtain the precursor PDB/MOl2/PDQBT files from the MD trajectory and input to Autodock Vina to obtain docking scores. These scripts are found in the Docking_Scores folder, as well as my processed scores for reference. After obtaining scores, Run_pca_docking_comparison.ipynb provides some functions to explore docking score distribution and correlation to residence time predictions.

Other Experiments: If you're interested in other formats for input data processing, the experiment_sweep folder contains a range of different scripts for setting up the following:
* Simple ligand-residue distance input
* Binarized Molecular interactions input
* Single/multiple frame input
* Average PCA input for a range of settings

## Acknowledgments

* Thanks to Apara Chakraborty for help with molecular interactions data processing and Anjali Dhar for help processing PDBQT files for docking with vina

