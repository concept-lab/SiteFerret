# SiteFerret
An ad-hoc hierarchical clustering algorithm able to extract and rank pockets. Extraction is based on geometrical primitives  generated by a customized version of the [NanoShaper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059744) software (see below, how to use the provided executable). The ranking is based on the [Isolation Forest](https://ieeexplore.ieee.org/document/4781136) anomaly detector.

## Details:
### Clustering:
SiteFerret generates putative pockets via a tailored hierarchical clustering procedure grouping probe-spheres extracted from several calls to the NanoShaper software. NanoShaper is called externally by the script. The probes are related to the geometrical primitives representing the reentrant (concave) regions of the SES molecular surface of the protein. The clustering process is detailed in *paper in preparation*.


**Free parameters**

DEFAULT: gamma =0, beta=0.9, rp_max=3 (Angstroms)

**gamma**: How "easy" is to cluster among them probes of the same radius (larger--> larger pockets, better for large shallow sites)

**beta**: How "easy" is to cluster among them probes of different radius (larger--> deeper ramified pockets, more putative pockets generated)

**rp_max**: Terminal probe radius. Starting from 1.4 (water molecule), successive calls to NanoShaper are realized at probe radius r = 1.4,..,rp_max with increments of 0.1.
### Ranking:
Ranking is based on Isolation Forest (IF) anomaly detector. IF is provided as a scikit-learn object previously trained and loaded from a  provided binary file (in siteFerret/trainedModels)

## Requirements:
**NOTE**: The program has been tested only on a linux machine (Ubuntu 20.04)
- python3 installed
 - install patchelf:
  - sudo apt get install patchelf (ubuntu)
  - or see https://gist.github.com/ruario/80fefd174b3395d34c14 
 - The NanoShaper executable is provided but must be linked to the libraries. To do so run the install_script within *install binaries* folder and follow the prompted instructions (type *./install_script*).
 - (Reccomended) Recompile locally the shared library. This is done by running the install_script and following the instructions (gcc required).
 
 To run the install script just move into *install binaries* folder and: *./install_script* (it might be necessary to change permissions: *chmod +x install_script*)
 
 ### Download trained model:
 **Using git lfs** (recomended)

 0. git clone the folder
 1. install git lfs
 2. run: git lfs pull 

**Without using git lfs**:
download from: https://istitutoitalianotecnologia-my.sharepoint.com/:f:/g/personal/luca_gagliardi_iit_it/ErrEE6yVBGpIt_f1z43nKxkBon5Rsd_OzadlasiGV-Xh3A?e=CecSIu

Contact me if the link expired (git lfs instead should always work)

Finally <ins>copy the content in siteFerret/trainedModels/</ins>

 ### Python modules:
- numpy
- scikit-learn

## Installation
First check **Requirements**

run within the folder *pip3 install .*

**CAREFUL**: In a virtual environment you might force pip to install the package in the same directory (default behavior is to copy to another location) to not miss correct pointing to libraries. If the option -e is given (develop mode) it should prevent this problem to happen.


Then the library should be available for import (see advanced use) or use it as an executable (recomended)
## Instructions:

### Standard use
**python3 -m siteFerret \<file.pqr\>**

### Change default clustering parameters
With the *config.txt* --> see example in the confFiles folder

**OUTPUTS**:

Note: the numbering reflects the ranking.

- logfile --> contains recap of info (same printed on stdout)
- output_\<pqr_name\>.txt --> summary of ranked pockets (scores, subpockets etc..)
- errorLog.txt --> errors and warnings
- Folder 6gj6_Pfiles: contains:
 - clusterPocket\<pocket_number\>.pqr --> dummy "atoms" to represent the probe spheres. Compatible with VMD.
 - p\<pocket_number\>.off --> the triangulation of the above. Compatible with VMD.
 - p\<pocket_number\>_atm.pqr --> the protein surface atoms belonging to the pocket envelope (Recomended for practical use). Compatible with VMD.
 - Similarly for sub<number> when subpockets are available.
 - infoPocket\<pocket_number>\.txt --> info on residues and (pseudo) mouths with relative normals.
 - \<structure_name\>.vert and .face for nice triangulation in VMD of the structure. This is a "classical" NanoShaper output.

### Train and test your own Isolation Forests on a custom dataset
 The following scripts are contained in the *scripts/* folder.
 
 <ins>**NOTE**</ins>: To create a database with structures in PQR, extracted ligands, and correct corresponding mapping file for SiteFerret, use the script provided in https://github.com/concept-lab/MOAD_ligandFinder .
This script is based on the online [bindingMOAD](http://bindingmoad.org/) database to establish which ligand(s) on a structure is(are) valid.
 
**Procedure**: 
 
 1. Generate the training data using the **getData.py**This scripts produces a *features_data.pkl* binary file containing all matching pockets and sub-pockets of the algorithm run over several clustering parameters and a given list of structure-ligands pairs. In addition, the script also produces a number of information whose main one is the *statANALYSIS.txt* file containing statistical information over the generated hitting pockets (as well as the total number of putative pockets ans sub-pockets produced). The script needs an **input.prm** file detailing the clustering parameters to be explored (an example is given in the confFiles folder). Furthermore, a structure folder containing all necessary files must be set-up. The folder contains PQRs of analysed structures, ligands in xyz format, and a txt file containing the map between structures and corresponding ligands (see above regarding the database creation and its correct formatting). An example of structure folder and the type of info it should contain is given in the *script/example* folder.
 
2. Train your custom Isolation Forests using the notebook in *IF_builder*. The notebook should contain enough details. Prior to running the notebook make sure the *features_data.pkl* file is displaced in the working directory.

3. Test the goodness of your model on the same or a different dataset containing the same info described in 1. and the above NOTE. This is done by running the **test.py** script. Again the **input.prm** file must be given (but you are free to test whatever clustering parameter even if different from the ones used to generate the training set). Prior to running, check the code parameters defined at the beginning of the test.py source code (MAIN: evaluation metrics--see paper, and name and location of the trained models) which will load the 4 forests that must be present in the working directory or in the specified PATH (to be changed within the source code) and produces statistics on the goodness of the ranking. The main output file produced is statTEST_\<model_name\>.txt 
 
If a new trained Isolation Forest want to be used by default by SiteFerret,
<ins>copy the forests in siteFerret/trainedModels/</ins>. The default name of the 4 isolation forests (see **isoForest.ipynb** Jupyter notebook) are IF_geometryL, IF_geometryS, IF_chemistryL, IF_chemistryS
 
### Run SiteFerret over many structures at once
The script **loop.py** in the *scripts/* folder can be launched to produce an *output* folder containing for each structure the top10 ranked putative pockets in their *atm* format (*p\<pocket_number\>_atm.pqr*). The script by default runs on all PQR structure files present in a folder names *structures*. If the user wants to use a different folder, this can be specified at calling: **loop.py \<folder_path\>**.

This scripts, as for standard single structure call to SiteFerret, can overwrite default clustering parameters via the **config.txt** file.
 
 ## Cite
 
 TODO


