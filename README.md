# TSR-3D representation of Proteins and Proteins Clustering
This repository contains the code files used for representing proteins 3D structure as integers called "keys" based on Triangular Spatial Relationship extended for 3D objects. This also includes the sample dataset generation and all the steps followed by key generation to proteins clustering and motif identification.


The steps are mentioned in the flow diagram below.
![](https://user-images.githubusercontent.com/16475732/100969938-9074b200-34f9-11eb-8c9f-45771abcd777.png)

Description of the folder structure:
* Classification: It contains the code files of Key Generation, Vector representation of the protein with keys as features and their frequencies as the values, Jaccard similarity calculation and visualization generation files. The visualization generates heatmaps, dendrgrams and clustermaps of the clustering results.
* Descritization: This is a precursor for key generation. It contains the code files of Adaptive Unsupervised Iterative Discretization algorithm used in generating bins and bin boundaries to use in key calculation.
* Helper Functions: It contains miscellaneous code files used in the project.
* Loni_Scripts: The key generation is computationally intensive. It also requires huge storage. Majority of our calculations were executed at LONI (Louisiana Optical Network Initiative) supercomputers. This folder contains sample shell scripts used.
* Sample Generation: This folder contains the code to download PDB files from PDB databank https://www.rcsb.org/.

* data: The data folder contains the datasets details.

## Instructions for Execution:
Installing the dependencies:
The following command will install the packages according to the configuration file requirements.txt

```
$ pip install -r requirements.txt
```

Put requirements.txt in the directory where the command will be executed. If it is in another directory, specify the path.

Execution:

Below are sample commands for performing Clustering.

Sample Generation:

```
$ python ../Sample_Generation/generate_samples.py --path <path_to_the_all_datasets_folder> --sample_name <dataset_or_sample_name>
```
Parallel Key Generation:
```
$ python2 ../Classification/lib/key_generation_parallel.py --path <path_to_the_all_datasets_folder> --sample_name <dataset_or_sample_name> --thetaBounds <theta_bin_boundaries>  --distBounds <distance_bin_boundaries> 
```
If thetaBounds and distBounds are not passed default bin boundaries of 29 bins for Theta and 35 bins for Distance are considered.

Clustering:
A single script will enable one to run steps from Key Genaration to Clustering. However key generation step using this script is executed sequentially. Hence consider running the previous parallel mode key generation script before running this. In this case omit '1' in the steps argument.

```
$ python2 ../Classification/__main__.py --path <path_to_the_all_datasets_folder> --sample_name <dataset_or_sample_name> --thetaBounds <theta_bin_boundaries>  --distBounds <distance_bin_boundaries>  --steps '2,3,4,5'
```
If thetaBounds and distBounds are not passed default bin boundaries of 29 bins for Theta and 35 bins for Distance are considered.

We used LONI servers for faster computaion and huge storage resources. If the dataset is smaller (<100 proteins or drugs) one can use their local computers. Hence above commands are writte in a shell script file for execution of all stages at once.

### Contributors

* [Sarika Kondra](mailto:sarika.vm35@gmail.com?subject=[GitHub] TSR-3D)
