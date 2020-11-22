# TSR-3D representation of Proteins and Proteins Clustering
This repository contains the code files used for representing proteins 3D structure as integers called "keys" based on Triangular Spatial Relationship extended for 3D objects. This also includes the sample dataset generation and all the steps followed by key generation to proteins clustering and motif identification.


The steps are mentioned in the flow diagram below.

Description of the folder structure:
* Classification: It contains the code files of Key Generation, Vector representation of the protein with keys as features and their frequencies as the values, Jaccard similarity calculation and visualization generation files.
* Descritization: This is a precursor for key generation. It contains the code files of Adaptive Unsupervised Iterative Discretization algorithm used in generating bins and bin boundaries to use in key calculation.
* Helper Functions: It contains miscellaneous code files used in the project.
* Loni_Scripts: The key generation is computationally intensive. It also requires huge storage. Majority of our calculations were executed at LONI (Louisiana Optical Network Initiative) supercomputers. This folder contains sample shell scripts used.
* Sample Generation: This folder contains the code to download PDB files from PDB databank https://www.rcsb.org/.
