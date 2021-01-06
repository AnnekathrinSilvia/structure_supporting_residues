# Identification of residues supporting protein structure in large protein sequence alignments
# Introduction
The aim of the framework is the determination of residues supporting the protein structure in large protein sequence alignments. In this connection, protein structure supporting residues are identified as conserved and correlated protein positions which are crucial for the protein structure. After the determination the residues are used in the process of protein design.

# Requirements
The framework has some technical and data requirements, which are listed below: 
* a recent C\++ compiler that can compile C\++ 14 standard
* CMake 2.8 or later
* Python 3
* Snakemake
* a dowloaded and executible CCMPred suite
* following Python packages: Argparse, Os, PrettyTable (0.7.2 or later), Numpy, Subprocess, Shutil, Time, Biopython,

Furthermore, the framework can only deal with a MSA of protein sequences in fasta format and a quadratic substitution matrix.

# Installation
0. Install the necessary requirements.
1. Clone the repository with `git clone --recursive https://wibi-git.helmholtz-hzi.de/asl20/structure_supporting_residues.git`.
2. Change into the folder with the cloned repository (`cd structure_supporting_residues`).
3. Change into the Pipeline folder (`cd Pipeline`).
4. Clone the CCMPred suite with `git clone --recursive https://github.com/soedinglab/CCMpred.git`.
5. Change into the CCMpred folder (`cd CCmpred`).
6. Compile CCmpred with `cmake CMakeLists.txt -DWITH_CUDA=OFF`.
7. Build CCMpred with `make`.
8. Change back into the previous folder (`cd ..`).
9. Compile the C\++ code for the framework with `cmake CMakeLists.txt`.
10. Build the C\++ code with `make`.

# Setup and Usage
Before the framework can be used, some adjustments to the input paths of the individual files must be applied. 

The framework takes the necessary information of the file locations from a config file. The config file can be found in the **Input** folder. In this folder all necessary input files should be placed. These files are: 

1. A file which contains a multiple sequence alignment of the sequences which should be analyzed with the framework.
2. A normalized substitution matrix (if no normalized substitution matrix is available see below **Normalize substitution matrix**)
3. A PDB file which can be used as reference protein structure (see **PDB reference** below for a more detailed explanation).
4. A file with the unaligned sequences of the multiple sequences alignment file.
5. A header of a sequences which should be used for sequences design in the last part of the pipeline (see **Reference sequence for design** below for a more detailed explanation). 

Before the framework can be used with the input placed in the **Input** folder, the config file must be adapted. In the config file, key words can be found which address the above mentioned files (MSA, SUBSTITUTION\_MATRIX, PDB, SEQUENCES and HEADER). For these key words, the placeholder behind the equal sign must be substituted by the path to the correct file to be used for each key word respectively. 

For example, if a multiple sequence alignment file called "MSA_sequences", which is a text file, should be used as MSA file, the line in the config file must be changed to the following: 

**MSA=../Input/MSA_sequences.txt**

Please note that the correct file type must be added to the line in order for the framework to find the file. Also for a correct execution of the framework, it is necessary that all input is placed in the **Input** folder and that the preceding **../Input/** in the config file is also necessary for the framework to find the input files. The only exception to this rule is the key word **HEADER** which only needs the name of the header after the equal sign.

Besides the above mentioned key words, the config file also includes the following key words:

1. **STRICTLY\_CONSERVED**: A threshold to characterize strictly conserved residues (see below **Thresholds**).
2. **CONSERVED**: A threshold to characterize conserved residues (see below **Thresholds**).
3. **CORRELATED**: A threshold to characterize correlated residues (see below **Thresholds**).

As the config file is provided with usable thresholds, it is not mandatory to change these values. However, the user can choose to change these thresholds according to the used data set and subtitution matrix. 

After the config file has been changed and saved, the framework can be used with the following command: 

* `snakemake --cores 1 evaluateSequences`

This command executes the whole framework. The results of the framework can be found in the **Results** folder. At the end of the framework, 6 files can be found in this folder.

1. **ConservationScores.txt**: This file contains the calculated conservation scores for the sequences in the input MSA.
2. **Contacts.txt**: This file contains the predicted contacts for the sequences in the input MSA.
3. **Mapping.txt**: This file contains a mapping of conservation cores and contact numbers onto residues of a reference protein. This reference protein is taken from the PDB file in the **Input folder**.
4. **generatedSequences.txt**: The newly generated sequences.
5. **failedSequences.txt**: Those sequences, which could not be scored with ComProDes.
6. **evaluatedSequences.txt**: The newly generated sequences sorted after their likelihood to be natural protein sequences.


# Usage of only parts of the framework
As it is not alsways necessary to run the complete framework, it is possible to run only parts of the framework. 

For example, if only a new set of sequences with a different reference sequence should be generated, only the last three rules of the pipeline can be executed. For this a user must ensure that the input files of all rules, which should be executed (in this example *rule\_generateNewSequences*, *rule\_comProDes* and *rule\_evaluateSequences*) exist and are placed in the correct folders. If now only the three mentioned rules above should be executed, the following command is needed: 

* `snakemake --cores 1 evaluateSequences`

With this command Snakemake tries to execute the rule evaluateSequences. If the input necessary for the execution of this rule is not available, Snakemake will search for the rule which will output the necessary input and execute this rule first. This scheme is recursively applied until Snakemake finds an executable rule. 

With this scheme, any subset of the framework can be executed.

#Normalized Substitution Matrix
The framework calculates the conservation scores for a given input multiple sequences alignment. The used conservation scoring scheme has a bound output space in [0, 1]. To ensure the boundness of this output space, each used substitution matrix must be normalized before used for the calculation of the conservation scores. 

This repository provides a script to normalize a given quadratic substitution matrix. This script can be found in the **Pipeline** folder executed with the following command:

* `python normalizeMatrx.py "path_to_matrix_file"`

Afterwards the normalized substitution matrix must be placed in the **Input** folder (see **Setup and Usage**)

# PDB reference
The presented framework can be used to determine structure supporting residues in a given multiple sequence alignment of protein sequences. In order to classify the structure supporting residues, conservation scores and contact numbers are calculated. However, to identify the important residues, a classification scheme has to be applied. 

To facilitate the development of such a classification scheme, the framework includes a mapping step in which the conservation scores and contact numbers are mapped onto secondary structure elements of a reference protein. This mapping facilitates the interpretation of the conservation scores and contact numbers.

To ensure reliable and interpretable mapping results, the chosen reference protein should ideally be a protein with resolved structure whose sequence is included in the input MSA. If there is no reference protein with a resolved structure included in the input MSA, the used reference protein should ideally be a protein which is homologous to the majority of sequences in the input MSA. 

However, it should be noted that with increasing evolutionary distance between the reference protein sequence and the sequences in the input MSA, the reliability of the additional information which is obtained from the reference protein decreases.

# Reference sequence for design

The presented framework uses protein structure supporting residues to design new protein sequences.

For this design of new protein sequences, a reference protein sequence is needed. This reference sequence is needed to ensure that the newly designed sequences have lengths and amino acid distributions similar to the protein sequences in the input MSA. 

During the design of new protein sequences, the references sequence is used to fill the parts of the new sequences, where the MSA contains in the majority of sequences gaps (unalignable regions). These regions are substituted by the corresponding amino acids of the reference protein. 

For this choice of reference protein sequence, any sequence of the given input MSA can be used. The only restriction for this reference sequence is that the header of the sequences in the config file (see **Setup and Usage**) must be identical to the header in the given MSA and has to be unique. If the header indicated in the config file is not equal to any header in the MSA, the first sequence is used.

# Thresholds

The config file provided with the framework provides the possibility to define thresholds which characterize conserved and correlated residues. These residues are in subsequent steps of the framework used to design new protein sequences. Hence, it is necessary to characterize the protein structure supporting residues. 

In order to characterize a conserved residue, a two threshold system is applied. The structure supporting residues which are conserved are divided into two types: strictly conserved residues and conserved residues. Strictly conserved residues are such residues which are conserved in the majority of sequences in the given input MSA. Conserved residues are less conserved than strictly conserved residues. However, they still show a considerable conservation. The config file provides examples for these two thresholds with 85% conservation for a strictly conserved residue and 50% for a conserved residue.

In the same manner the correlated residues which support the protein structure have to be distinguished from those residue which are not important. For correlated residues a single threshold is applied. This threshold characterizes residues in structure supporting correlated residues and those which do not support the protein structure. The example for such a threshold provided by the config file is 50% likelihood of the correlation of two residues.

Although the config file provides example thresholds, a user is recommended to evaluate the thresholds according to the used data and substitution matrix and adapt the thresholds to the data.