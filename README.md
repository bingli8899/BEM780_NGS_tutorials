# BEM 780 Class Project: Basic NGS pipelines for SNP calling 

# Goals of this repo: 
1) This repo decouments a very basic pipeline for NGS analysis from quality control to downstream analysis including variant calling and SNP detection. 
2) This pipeline should be designed as well-established and could be applicable to wide range of NGS data with modification. 
3) This pipeline assumes that the users have no experience in NGS analysis but have some knowledge in Git and Terminal. 

# Before we get started: 
The tutorial is designed as a github repository, where all scripts will be stored in $ROOT/script/ and example data set will be stored in $ROOT/example. To find this code, run the below code in terminal. If git is not pre-installed, visit this website: https://github.com/git-guides/install-git 
```
git clone https://github.com/bingli8899/BEM780_NGS_tutorials.git
```
## Description of the example data 
Below section is a basic description of how I got the example data. For this tutorial, I will use example/example_data/*.fastq.gz . Those are simulated data based on 6 accessions/taxa/samples. I attached the detials about how I simulated the data to the $ROOT/example/README.md (). Basically if simphy and seqgen got installed in one executable folder, running the below codes from $ROOT will generate the data: 
```
chmod +x ./example/simulation/simulation_seqs.sh
./example/simulation/simulation_seqs.sh <executable_folder> <configuration_folder> <output_folder>
```
Here, I attached the data in the $ROOT/example/example_data already to avoid the processes of installing simphy and seqgen and conducting data format transformation. The six fastq.gz files are ready to use. 

# Example codes run in the terminal 
./script/assembly_Nanodrop.py --input INPUT_DATA --output OUTPUT_DATA --others OTHER_ARGUMENTS

I will design this tutorial assuming the audience has no background knowledge in linux or perl and the user is new to NGS analysis. The basic workflow is easy. If the user git clone the git repository and run the designed script in a Linux-based system following the steps listed on README.md file. This will automatically analyze any sequencing data the user have. This tutorial will have pre-written scripts which allow the users to simply run the scripts with the input, output folders, and parameters, which allow the user to adjust the settings based on their own data while allowing for a easy workflow. 

Running this tutorial will need to download github. Then, user could git clone the repo and the example and script folders are ready to use: 
```
git clone https://github.com/bingli8899/NGS_pipeline_tutorials.git # Now it is private. 
cd - # move to the GitHub repo, which I will call it $ROOT directory 
ls ./script/ # This will show scripts in the script folder
ls ./example/ # This will show all example data in the example folder
```
All codes are designed to be run in the $ROOT directory. 

The structure of this Gitub repo so far: 

```
$ROOT 
---script: 
    dependency_installation.sh
    software_installation.sh
    simulation.jl # script to simulation DNA sequences from Seq-Gen -- not useful to users but keep records
    concatenate.sh
---example:
    simulationed_dna:
        A_R1.fastq, A_R2.fastq
        B_R1.fastq, B_R2.fastq
        C_R1.fastq, C_R1.fastq
        D_R1.fastq, D_R2.fastq
    grand_truth.tre # ((A,B),C),D) 
--executable: 
    Any binary executable that could not be installed via Conda

```

# Basic workflow for NGS analysis
The basic workflow for NGS analysis including: 1. Quality Control; 2. Assembly; 3. Alignment; 4. Any downstream analysis. For this repo, I will specifically address variant calling (SNP calling) as the major downstream analysis. This main tutorial will be splitted into those four parts. Before that, let's download all required softwares and set up directories. 

# Directory set-up and software installation: 
First, log into any Linux-based server or using local terminal, and create a directory called ,NGS_tutorials, all operations will be conducted within this ROOT folder. 

## Method 1: Use Conda 
The easiest way to install all softwares is using Conda with instruction here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html depending on the operating system: 

Then, run the following codes to create a conda enviroment to run scripts: 

```
conda create -n "ngs_tutorial"
conda activate ngs_tutoria
```

Then, run the following codes to install necessary softwares. In real reseach, it would be better to keep a good record of how to download the softwares and the software versions:  
```
conda install bioconda::fastqc
```
Note: I am still thinking of what pieplines to include here. There are so many softwares to do the same stuff. This software installation procedure should be modified later. 

## Method 2: Use Souce codes 

If conda doesn't work for any reasons, softwares could always be installed through source code. To explain the basic workflow, create a executables folder first and later we will move all binary executables we need in this folder. Then, download all essemtial softwares and move the binary executables into the executables file. 

To run any softwares, several dependencies should be installed and it depends on the operating system of the users:

On my MacBook, to download and install dependencies, download homebrew from here: https://brew.sh/ first and then run the below codes on terminals: 

```
gcc --version # Check the gcc version. Need to be version 4.8 or later. If not, check gcc website and download the latest version  
brew install isa-l
brew install libdeflate
``` 
Note: Since MAC could easily install conda (method 1), so I will change the above methods to focus on Linux system without using homebrew. This requires me to change the software_installation.sh later.  

On Linux system, run the codes below, run the below code first in the root folder that the users want to store the dependencies. The below script should be revised to the directory the user wants to install the dependencies. 

```
./script/dependency_installation.sh
```

Then, the installation script using source code (script/software_installation.sh) could be run in the $ROOT directory to install all softwares: 

```
./script/software_installation.sh
``` 
If there is any problem of running the code, I would suggest to run separate chunks of codes as separated by the #### in either script to figure out the issues. 

# Understand the dataset
Sequencing raw reads are typically stored in fastq (fq, fq.gz) or sometimes fasta (.fa, .fa.gz) files. Most early stage NGS analysis softwares could recognize the zipped (.gz) format. For paried-read sequences, each sample has two files .R1 (forward sequencing) and .R2 (backward sequencing). 

Here, I present a pipeline with 4 illumina sequencing samples (Reminder: Need to simulate the seuqence first).  

Sometimes, each .R1 and .R2 might have multiple files such as .R1_1, .R1_2, or .R2_1, R2)2 etc. In the examplem data, sample B has two files for each paired-end read (Reminder: need to simulate the data later). If this is the case, run the below procedure to concatenate different files: 

```
./script/concatenate.sh --input INPUT_FILE --output OUTPUT_DIR
```

# Quality control 
After downloading necessary softwares, there are several tools for quality control, including FastQC & MultiQC, Fastp, Trimommatics. This step trim low-quality reads, remove short reads, remove adaptor sequences, etc. After read trimming, a check procedure is recommended to double check if the trimmed reads have any issues. 

I will trim low-quality reads using fastp and checked the quality using FastQC and aggregate the results using MultiQC. (Reminder: This should be written is three bash script which takes arguments from the softwares but can process multiple samples through loops. The reason to design those scripts is to allow user with no experience in bash loop to multi-process bioinformatic pipeline in parallel).

# Assembly and Alignment 
There are two many pipelines and methods to assemble genomes, and I will present the most commonly used method, which used BCFtools, BWA, SamTools and Bowtie to assemble and align the genomes. 

First, assembly could be classified into three ways: 1) De-novo; 2) Reference-based method, 3) A mixture of de-novo and reference based method. 

De-novo requires specific softwares designed for different types of dataset, so I will present a reference-based assembly methods here by mapping reads to the reference genome via bwa. 

Note to myself: Need to write the pipeline which could be universal to most dataset. Check and modify my previous bwa pieplines. 

# Variant Calling
Here, I will use BCFtools to conduct variant calling and freebayes for SNP detection. 


# Reference:


