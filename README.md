# Basic NGS pipelines (for self-use) 

# Goals of this repo: 
1) This repo decouments a very basic pipeline for NGS analysis from quality control to downstream analysis including variant calling and SNP detection. 
2) This pipeline is well-established and could be applicable to wide range of NGS data with modification. However, the example data here presented will use Illumina paired-read with 150 base pair per read. 
PS: I will simulate the data with a grand truth relationship ((A,B),C),D) with 100 genes. I will concatenate the genes and add into random regions in between and suffle the sequences into 150 sequences. By this, I know the actual genes and SNP and can compare the results. Alternatively, if there is not enough time for me to simulate the sequences, I will find online resources for a demo.  

# Description: 
I worked on this as a private github repo. The major tutorial is listed as the README.md file. In the final draft, all scripts will be stored in $ROOT/script/ and example data set will be stored in $ROOT/example. To work through the tutorial, user could simphy specify the directory to their input path to run the example codes and any of their real data: 

```
# Example codes run in the terminal 
./script/assembly_Nanodrop.py --input INPUT_DATA --output OUTPUT_DATA --others OTHER_ARGUMENTS
```
I will design this tutorial assuming the audience has no background knowledge in linux-based system and is new to NGS analysis. The basic workflow is easy. If the user git clone the git repository and run the designed script in a Linux-based system following the steps listed on README.md file. This will automatically analyze any sequencing data the user have. This tutorial will have pre-written scripts which allow the users to simply run the scripts with the input and output folders specifoed to analyze their own data. 

Running this tutorial will need to download github. Then, user could git clone the repo and the example and script folders are ready to use: 
```
git clone https://github.com/bingli8899/NGS_pipeline_tutorials.git
cd - # move to the GitHub repo, which I will call it $ROOT directory 
ls ./script/ # This will show scripts in the script folder
ls ./example/ # This will show all example data in the example folder
```
All codes are designed to be run in the $ROOT directory. 

# Basic workflow for NGS analysis
The basic workflow for NGS analysis including: 1. Quality Control; 2. Assembly; 3. Alignment; 4. Any downstream analysis. For this repo, I will specifically address variant calling as the major downstream analysis. This tutorial will be splitted into those four parts. 

# Directory set-up and software installation: 
First, log into any Linux-based server or using local terminal, and create a directory called ,NGS_tutorials, all operations will be conducted within this ROOT folder. 

## Method 1: Use Conda 
The easiest way to install all softwares is using Conda with instruction here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html depending on the operating system: 

Then, run the following codes to create a conda enviroment to run scripts: 
conda create -n "ngs_tutorial"
conda activate ngs_tutoria


Then, run the following codes to install necessary softwares: 
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
If there is any problem of running the code, I would suggest to run chunks of codes as separated by the #### in either script to figure out the issues. 

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



