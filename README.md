# BEM 780 Class Project: Basic NGS pipelines for SNP calling 

# Goals of this repo: 
1) This repo decouments a very basic pipeline for NGS analysis from quality control to downstream analysis including variant calling and SNP detection. 
2) This pipeline should be designed as well-established and could be applicable to wide range of NGS data with modification. 
3) This pipeline assumes that the users have no experience in NGS analysis but have some knowledge in Git and Terminal. 

# Before we get started: 
The tutorial is designed as a github repository, where all scripts will be stored in $ROOT/script/ and example data set will be stored in $ROOT/example. To find this code, run the below code in terminal. If git is not pre-installed, visit this website: https://github.com/git-guides/install-git 
```
git clone https://github.com/bingli8899/BEM780_NGS_tutorials.git
cd - # move to the GitHub repo, which I will call it $ROOT directory 
ls ./script/ # This will show scripts in the script folder
ls ./example/ # This will show all example data in the example folder
```

All codes are designed to be run in the $ROOT directory. The structure of this Gitub repo is described: 

```
$ROOT 
---script: 
    dependency_installation.sh
    software_installation.sh
---example:
    README.md # This documents how I simulated the example data
    example_data:
        A_R1.fastq, A_R2.fastq
        B_R1.fastq, B_R2.fastq
        C_R1.fastq, C_R1.fastq
        D_R1.fastq, D_R2.fastq
```

## Description of the example data 
Below section is a basic description of how I got the example data. For this tutorial, I will use example/example_data/*.fastq.gz . Those are simulated data based on 6 accessions/taxa/samples. I attached the detials about how I simulated the data to the $ROOT/example/README.md (). Basically, if simphy (https://github.com/adamallo/SimPhy; Mallo et al. 2016) and seqgen (https://github.com/rambaut/Seq-Gen; Rambaut et al. 1997) got installed in one executable folder, running the below codes from $ROOT will generate the data: 

```
chmod +x ./example/simulation/simulation_seqs.sh
./example/simulation/simulation_seqs.sh <executable_folder> <configuration_folder> <output_folder>

```
Here, I attached the data in the $ROOT/example/example_data already to avoid the processes of installing simphy and seqgen and conducting data format transformation. The six fa.gz files are ready to use.

## Understand the dataset
Sequencing raw reads are typically stored in fastq (fq, fq.gz) or sometimes (very rarely) in fasta (.fa, .fa.gz) files. Most early stage NGS analysis softwares could recognize the zipped (.gz) format. For paried-read sequences, each sample has two files .R1 (forward sequencing) and .R2 (backward sequencing). The difference between .fasta file and .fastq file is that fasta files only contain the nucleotide or protein sequences, while fastq files contain more sequencing information such as quality score. Here, because our example data were simulated without quality score. We will start with fa.gz files. 

```

```


# Basic workflow for SNP calling 
The basic workflow for NGS analysis including: 1. Quality Control; 2. Alignment; 3. SNP calling. To run that, let's download the required softwares. 

## Directory set-up and software installation: 
First, log into any Linux-based server or using local terminal, and create a directory called ,NGS_tutorials, all operations will be conducted within this $ROOT folder. 

There are two methods to install softwares: 
### Method 1: Use Conda 
The easiest way to install all softwares is using Conda with instruction here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html depending on the operating system. Here, I intentionally chose the softwares so all required softwares could be installed through conda. 

Then, run the following codes to create a conda enviroment to run scripts. This creates a conda environment to work on and requires to re-activate the conda enviroments each time to resume the work: 

```
conda create -n "snp_tutorial"
conda activate snp_tutoria
```

Then, run the following codes to install necessary softwares. In real reseach, it would be better to keep a good record of how to download the softwares and the software versions:  
```
conda install bioconda::samtools


conda install bioconda::fastqc
conda install bioconda::fastp

```

### Method 2: Use Source codes 
If conda doesn't work for any reasons, softwares could always be installed through source code. To explain the basic workflow, create a executables folder first and later we will move all binary executables we need in this folder. Then, download all essemtial softwares and move the binary executables into the executables file. 

This method is not recommended for first-time SNP calling. I presented my method of downloading through source dependencies in $ROOT/script/dependency_installation.sh and $ROOT/script/software_installation.sh. 

```
./script/dependency_installation.sh
./script/software_installation.sh
``` 
If there is any problem of running the code, I would suggest to run separate chunks of codes as separated by the #### in either script to figure out the issues. 

# Quality control 
After downloading necessary softwares, there are several tools for quality control, including FastQC & MultiQC, Fastp, Trimommatics. This step trim low-quality reads, remove short reads, remove adaptor sequences, etc. After read trimming, a check procedure is recommended to double check if the trimmed reads have any issues. 

I will trim low-quality reads using fastp and checked the quality using FastQC and aggregate the results using MultiQC. (Reminder: This should be written is three bash script which takes arguments from the softwares but can process multiple samples through loops. The reason to design those scripts is to allow user with no experience in bash loop to multi-process bioinformatic pipeline in parallel).

# Alignment 
There are two many pipelines and methods to assemble genomes, and I will present the most commonly used method, which used BCFtools, BWA, SamTools and Bowtie to assemble and align the genomes. First, assembly could be classified into three ways: 1) De-novo; 2) Reference-based method, 3) A mixture of de-novo and reference based method. 

Here, I will present the most common and computationally cheap method, reference-based method, to call SNPs. The basic 

# Variant Calling
Here, I will use BCFtools to conduct variant calling and freebayes for SNP detection. 


# Reference:
Diego Mallo, Leonardo De Oliveira Martins, David Posada, SimPhy : Phylogenomic Simulation of Gene, Locus, and Species Trees , Systematic Biology, Volume 65, Issue 2, March 2016, Pages 334–344, https://doi.org/10.1093/sysbio/syv082
Andrew Rambaut, Nicholas C. Grass, Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees, Bioinformatics, Volume 13, Issue 3, June 1997, Pages 235–238, https://doi.org/10.1093/bioinformatics/13.3.235



