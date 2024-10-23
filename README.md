# Basic NGS pipelines (for self-use) 

Description: 
I worked on this as a private github repo. The major tutorial is listed as the README.md file. In the final draft, all scripts will be stored in $ROOT/script/ and example data set will be stored in $ROOT/example. To work through the tutorial, I will design all script to take arguments, so user could simphy specify the following format to run the example codes and any of their real data: 

```
# Example codes run in the terminal 
./script/assembly_Nanodrop.py --input INPUT_DATA --output OUTPUT_DATA --others OTHER_ARGUMENTS
```
I designed this tutorial assuming the audience has no background knowledge in linux-based system and is new to the subject. The basic workflow is easy. If the user git clone the git repository and run the designed script in a Linux-based system following the steps listed on README.md file. This will automatically analyze any sequencing data the user have. This tutorial will have pre-written scripts which allow the users to simply run the scripts with the input and output folders specifoed to analyze their own data. 

To-do list for me to finish this tutorial 
1) I will use my real data for the demo now but it would be better to simulate some sequences (different types including single-copy genes, WGS, long-reads) and include it as examples/data
2) R1 and R2 might have multiple raw data files. Include a script if this is the case for the users.
3) Need to add in the left details. Make sure I document the steps for software installation since that might be the hardest part for the user. 

# Description of this repo: 
1) This repo decouments the basic pipeline for NGS analysis from quality control to downstream analysis including variant calling and annotation. 
2) Address different pipelines by using different types of data (NanoDrop long reads, Target-Capture, Illumina whole genome assembly, etc) 
3) Potentially have different folders with README files documenting the steps through Quality Control, Assembly, Alignment, and Variant Calling. This steps will be the four major sections in my final draft. 

# Basic workflow for NGS analysis
The basic workflow for NGS analysis including: 1. Quality Control; 2. Assembly; 3. Alignment; 4. Any downstream analysis. For this repo, I will specifically address variant calling as the major downstream analysis. This tutorial will be splitted into those four parts. 

# Directory set-up and software installation: 
First, log into any Linux-based server or using local terminal, and create a directory called ,NGS_tutorials, all operations will be conducted within this ROOT folder. 

## Method 1: Use Conda 
The easiest way to install all softwares is using Conda with instruction here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html depending on the operating system: 

```


```

# Method 2: Use Souce codes 

If conda doesn't work for any reasons, softwares could always be installed through source code. To explain the basic workflow, create a executables folder first and later we will move all binary executables we need in this folder. Then, download all essemtial softwares and move the binary executables into the executables file. 

To run any softwares, several dependencies should be installed and it depends on the operating system of the users:

On my MacBook, to download and install dependencies, download homebrew from here: https://brew.sh/ first and then run the below codes on terminals: 
```
gcc --version # Check the gcc version. Need to be version 4.8 or later. If not, check gcc website and download the latest version  
brew install isa-l
brew install libdeflate
``` 
On Linux system, run the codes below, run the below code first in the root folder that the users want to store the dependencies. The below script should be revised to the directory the user wants to install the dependencies. 

```
./script/dependency_installation.sh
```

Then, the installation script using source code (script/software_installation.sh) could be run in the $ROOT directory to install all softwares: 
```
./script/software_installation.sh
``` 

On Oct 22 2024, the following softwares were downloaded. If the user want to specify the latest version, this should be modified in the source code of $ROOT/script/software_installation.sh. 
fastQC version 0.11.9. This version is not the latest but it's pre-compiled so it is easier for installation 
fastp version 0.23.4. Latest version on Oct 22. 2024. 

```

```

# Quality control 
There are several tools for quality control, including FastQC & MultiQC, Fastp, Trimommatics. Here, I present a pipeline with 4 illumina sequencing samples. I will trim low-quality reads using fastp and checked the quality using FastQC and aggregate the results using MultiQC. 


```
```

# Assembly 

# Alignment 

# Variant Calling



