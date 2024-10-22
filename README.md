# Basic NGS pipelines (for self-use) 

Notes: 
I worked on this as a private github repo. Now, I will list this as a draft to BEM780 tutorial assignment. I designed this tutorial assuming the audience has no background knowledge in linux-based system and is totally new to the subject. In the final draft, depending on the knowledge level of my classmates, I might omit some unnecessary details. 

The basic workflow is easy. If the user git clone the git repository and run the designed script in a Linux-based system following the steps listed on README.md file. This will automatically analyze any sequencing data the user have. 

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

First, log into any Linux-based server or using local terminal, and create a directory called ,NGS_tutorials, all operations will be conducted within this ROOT folder. Then, create a executables folder and install all binary executables we need in this folder. 

Then, download all essemtial softwares and move the binary executables into the executables file. 
Most of the softwares could be downloaded via conda. Since this tutorial is designed for beginner level and not all system could build conda environment, I only present the installation through source code and sumamrized all installation steps through one script $ROOT/script/software_installation.sh. User could simply run ./script/software_installation.sh to get everything downloaded. 

To run those softwares, several dependencies should be installed and it depends on the operating system of the users:
On my MacBook, to download and install dependencies, download homebrew from here: https://brew.sh/ first and then run the below codes on terminals: 

```
gcc --version # Check the gcc version. Need to be version 4.8 or later 
brew install isa-l
brew install libdeflate
``` 
On Linux system, run the codes below, run the code script/dependency_installation.sh first on a folder that the users want to store the dependencies. 

On Oct 22 2024, the following softwares were downloaded. If the user want to specify the latest version, this should be modified in the source code of $ROOT/script/software_installation.sh. 
fastQC version 0.11.9. This version is not the latest but it's pre-compiled so it is easier for installation 
fastp version 0.23.4. Latest version on Oct 22. 2024. 
```

```

# Quality control 
There are several tools for quality control, including FastQC & MultiQC, Fastp, Trimommatics. Here, I present a pipeline with 4 illumina sequencing samples. I will trim low-quality reads using fastp and checked the quality using FastQC and aggregate the results using MultiQC. 

# Assembly 

# Alignment 

# Variant Calling



