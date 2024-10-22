# NGS_pipeline_tutorials
A detailed description of very basic NGS pipelines (for self-use)

Notes: 
This will be listed as my provate github repo later for my self use in the future. Now, I will list this as a draft to BEM780 tutorial assignment. I designed this tutorial assuming the audience has no background knowledge in linux-based system and is totally new to the subject. In the final draft, depending on the knowledge level of my classmates, I might omit some unnecessary details. 

# Description of this repo: 
1) This repo decouments the basic pipeline for NGS analysis from quality control to downstream analysis including variant calling and annotation. 
2) Address different pipelines by using different types of data (NanoDrop long reads, Target-Capture, Illumina whole genome assembly, etc) 
3) Potentially have different folders with README files documenting the steps through Quality Control, Assembly, Alignment, and Variant Calling. This steps will be the four major sections in my final draft. 

# Basic workflow for NGS analysis
The basic workflow for NGS analysis including: 1. Quality Control; 2. Assembly; 3. Alignment; 4. Any downstream analysis. For this repo, I will specifically address variant calling as the major downstream analysis. This tutorial will be splitted into those four parts. 

First, log into any Linux-based server or using local terminal, and create a directory called ,NGS_tutorials, all operations will be conducted within this ROOT folder. Then, create a executables folder and install all binary executables we need in this folder. 

'''
cd # to the directory to create the ROOT folder, NGS_tutorials 
mkdir -p executables && cd executables 

'''

# Quality control 
There are several tools for quality control, including FastQC & MultiQC, Fastp, Trimommatics. Here, I present a pipeline with 4 illumina sequencing samples. I will trim low-quality reads using fastp and checked the quality using FastQC and aggregate the results using MultiQC. 



