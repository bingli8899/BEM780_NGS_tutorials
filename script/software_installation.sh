#!/bin/bash

# The script that download and install all softwares 
# Run this script by ./script/software_installation.sh in the $ROOT folder 
# It will move all executables to a executables folder after installation 

mkdir -p executables
mkdir -p downloaded_softwares && cd downloaded_softwares

### Download SamTools 
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 
# This downloades the newest samtool (v.1.21) by Nov 27 2024, but this version could be updated 

tar -xf samtools*
rm *.tar* 
cd samtools*
./configure
make
cd ../../executables
ln -s ../downloaded_softwares/samtools*/samtools

./script/check_command.sh samtools

### Download bcftools 
cd ../downloaded_softwares
wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 
# Again, the version could be changed to the newest version 
tar -xf bcftools*
rm *tar*
cd bcftools*
./configure
make
cd ../../executables
ln -s ../downloaded_softwares/bcftools*/bcftools 

## Download Bowtie2 








######################### Fastp and FastQC########### 
# This script assumes the user has no git installed so use wget instead 
# Download FastQC: 
# wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
# unzip fastqc*
# rm *.zip
# cd FastQC 
# chmod +x fastqc
# cd ../../executables/
# ln -s ../downloaded_softwares/FastQC/fastqc fastqc
# Print a warning if fastqc is not installed correctly 
# if ./fastqc --help > /dev/null 2>&1; then
    echo "FastQC is installed and working correctly. Stop and check for installation."
# else
#    echo "Warning: FastQC is not installed or not working correctly. Move forward."
# fi
# cd ../downloaded_softwares

# Download FastP: 
# wget https://github.com/OpenGene/fastp/archive/refs/tags/v0.23.4.zip
# unzip *.zip 
# cd fastp*
# make -j # User might want to make sure they have gcc+ installed 

# This could be changed to trimmomatic if easier 









