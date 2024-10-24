#!/bin/bash

# The script that download and install all softwares 
# Run this script by ./script/software_installation.sh in the $ROOT folder 
# It will move all executables to a executables folder after installation 

mkdir -p executables
mkdir -p downloaded_softwares && cd downloaded_softwares

# This script assumes the user has no git installed so use wget instead 
# Download FastQC: 
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc*
rm *.zip
cd FastQC 
chmod +x fastqc
cd ../../executables/
ln -s ../downloaded_softwares/FastQC/fastqc fastqc
# Print a warning if fastqc is not installed correctly 
if ./fastqc --help > /dev/null 2>&1; then
    echo "FastQC is installed and working correctly. Stop and check for installation."
else
    echo "Warning: FastQC is not installed or not working correctly. Move forward."
fi
cd ../downloaded_softwares

# Download FastP: 
wget https://github.com/OpenGene/fastp/archive/refs/tags/v0.23.4.zip
unzip *.zip 
cd fastp*
make -j # User might want to make sure they have gcc+ installed 

# This could be changed to trimmomatic if easier 









