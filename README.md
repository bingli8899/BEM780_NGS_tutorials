# BEM 780 Class Project: Basic NGS pipelines for variant calling 
Variant calling is the process of identifying genetic variations, such as single nucleotide polymorphisms (SNPs) or insertions/deletions (indels), in an individual's DNA by comparing sequencing data to a reference genome. This tutorial/repo present a basic pipeline of variant calling from DNA raw reads to visualization of the variant calling. 

# Basic overview: 
What are the goals of the tutorial? 
1) This repo decouments a very basic pipeline for NGS analysis from quality control to variant calling. 
2) This pipeline assumes that the users have no experience in NGS analysis but have some knowledge in Linux/Bash and Git. 
3) The tutorial presents the basic workflow for SNP calling including: 1. Quality Control; 2. Assembly and Alignment; 3. SNP calling. 4. Downstream analysis

How to use this tutorial? 
The tutorial is designed as a github repository, where all scripts will be stored in $ROOT/script/ and example data set will be stored in $ROOT/example. $ROOT is the top directory of this repo. *It is noted that all codes are designed to be run in the $ROOT directory.* If git is not pre-installed, visit this website: https://github.com/git-guides/install-git 

To start this tutorial, run the below code in Linux/Terminal: 
```
git clone https://github.com/bingli8899/BEM780_NGS_tutorials.git
cd - # move to the GitHub repo, which I will call it $ROOT directory 
ls ./script/ # This will show scripts in the script folder
ls ./example/ # This will show all example data in the example folder
```

What's the basic workflow? 
The tutorial contains several parts: 
    1. Software installation: This step walked through how to install the softwares. 
    2. Description of the example data: This part presents two datasest (simulated and real human datasets) to be used. 
    3. Quality control
    4. Assembly and Alignment 
    5. Visualize variant calling results 
    6. References 

# Software installation: 

There are two methods to install softwares:

### Method 1: Use Conda 
The easiest way to install all softwares is to use Conda (Instruction here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). The below code creates a conda environment to work on and requires to re-activate the conda enviroments each time to resume the work: 

```
conda create -n "snp_tutorial"
conda activate snp_tutorial 
```

Then, run the following codes to install necessary softwares. In real reseach, it would be better to keep a good record of how to download the softwares and the software versions. Here, I mainly used the software samtools (https://github.com/samtools/samtools; Danecek et al. 2021), bowtie2 (https://github.com/BenLangmead/bowtie2; Langmead et al. 2012),  bcftools (https://github.com/samtools/bcftools; Danecek et al. 2021) for variant calling, fastp, multiqc and fastqc for quality trimming, sra-tools for data downloading, and seqtk for format transitioning. 

```
conda install bioconda::samtools
conda install bowtie2
conda install bioconda::bcftools
conda install bioconda::fastqc 
conda install bioconda::multiqc  
conda install bioconda::fastp
conda install bioconda::sra-tools 
conda install bioconda::seqtk
```
If any of those download failed, I would suggest to trouble shoot through conda or install with source codes (see below). 

### Method 2: Use Source codes 
If conda doesn't work for any reasons, softwares could always be installed through source codes. To explain the basic workflow, create a executables folder first and later we will move all binary executables we need in this folder. Then, download all essemtial softwares and move the binary executables into the executables file. 

If conda failed and softwares need to installed through source codes, I would recommend to check dependencies first. Run the below script to check if all dependencies are installed: 
```
chmod +x script/check_command.sh
# Usage: check_command.sh <command>
script/check_command.sh gcc
script/check_command.sh cmake
``` 
I also presented part of my method to download through source dependencies in $ROOT/script/software_installation.sh. It is noted that I tested this code in a online server with Linux system. Depending on the system, the codes (./script/software_installation.sh) may need to be modified. If there is any problem of running the code, I would suggest to run separate chunks of codes as separated by the ### in either script to figure out the issues. 

# Description of the example data 
This tutorial uses two sets of data, simulated data and real human data. 

Sequencing raw reads are typically stored in fastq (fq, fq.gz). Later, such files are processed and potentially transferred to fasta files (.fa, .fa.gz). The difference between .fasta file and .fastq file is that fasta files only contain the nucleotide or protein sequences, while fastq files contain more sequencing information such as quality score. Here, because our example data were simulated without quality score. We will start with fasta files. The real human data are stored in .fastq.gz with R1 and R2 for each sample. For paried-read sequences, each sample has two files .R1 (forward sequencing) and .R2 (backward sequencing). 

#### Simulated data 
Below section is a basic description of how I got the simulated data under $ROOT/example. This simulated_data is small in size and could be run to test if the scripts is correctly implemented and the softwares are correctly downloaded. Data are stored in $ROOT/example/example_data/*.fasta. 

Those are simulated data based on 5 accessions/taxa/samples with one simulated reference.fasta. I attached the detials about how I simulated the data to the $ROOT/example/README.md. Basically, if simphy (https://github.com/adamallo/SimPhy; Mallo et al. 2016) and seqgen (https://github.com/rambaut/Seq-Gen; Rambaut et al. 1997) got installed in one executable folder, running the below codes from $ROOT will generate the data. 
```
chmod +x ./example/simulation/simulation_seqs.sh
./example/simulation/simulation_seqs.sh <executable_folder> <configuration_folder> <output_folder>
```
Here, the five simulated fasta files are ready to use in $ROOT/example/example_data. 

#### Real Huma data 
Let's work on some real human data as well. The below steps show how to download read human data from online resouces: 

Download the reference: 
```
mkdir -p example/human_data && cd example/human_data
wget ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo*
```
Download the raw reads:
``` 
prefetch SRR098401
fastq-dump --split-files --gzip SRR098401ra 
```
The above code might be slow. If that's the case, run the below code after "prefetch" to only extract 100000 bases: 
```
fastq-dump --split-files --gzip -N 1 -X 100000 SRR098401ra 
cd ../../ 
```
Now, the real human data are downloaded in $ROOT/example/human_data

# Quality control 
The first step of most NGS pipeline is quality control. There are several tools for quality control, including Fastp, Trimommatics, etc. This step trims low-quality reads, removes short reads, removes adaptor sequences, etc. After read trimming, a check procedure using tools such as FastQC and MultiQC is recommended to double check if the trimmed reads have any issues. FastQC could check the quality of individual sequencing files, while MultiQC could aggregate the quality from multiple FastQC results. 

In $ROOT/script/quality_control.sh, the pipeline of using fastp, fastqc and multiqc is integrated. Here, because the real human data is in fastq.gz files, so I will only run quality control on $ROOT/script/human_data/*.fastq.gz. 

```
# Usage: quality_control.sh <input directory> <number of threads> <output directory>
chmod +x script/quality_control.sh 
script/quality_control.sh example/human_data 3 ./human_results 
```
The above code should output "Oyay! Outputs are in ./human_results.", which means that the code is successfully run. The first argument in quality_control.sh takes the input directory which stores the *.fastq.gz or *.fastq files. The second argument takes the number of threads to be used. The third argument specifies the output directory. This above code will generate a $ROOT/human_results directory and stores all the results. In the output directory $ROOT/human_results, there are three output directories. $ROOT/human_results/cleaned contain the cleaned reads trimmed from fastp. $ROOT/human_results/fastqc and $ROOT/human_results/multiqc stored the results of quality control files. Intermediate files from fastp are stored in $ROOT/human_results. 

It is highly recommended to check the *.html results from  $ROOT/human_results/fastqc and $ROOT/human_results/multiqc. Details about how to check the results from fastqc and multiqc are attached here: https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html. 

To keep data in the correct place, let's move the cleaned data to our example data folder: 
```
mv human_results/cleaned/*.fastq.gz 
```

# Assembly and Alignment 
There are many pipelines and methods to assemble genomes, and I will present a commonly used method, which used BcfTools , SamTools and Bowtie to assemble and align the genomes. First, assembly could be classified into three ways: 1. De-novo; 2. Reference-based method, 3. A mixture of de-novo and reference based method. 

I will present the most common method, reference-based method, to call variants. The first step is to map the reads to a reference genome with Bowtie2. Bowtie2 maps the sequencing reads to a reference genome, which in general is composed of two steps: 1. create a bowtie index based on reference genome (bowtie2-build); 2. align reads to reference (bowtie2 -x). Bowtie2 generates alignment file (*.sam)containing the information about where to map the reads to a reference. Details of the code are in the comments of $ROOT/script/alignment_variant_calling.sh. 

# Variant Calling
For NGS pipeline, converting between file formats is an essential step, and SamTools is useful to transfering between file formats and conduct a series of bioinformatic applications. To conduct variant calling, SamTools will be used to convert the output from Bowtie2 (*.sam) to *.bam files, which will be input to BcfTools for variant calling. 

Below I wrote the alignment and variant calling into one major scipts $ROOT/script/alignment_variant_calling.sh. Running the below script can conduct alignment and variant calling at the same time. The script contains comment to explain each step. 

Usage of the major script: 
```
alignment_variant_calling.sh <input_fasta_directory> <input_reference_file> <number of threads> <path to the output directory>
```
Run the major script with the simulated data first. In the below code, the first argument is the directory of the input fasta files, and the second argument is the path and file name of the reference fasta file. The third argument is the number of threads to run the script, which could be checked with "nproc", and the forth argument is the directory to the output file. 

```
chmod +x ./script/alignment_variant_calling.sh
./script/alignment_variant_calling.sh example/simulated_data example/simulated_data/reference.fasta 3 ./results_simulated
```
The above code generates a output directort in $ROOT or the current directory, which contain the bowtie2 index files (*.bt2.tmp), alignment output (simulated_data.bam), sorted BAM files (simulated_data_sorted.bam), BAM index file (simulated_data_sorted.bam.bai). The sorted bam files could be used to visualize the results from variant calling. 

If no issue is detected, let's run the script with the real human data: 
```
./script/alignment_variant_calling.sh example/human_data/ example/human_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa 4 ./results_human_variant_calling
```
The above code will output the results into ./results_human_variant_calling. 

# Visualize variant calling results: 
After variant calling, it would be good to visualize the results first. There are many ways to visualize the results from SNP calling. Let's use a text-based tool (samtool tview) to visualize the results. 

```
samtools tview ./results_human_variant_calling/human_data_sorted.bam ./example/human_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```
The above code will give a text view to view the actual results from variant calling. Type "g" could jump to anywher in the alignment to view the results. 


# Reference:
Andrew Rambaut, Nicholas C. Grass, Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees, Bioinformatics, Volume 13, Issue 3, June 1997, Pages 235–238, https://doi.org/10.1093/bioinformatics/13.3.235

Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li
GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

Diego Mallo, Leonardo De Oliveira Martins, David Posada, SimPhy : Phylogenomic Simulation of Gene, Locus, and Species Trees , Systematic Biology, Volume 65, Issue 2, March 2016, Pages 334–344, https://doi.org/10.1093/sysbio/syv082

Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012). https://doi.org/10.1038/nmeth.1923





