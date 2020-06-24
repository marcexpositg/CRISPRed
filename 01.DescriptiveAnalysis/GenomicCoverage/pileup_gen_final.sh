#!/bin/bash
#SBATCH -p normal 
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu 1000 
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marc.exposit@upf.edu


#SBATCH -e stderr_filt_%j.err
#SBATCH -o stdout_filt_%j.out

### LOAD INITIAL MODULES ###
module load SAMtools/1.9-foss-2016b


### SAVE VARIABLES ###
GenomicReads=$1
OutFile=$2

### PIPELINE ###

samtools mpileup  -f ../Reference_Genomes/GCA_001632575.1_C3H_HeJ_v1.fa -l ./C3H_gRNA-coordinates.bed $GenomicReads -o $OutFile -a

#### Command:
#### execute inside /Coverage/ folder
#### sbatch pileup_prova.sh <genomicReads.bam> <outputFile.pileup>
