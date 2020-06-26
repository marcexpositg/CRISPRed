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
module load Python/3.6.6-foss-2018b

### SAVE VARIABLES ###
plasmid_seq_file=$1
gRNA_library_file=$2
output_file_name=$3

### PIPELINE ###

# Count distribution of sequenced gRNA library
python gRNA_dist_count.py ${plasmid_seq_file} ${gRNA_library_file} ${output_file_name}


#### Command:
#### sbatch gRNA_dist_count.sh <plasmid_seq_file> <gRNA_library_file> <output_file_name.csv>