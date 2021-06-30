#!/bin/bash
#SBATCH --job-name=elp_mergelanes
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=12:00:00
## Memory per node
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ##########
DATADIR=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3 #directory with fastq files
OUTDIR=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3_merged_fastq

# Navigate to directory with fastq files and create a text file that contains of the unique individual id basenames (naming convention = POP_INDIV_stuff_LaneID_ReadID_001.fastq)
# cut -d specifies that the delimiter (in our case, and underscore)
# -f 1,2 will retain the first two fields (POP_INDIV) 

cd $DATADIR
ls *R1* | cut -d _ -f 1,2 | sort | uniq > sample_basenames.txt

# sanity check
#cat sample_basenames

# Run loop to merge fastq files
# For each sample, print content to standard output for both lane 1 and lane 2 fastq files.
# Do this separately for the forward  reads (*R1) and the reverse reads (*R2).
# Save the standard output to a new file (thus merging the fastq from Lane 1 and Lane 2 for each sample and read orientation)


for sample in `cat sample_basenames.txt`
do
	cat ${sample}*R1_001.fastq > ${sample}_merged_R1_001.fastq
	cat ${sample}*R2_001.fastq > ${sample}_merged_R2_001.fastq
done

# Move merged files to new directory
mv *merged*.fastq $OUTDIR
