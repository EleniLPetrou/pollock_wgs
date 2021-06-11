#!/bin/bash
#SBATCH --job-name=pollock_gunzip
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-12:00:00
## Memory per node
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

# Move the folders into a new directory
PLATE3=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3
LANE1=/gscratch/scrubbed/elpetrou/pollock/fastq/Lane1
LANE2=/gscratch/scrubbed/elpetrou/pollock/fastq/Lane2

# unzip the fastq files in each directory
cd $PLATE3
gunzip *fastq.gz

cd $LANE1
gunzip *fastq.gz

cd $LANE2
gunzip *fastq.gz
