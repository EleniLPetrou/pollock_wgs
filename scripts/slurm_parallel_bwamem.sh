#!/bin/bash

## Job-specific Variables
###############################################
# Specify the path to your scipt and its name
MYSCRIPT=/mmfs1/home/elpetrou/pollock_scripts/bwa_mem.sh

# Path to directory with fastq files
DATADIR=/gscratch/scrubbed/elpetrou/pollock/fastq_trimmed/fastq_PWS20

# Specify suffix of files to analyze
SUFFIX1=_R1_001.trim.fastq

################################################
# Make your script executable
chmod +x $MYSCRIPT

# Run the script on each file in the current directory
for i in $DATADIR'/'*$SUFFIX1
do
	echo $i
	sbatch $MYSCRIPT $i
	sleep 1
done
