#!/bin/bash
#SBATCH --job-name=pollock_fastqc
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

##### ENVIRONMENT SETUP ##########
DATADIR=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3
OUTDIR=/mmfs1/gscratch/scrubbed/elpetrou/pollock/fastqc/Plate3
MYSUFFIX=.fastq

#### CODE FOR JOB #####
## make output directory
mkdir $OUTDIR

## navigate into directory holding fastq files
cd $DATADIR 

# Run fastqc on all files that end in a particular suffix ( in this case, .fastq)
for SAMPLEFILE in *$MYSUFFIX
do
    fastqc -f fastq --extract \
    $SAMPLEFILE
done

## Move all the fastqc results files (ending in _fastqc, .zip, .html) to the output directory
mv *.html $OUTDIR
mv *.zip $OUTDIR
mv *_fastqc $OUTDIR
