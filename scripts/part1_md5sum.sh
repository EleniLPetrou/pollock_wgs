# Eleni 20210610
# This file documents how I copied the pollock sequencing data, checked it for file corruption, and extracted the archives.

#request compute node
srun -p compute-hugemem -A merlab --nodes=1 --ntasks-per-node=1 --time=1-02:00:00 --mem=100G --pty /bin/bash

## Copy all raw sequencing files from Ellie's directory to mine
ELLIEDIR=/gscratch/scrubbed/ekbors/raw_sequences
MYDIR=/gscratch/scrubbed/elpetrou/pollock/raw_sequences

cp $ELLIEDIR'/'* $MYDIR

##########################################################

# Navigate to directory with raw pollock sequences
cd /gscratch/scrubbed/elpetrou/pollock/raw_sequences

# check raw sequence files against md5sum files to make sure that they are not corrupted
md5sum EKB_plate3_93plex.tar
#abdefde38d688ef8d89d2a29cfb7ba31  EKB_plate3_93plex.tar
#abdefde38d688ef8d89d2a29cfb7ba31  EKB_plate3_93plex.tar # It matches! Hurray!

md5sum Ellie_Lane1_plate_1_5_7.tar
#cc5c7c332fdc354a829580af73d0e14b  Ellie_Lane1_plate_1_5_7.tar
#cc5c7c332fdc354a829580af73d0e14b  Ellie_Lane1_plate_1_5_7.tar # It matches! Hurray!

md5sum Ellie_Lane2_plates_2_4_6.tar
#b7e3710aeae971aadc4c55e7e9a0e6c8  Ellie_Lane2_plates_2_4_6.tar
#b7e3710aeae971aadc4c55e7e9a0e6c8  Ellie_Lane2_plates_2_4_6.tar # It matches! Hurray!

##############################################################
# extract the tar files
tar -xvf EKB_plate3_93plex.tar

tar -xvf Ellie_Lane1_plate_1_5_7.tar

tar -xvf Ellie_Lane2_plates_2_4_6.tar

#############################################################
# Rename and move each extracted folder into a new directory: 

#Rename the folders
mv Ellie_Bors_Plate3_done Plate3
mv Ellie_Lane1_done Lane1
mv Ellie_Lane2_done Lane2

# Move the folders into a new directory
FASTQDIR=/gscratch/scrubbed/elpetrou/pollock/fastq

mkdir $FASTQDIR

mv Plate3 $FASTQDIR
mv Lane1 $FASTQDIR
mv Lane2 $FASTQDIR


