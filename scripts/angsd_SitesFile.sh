# How to make a sites file for angsd from a .mafs file

# Request interactive node
srun -p compute-hugemem -A merlab --nodes=1 --ntasks-per-node=1 --time=01:00:00 --mem=80G --pty /bin/bash

# activate conda angsd

conda activate angsd_env

# Specify paths and file names
DATADIR=/gscratch/scrubbed/elpetrou/pollock/angsd
MAFS_FILE=samples617_maf0.05_miss0.3.nuclear.mafs
SITES_FILE=$MAFS_FILE'.sites'

# Make the sites file
cd $DATADIR
gunzip $MAFS_FILE'.gz'

cut -f 1,2,3,4 $MAFS_FILE > $SITES_FILE 

# index the sites file
angsd sites index $SITES_FILE

