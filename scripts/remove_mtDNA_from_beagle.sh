# Script name: remove_mtDNA_from_beagle.sh
# The purpose of this script is to remove SNPs in mtDNA from a beagle (or mafs) file created by angsd.

# Request compute time on an interactive node
srun -p compute-hugemem -A merlab --nodes=1 --ntasks-per-node=1 --time=02:00:00 --mem=40G --pty /bin/bash

# Specify the working directory and the data files 
DATADIR=/gscratch/scrubbed/elpetrou/pollock/angsd
MYINFILE=samples617_maf0.05_miss0.3.mafs # beagle file that contains the genotype likelihoods (without extension) 
MYOUTFILE=samples617_maf0.05_miss0.3.nuclear.mafs # name of the output file
MTDNA=NC_002081.1 #name of mtDNA scaffold from RefSeq Genome (GadMor3)

################################################################
# unzip the beagle file
cd $DATDIR
gunzip $MYINFILE'.gz'

# The mtDNA loci are at the end of the beagle file because of the order in which the chromosomes appear in the reference genome. 
# Count the number of lines in the beagle file that contain MtDNA
LINECOUNT=$(grep $MTDNA* $MYINFILE | wc -l)
echo $LINECOUNT

# Create a new beagle file that does not contain the lines at the end with the mtDNA loci
head -n -$LINECOUNT $MYINFILE > $MYOUTFILE

# Check to make sure there are no mtDNA loci in the output file
grep $MTDNA* $MYOUTFILE #ok, all clear, no hits!

# Gzip the beagle file you just created
gzip $MYOUTFILE