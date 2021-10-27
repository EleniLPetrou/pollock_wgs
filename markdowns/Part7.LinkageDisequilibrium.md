# Analyses of linkage disequilibrium using ngsLD

I estimated the extent of linkage disequilibrium  using the genotype likelihoods estimated by angsd and the program ngsLD.

## Install ngsLD on Klone cluster
The following commands document how I installed ngsLD version 1.1.1 and its dependencies in a conda environment on our lab's Klone node on 20210726. 

``` bash
# Create a conda environment for ngsLD
conda create -n ngsLD_env

# enter the conda environment and download ngsLD from GitHub
conda activate ngsLD_env
cd ngsLD_env
git clone https://github.com/fgvieira/ngsLD.git

# install dependencies
#gcc: >= 4.9.2 tested on Debian 7.8 (wheezy)
conda install -c anaconda gcc_linux-64 # version 9.3.0 was available on conda and that is what I installed

#zlib: v1.2.7 tested on Debian 7.8 (wheezy). Version 1.2.11 was available on conda and that is what I installed
conda install -c conda-forge zlib

#gsl : v1.15 tested on Debian 7.8 (wheezy)
conda install -c conda-forge gsl ## version 2.7 was available on conda and that is what I installed

# PKG_CONFIG_PATH is a environment variable that specifies additional paths in which pkg-config will search for its .pc files.
#The pkg-config program is used to retrieve information about installed libraries in the system. 
#The primary use of pkg-config is to provide the necessary details for compiling and linking a program to a library

# To set the PKG_CONFIG_PATH value use the following command:
export PKG_CONFIG_PATH=/gscratch/merlab/software/miniconda3/envs/ngsLD_env/lib/pkgconfig

# Check
echo $PKG_CONFIG_PATH

# Compile ngsLD. The purpose of the make utility is to determine automatically which pieces of a large program need to be recompiled, and issue the commands to recompile them.
cd /gscratch/merlab/software/miniconda3/envs/ngsLD_env/ngsLD
make 

#Test the installation:
cd /gscratch/merlab/software/miniconda3/envs/ngsLD_env/ngsLD
./ngsLD

# Got this error message: ./ngsLD: error while loading shared libraries: libgsl.so.25: cannot open shared object file: No such file or directory

# Thus, I have to specify path to libgsl.so.25. NB: You will probably have to include this in all of the future scripts running ngsLD.

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gscratch/merlab/software/miniconda3/envs/ngsLD_env/lib 
export LD_LIBRARY_PATH

# test out the program again:
./ngsLD

# Success! It works!!
```

## Run LD analysis
I used ngsLD to estimate linkage disequilibrium in the herring data set. 

Important ngsLD parameters:

  --probs: specification of whether the input is genotype probabilities (likelihoods or posteriors)?
  
  --n_ind INT: sample size (number of individuals).
  
  --n_sites INT: total number of sites.
  
  --max_kb_dist DOUBLE: maximum distance between SNPs (in Kb) to calculate LD. Set to 0(zero) to disable filter. [100]
  
  --max_snp_dist INT: maximum distance between SNPs (in number of SNPs) to calculate LD. Set to 0 (zero) to disable filter. [0]
  
  --n_threads INT: number of threads to use. [1]
  
  --out FILE: output file name. [stdout]

Update as of 20210803: If you set --max_kb_dist to zero, this program will do all pairwise comparisons across all snps, which will create an absolutely massive file and take an eon to complete. I looked up the maximum chromosome length for the Atlantic herring genome on NCBI, https://www.ncbi.nlm.nih.gov/assembly/GCF_900700415.2#/st: and it is 33,084,258. I will use this value as the maximum distance between SNPs, so LD between SNPs on different chromsomes (with "Inf" distance separating them) does not get calculated.

``` bash
#!/bin/bash
#SBATCH --job-name=pollock_ngsLD
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=11-6:00:00
## Memory per node
#SBATCH --mem=300G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ##########
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=ngsLD_env #name of the conda environment containing samtools software. 

## Specify file names and data directories
DATADIR=/gscratch/scrubbed/elpetrou/pollock/angsd #path to input files from angsd
OUTDIR=/gscratch/scrubbed/elpetrou/pollock/ngsld #path to output files from ngsld
MYFILE=samples617_maf0.05_miss0.3.nuclear # name of input beagle file created by angsd without the .beagle.gz extension

## Specify the variables
N_IND=617 # number of individuals in beagle file
N_SITES=5315997 #number of variant sites in beagle file
N_THREADS=16 #make sure this matches the value in line 6 (--ntasks-per-node)
MAX_DIST=43798135 #maximum length of chromosome from here: https://www.ncbi.nlm.nih.gov/assembly/GCF_902167405.1/#/st

## Specify the path to required program (ngsLD) as variables
NGSLD=/gscratch/merlab/software/miniconda3/envs/ngsLD_env/ngsLD/

##################################################################
## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV

## Specify path to libgsl.so.25 file for ngsLD (do not change this!!)
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gscratch/merlab/software/miniconda3/envs/ngsLD_env/lib 
export LD_LIBRARY_PATH

################################################################
## Prepare the input files - ngsLD requires two input files:

##1)--geno FILE : a beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) can be inputted into ngsLD
## after the header row and the first three columns (i.e. positions, major allele, minor allele) are removed.

## Prepare the geno file

zcat $DATADIR/$MYFILE.beagle.gz | cut -f 4- | sed '1d' | gzip  > $OUTDIR/$MYFILE.formatted.beagle.gz

##2)--pos FILE: input file with site coordinates (one per line), where the 1st column stands for the chromosome/contig and 
## the 2nd for the position (bp). One convenient way to generate this is by selecting the first two columns of the mafs file outputted by ANGSD, 
## with the header removed. 

## Prepare a pos file - viera says to take the first two columns of mafs file and remove header: https://github.com/fgvieira/ngsLD/issues/4. That should do it.

zcat $DATADIR/$MYFILE.mafs.gz | cut -f 1,2 | sed '1d' | gzip > $OUTDIR/$MYFILE.formatted.pos.gz


################################################################
## Run ngsLD

$NGSLD/ngsLD \
--geno $OUTDIR/$MYFILE.formatted.beagle.gz \
--pos $OUTDIR/$MYFILE.formatted.pos.gz \
--probs \
--n_ind $N_IND \
--n_sites $N_SITES \
--max_kb_dist $MAX_DIST \
--n_threads $N_THREADS \
--out $OUTDIR/$MYFILE.ld

#Important ngsLD parameters:

#--probs: specification of whether the input is genotype probabilities (likelihoods or posteriors)?
#--n_ind INT: sample size (number of individuals).
#--n_sites INT: total number of sites.
#--max_kb_dist DOUBLE: maximum distance between SNPs (in Kb) to calculate LD. Set to 0(zero) to disable filter. [100]
#--max_snp_dist INT: maximum distance between SNPs (in number of SNPs) to calculate LD. Set to 0 (zero) to disable filter. [0]
#--n_threads INT: number of threads to use. [1]
#--out FILE: output file name. [stdout]


# Update as of 20210803: If you set --max_kb_dist to zero, this program will do all pairwise comparisons across all snps, which will create an absolutely massive file and take an eon to complete. I looked up the maximum chromosome length for the Atlantic herring genome on NCBI, https://www.ncbi.nlm.nih.gov/assembly/GCF_900700415.2#/st: 33,084,258. I will use this value as the maximum distance between SNPs, so LD between SNPs on different chromsomes, with "Inf" distance separating them, does not get calculated. 
```
