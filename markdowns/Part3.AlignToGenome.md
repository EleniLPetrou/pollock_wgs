# Align to genome, remove PCR duplicates, clip overlapping sequences, and realign around indels 

## Download latest version of Atlantic cod genome from NCBI

I navigated to the NCBI RefSeq website and found the Atlantic cod genome, GadMor3: https://www.ncbi.nlm.nih.gov/genome/2661?genome_assembly_id=638047

  - Assembly information: GCA_902167405.1 gadMor3.0 scaffolds: 227 contigs: 1,442 N50: 1,015,663 L50: 169


I downloaded all files (fasta =GCF_902167405.1_gadMor3.0_genomic.fna.gz, gff = GCF_902167405.1_gadMor3.0_genomic.gff.gz) and uploaded them to Klone. I was working from a Windows laptop, so I used putty and the windows command terminal to do this:

```
pscp C:\Users\elpet\OneDrive\Documents\pollock\genome\* elpetrou@klone.hyak.uw.edu:/gscratch/merlab/genomes/atlantic_cod

```
 
## Install bwa-mem on Klone

I installed bwa v0.7.17 on Klone in a conda environment called bwa_env using the following commands:

``` bash

cd /gscratch/merlab/software/miniconda3/envs
conda create -n bwa_env
conda activate bwa_env
conda install -c bioconda bwa

```

## Index the reference genome with bwa
The first and basic step of running bwa is to index the genome. The basic usage of this command is:

```
bwa index [-a bwtsw|is] input_reference.fasta index_prefix

```
input_reference.fasta is an input file of sequence reads in fasta format, and index_prefix is the prefix of the generated index files. The option -a is required and can have two values: bwtsw (does not work for short genomes) and is (does not work for long genomes). Therefore, this value is chosen according to the length of the genome. I think we use the value -a bwtsw.


``` bash

#!/bin/bash
#SBATCH --job-name=pollock_bwaindex
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=2-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=bwa_env #name of the conda environment containing samtools software. 

## Specify directories with data
REFGENOME=/gscratch/merlab/genomes/atlantic_cod/GCF_902167405.1_gadMor3.0_genomic.fna #path to fasta genome
MYPREFIX=GCF_902167405.1_gadMor3.0_genomic
##################################################################
## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV

## Index the genome using bwa

bwa index -a bwtsw $REFGENOME $MYPREFIX

```
The files ending in ".bt2" are the indexed genome.


##################################################################################################
## Align sequences to reference genome with bowtie2

The command bowtie2 takes a Bowtie2 index and set of sequencing read files and outputs set of alignments in SAM format. I wrote a bash script containing the bowtie2 commands for aligning paired end data. Here it is:

```
!/bin/bash
##### ENVIRONMENT SETUP ##########
GENOMEDIR=/gscratch/merlab/genomes/atlantic_herring #location of genome
GENOME_PREFIX=GCF_900700415.1_Ch_v2.0.2 #prefix of .bt2 files made by bowtie2
SUFFIX1=_R1_001.trim.fastq # Suffix to trimmed fastq files. The forward reads with paired-end data.
SUFFIX2=_R2_001.trim.fastq # Suffix to trimmed fastq files. The reverse reads with paired-end data.

###################################

# Save the base name of each input file
MYBASE=$(basename --suffix=$SUFFIX1 "$1")

# Sanity check
echo "$1"
echo $MYBASE
echo ${MYBASE}$SUFFIX1
echo ${MYBASE}$SUFFIX2
echo ${MYBASE}.sam

# Run bowtie
bowtie2 -x $GENOMEDIR'/'$GENOME_PREFIX \
--phred33 -q \
-1 ${MYBASE}$SUFFIX1 \
-2 ${MYBASE}$SUFFIX2 \
-S ${MYBASE}.sam \
--very-sensitive \
--minins 0 --maxins 1500 --fr \
--threads ${SLURM_JOB_CPUS_PER_NODE} \
--rg-id ${MYBASE} --rg SM:${MYBASE} --rg LB:${MYBASE} --rg PU:Lane1 --rg PL:ILLUMINA

```

I did some tests and found that it took ~30 min to align a single sample to the genome, using 20 or 32 threads. Additionally, memory is not an issue for bowtie2 because it has a small memory footprint (3-4 Gb). Given that our node on Klone has 40 cores, I realized that I should stack "2 jobs" with 20 cores each and run them in parallel. To do this, I divided the fastq files into two folders and used a separate sbatch script to submit the job for each folder, like this : ` sbatch slurm_parallel_bowtie2.sh`. It's not the most elegant solution in the world but it reduced overall analysis time to 15 min per sample. 

Here is the sbatch script to do this:

```
#!/bin/bash
#SBATCH --job-name=elp_bowtie2_AK
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=6-8:00:00
## Memory per node
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/mmfs1/gscratch/scrubbed/elpetrou/fastq_trimmed/AK

## Job-specific Variables
###############################################
# Specify the path to your scipt and its name
MYSCRIPT=/mmfs1/home/elpetrou/scripts/parallel_bowtie2.sh

# Specify suffix of files to analyze
SUFFIX1=_R1_001.trim.fastq

################################################
# Make your script executable
chmod +x $MYSCRIPT

# Run the script on each file in the current directory
for i in *$SUFFIX1
do
	echo $i
	$MYSCRIPT $i
done

```
Next, I parsed the bowtie2 log files to estimate the average alignment rate across all samples. I did this quickly like this:

``` bash
LOGFILE=bowtie2_.out
OUFILE=bowtie2_unique_alignments.txt

grep -w "aligned concordantly exactly 1 time" $LOGFILE >> $OUTFILE
```

On average, it looks like ~86% of sequences are aligning to the Atlantic herring genome and 56% are aligning uniquely. I hope this translates into some nice data downstream!

## Use Samtools to convert sam to bam, format filter for quality & alignment, and sort bame files

I first had to install samtools on Klone, and that was a pain because I was having trouble with the version of openssl that samtools uses. To get around this issue, I decided to create a special conda environment for samtools v1.12. This how I did it:

```
# Install samtools using a conda environment (to get around the open ssl bug):
#The syntax for the first command says “conda” runs conda, “create -n” creates a new environment
conda create -n samtools_env

# To activate or enter the environments you just created simply type:
conda activate samtools_env

# Once in the environment, install the versions of htslib samtools and openssl that you want:
conda install -c bioconda htslib samtools openssl=1.0

```
I was able to run samtools, using this script (samtools_sbatch.sh):


```
#!/bin/bash
#SBATCH --job-name=elp_samtools
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
## Specify the directory containing data
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/bam #directory with sam files
SUFFIX1=.sam #file suffix
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=samtools_env #name of the conda environment containing samtools software. 

## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV


###################################################################################################################
## Move into the working directory and run script
cd $DATADIR

## Run samtools commands. This takes about 5 min per sample (so like 2 days total for whole data set?)
for MYSAMPLEFILE in *$SUFFIX1
do
    echo $MYSAMPLEFILE
    MYBASE=`basename --suffix=$SUFFIX1 $MYSAMPLEFILE`
    samtools view -bS -F 4 $MYBASE'.sam' > $MYBASE'.bam'
    samtools view -h -q 20 $MYBASE'.bam' | samtools view -buS - | samtools sort -o $MYBASE'_minq20_sorted.bam'
    samtools index $MYBASE'_minq20_sorted.bam'
done

## Flag explanations for samtools view:
## -b       output BAM
## -h       include header in SAM output
## -q INT   only include reads with mapping quality >= INT [0]
##-F INT   only include reads with none of the bits set in INT set in FLAG [0] (aka when this is set to 4, you remove unmapped reads)

## deactivate the conda environment
conda deactivate

```
## Use picard to remove PCR duplicates; use bamUtil to clip overlapping read pairs

In this step, I removed the PCR duplicates and trimmed the overlapping part of each read pair in pair-end data.

Install the software on Klone:(picard v2.18.7 and bamUtil v1.0.15)

``` 
conda create -n picard_env
conda activate picard_env
conda install -c bioconda picard bamutil
```

Run an sbatch script to conduct the analysis:

``` bash
#!/bin/bash
#SBATCH --job-name=elp_dedup_clip
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ####################################################
## Specify the directory containing data
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/bam #directory containing bam files
SUFFIX1=_minq20_sorted.bam # suffix to sorted and quality-filtered bam files produced in previous step of pipeline.
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=picard_env #name of the conda environment containing picard software. 

## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV

###### CODE FOR ANALYSIS ####################################################
## Move into the working directory
cd $DATADIR

## Run picard and bamutils (remove pcr duplicates and clip overlapping reads). These analyses take ~4 min to run per sample. Please note that picard can't handle newline breaks in the code (\) - so that is why all the commands are hideously written on one line.

for MYSAMPLEFILE in *$SUFFIX1
do
    echo $MYSAMPLEFILE
    MYBASE=`basename --suffix=$SUFFIX1 $MYSAMPLEFILE`
    echo $MYBASE
    picard MarkDuplicates INPUT=$MYBASE'_minq20_sorted.bam' OUTPUT=$MYBASE'_minq20_sorted_dedup.bam' METRICS_FILE=$MYBASE'_minq20_sorted_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
    bam clipOverlap --in $MYBASE'_minq20_sorted_dedup.bam' --out $MYBASE'_minq20_sorted_dedup_overlapclipped.bam' --stats
done 

## Leave the picard conda environment
conda deactivate

```
## Use GATK3 to realign bam reads around indels

Local realignment around indels allows us to correct mapping errors made by genome aligners and make read alignments more consistent in regions that contain indels.
To learn more about why we do indel realignment, follow this tutorial: https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md#section2

New versions of GATK take care of indel realignment during variant calling. Thus, I had to download GATK3 (version 3.8-1) from an online reporsitory to be able to use the IndelRealigner tool independently from variant calling (I want to do variant calling in angsd, and that program does not realign around indels). This is how I installed GATK3 on Klone:

``` bash
#Create a conda environment for GATK3
conda create -n gatk3_env

# enter the conda environment and download GATK3 from the web
conda activate gatk3_env
cd gatk3_env
wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
tar xjf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
 
#Test the installation:
GATK3=/gscratch/merlab/software/miniconda3/envs/gatk3_env/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
java -jar $GATK3 -h 

#This should print out some version and usage information, as well as a list of the tools included in the GATK. As the Usage line states, to use GATK you will always build your command lines like this: java -jar GenomeAnalysisTK.jar -T  [arguments] 


```

Next, I ran the pipeline for indel realignment. The first step is to identify potential indels using infromation from all bam files. This took ~72 hours to finish for all samples.  

``` bash
#!/bin/bash
#SBATCH --job-name=elp_RealignerTargetCreator
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=4-12:00:00
## Memory per node
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ####################################################
## Specify the directories and file names containing your data (edit lines 16-20 as needed)
DATADIR=/gscratch/scrubbed/elpetrou/bam #path to the bam files that you want to analyze with GATK3
GENOMEDIR=/gscratch/merlab/genomes/atlantic_herring #directory containing the genome
REFERENCE=GCF_900700415.1_Ch_v2.0.2_genomic.fna # Name of genome
BASEREFERENCE=GCF_900700415.1_Ch_v2.0.2_genomic #Name of genome without file extension
SUFFIX1=_minq20_sorted_dedup_overlapclipped.bam #Suffix of the bam files that you would like to analyze using GATK3

## Specify some information about the conda environments, singularities, and names of intermediate files. You probably do NOT need to edit this information.
BAMLIST=bam_list_dedup_overlapclipped.list # A list of merged, deduplicated, and overlap-clipped bam files. This file has to have a suffix of ".list"!! This list will be made in line 55 and will be saved to the $DATADIR
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
GATK3=/gscratch/merlab/software/miniconda3/envs/gatk3_env/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar # the path to the gatk3 jarfile
SAMTOOLS_ENV=samtools_env #name of the conda environment running samtools
PICARD_ENV=picard_env #name of the conda environment running picard
GATK3_ENV=gatk3_env #name of the conda environment running gatk3

###############################################################################
## Clean the environment before starting
module purge

## Source command will allow us to execute commands from a file in the current shell (conda)
source $MYCONDA

###### CODE FOR ANALYSIS ####################################################
## Use samtools to index the genome
conda activate $SAMTOOLS_ENV
samtools faidx $GENOMEDIR'/'$REFERENCE

# Make a text file containing a list of all the bam files you want to analyze
for MYSAMPLEFILE in $DATADIR'/'*$SUFFIX1
do
echo $MYSAMPLEFILE >> $BAMLIST
done

# Use samtools to index each bam file - this works!!
for MYSAMPLEFILE in $DATADIR'/'*$SUFFIX1
do
samtools index $MYSAMPLEFILE
done

#leave the samtools conda environment
conda deactivate 

###########################################
# activate the picard conda environment
conda activate $PICARD_ENV

# create a sequence dictionary for the reference genome (for some ridiculous reason, GATK3 needs this file)

picard CreateSequenceDictionary --REFERENCE $GENOMEDIR'/'$REFERENCE --OUTPUT $GENOMEDIR'/'$BASEREFERENCE.dict
#leave the picard conda environment
conda deactivate

##############################################
# activate the GATK3 conda environment
conda activate $GATK3_ENV
cd $DATADIR

# Create a list of potential indels
java -jar $GATK3 \
-T RealignerTargetCreator \
-R $DATADIR'/'$REFERENCE \
-I $DATADIR'/'$BAMLIST \
-nt ${SLURM_JOB_CPUS_PER_NODE} \
-o $DATADIR'/'all_samples_for_indel_realigner.intervals \
-drf BadMate

```

Next, each bam file was realigned around potential indels. This was slow (took 10 days to finish), as multithreading is not supported by IndelRealigner.

``` bash
#!/bin/bash
#SBATCH --job-name=elp_IndelRealigner
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=15-12:00:00
## Memory per node
#SBATCH --mem=300G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ####################################################
## Specify the directories and file names containing your data (edit lines 16-20 as needed)
DATADIR=/gscratch/scrubbed/elpetrou/bam #path to the bam files that you want to analyze with GATK3
GENOMEDIR=/gscratch/merlab/genomes/atlantic_herring #directory containing the genome
REFERENCE=GCF_900700415.1_Ch_v2.0.2_genomic.fna # Name of genome
BAMLIST=bam_list_dedup_overlapclipped.list # A list of merged, deduplicated, and overlap-clipped bam files. This file has to have a suffix of ".list"

## Specify some information about the conda environments and names of intermediate files. You probably do NOT need to edit this information.
GATK3=/gscratch/merlab/software/miniconda3/envs/gatk3_env/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar # the path to the gatk3 jarfile
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
GATK3_ENV=gatk3_env #name of the conda environment running gatk3


## Source command will allow us to execute commands from a file in the current shell (conda)
module purge
source $MYCONDA
conda activate $GATK3_ENV


########## Code for Indel Realignment ###################

## Move into the data directory so GATK3 writes output files here
cd $DATADIR 

## Run GATK3 IndelRealigner

java -jar $GATK3 \
-T IndelRealigner \
-R $GENOMEDIR'/'$REFERENCE \
-I $DATADIR'/'$BAMLIST \
-targetIntervals $DATADIR'/'all_samples_for_indel_realigner.intervals \
--consensusDeterminationModel USE_READS  \
--nWayOut _realigned.bam

## Leave gatk3 conda environment
conda deactivate

```

## Back up realigned bam files to LOLO Archive

Since bam realignment took such a long time, I decided to back up a copy of the realigned bam files (and the associated .bai files) to the LOLO archive: /archive/merlab/herring_wgs/realigned_bam

Explanation of tar command:

``` tar -zcvf archive-name.tar.gz source-directory-name ```

-z : Compress archive using gzip program in Linux or Unix

-c : Create archive on Linux

-v : Verbose i.e display progress while creating archive

-f : Archive File name

``` bash

# Compress all realigned bam files and transfer them to LOLO
# I think this will take around 24 hours to do
srun -p compute-hugemem -A merlab --nodes=1 --ntasks-per-node=1 --time=2-10:00:00 --mem=20G --pty /bin/bash

cd /gscratch/scrubbed/elpetrou/realigned_bam

# tar all the files in this directory
tar -zcvf realigned_bam_archive.tar.gz /gscratch/scrubbed/elpetrou/realigned_bam

# Transfer tarred archive from a hyak login node to LOLO
MYFILE=/mmfs1/gscratch/scrubbed/elpetrou/realigned_bam/realigned_bam_archive.tar.gz
TARGETDIR=elpetrou@lolo.uw.edu:/archive/merlab/herring_wgs/realigned_bam

scp $MYFILE $TARGETDIR



```
