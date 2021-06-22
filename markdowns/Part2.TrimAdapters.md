# Trim adapters using Trimmomatic

I will trim adapters using [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

Here are some highlights from the user manual:

Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop Illumina (FASTQ) data as well as to remove adapters.Trimmomatic works with FASTQ files (using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used). 

## Running Trimmomatic

### Processing Order
The different processing steps occur in the order in which the steps are specified on the command line. It is recommended in most cases that adapter clipping, if required, is done as early as possible.

### Paired End Mode Input and Output files
For paired-end data, two input files, and 4 output files are specified, 2 for the 'paired' output where both reads survived the processing, and 2 for corresponding 'unpaired' output where a read survived, but the partner read did not.

#### Syntax: 
See [user manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for explanation of terms
``` 
PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] #[-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>]
```

### Trimming commands:

  - ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
  - MINLEN: Drop the read if it is below a specified length


#### Syntax:
See [user manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for explanation of terms
```
ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
```

### File with Nextera adapters, as specified in Physalia course and Illumina Adapter Sequences document:

I made a [fasta file containing Illumina Nextera adapter sequences](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/main/Scripts/NexteraPE_EP.fa). These are the sequences that I will trim away from my raw sequencing data using Trimmomatic. 

### Script to run Trimmomatic over a directory of fastq files on Klone
 
Next, I wrote an sbatch script to run trimmomatic on Klone

This script removes Illumina adapters, trims sequences if their phred score drops below 20 (SLIDINGWINDOW:4:20), and removes sequences that are shorter than 40 bp (MINLEN:40):

```
#!/bin/bash
#SBATCH --job-name=pollock_trim
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
DATADIR=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3
SAMPLELIST=trimmomatic_list.txt
SUFFIX1=_R1_001.fastq # Suffix to raw fastq files. The forward reads with paired-end data.
SUFFIX2=_R2_001.fastq # Suffix to raw fastq files. The reverse reads with paired-end data. 
ADAPTERFILE=/mmfs1/home/elpetrou/scripts/NexteraPE_EP.fa # File with adapter sequences. This is based on the file that was provided in the Physalia course, and double-checked using the Illumina Adapter Sequences pdf
OUTDIR=/gscratch/scrubbed/elpetrou/pollock/fastq_trimmed #where to store output files


##############################################################################
## Save the sample names of the forward reads into a text file (for looping thru samples later)
mkdir $OUTDIR

cd $DATADIR
ls *$SUFFIX1 > $SAMPLELIST

##############################################################################
## Run Trimmomatic. For an explanation of trimmomatic commands, see: 
## http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
## https://datacarpentry.org/wrangling-genomics/03-trimming/index.html


for infile in `cat $SAMPLELIST`
do
        base=`basename --suffix=$SUFFIX1 $infile`
        trimmomatic PE \
        -threads ${SLURM_JOB_CPUS_PER_NODE} -phred33 \
        ${infile} \
        ${base}$SUFFIX2 \
        ${base}_R1_001.trim.fastq \
        ${base}_R1_001.unpaired.trim.fastq \
        ${base}_R2_001.trim.fastq \
        ${base}_R2_001.unpaired.trim.fastq \
        ILLUMINACLIP:$ADAPTERFILE:2:30:10:1:true \
        SLIDINGWINDOW:4:20 \
        MINLEN:40
done

## Yaaay! This code is running!!
## Here is an example of some output: Input Read Pairs: 6825779 Both Surviving: 5932557 (86.91%) Forward Only Surviving: 468865 (6.87%) Reverse Only Surviving: 215887 (3.16%) Dropped: 208470 (3.05%)
## For one set of fastq samples, I think it took about 5 min. So it will take about a day to finish trimming one lane on Klone (give it two days as buffer)
#############################################################################
## Move the results files to the output directory

mv *_R1_001.trim.fastq $OUTDIR
mv *_R2_001.trim.fastq $OUTDIR
mv *_R1_001.unpaired.trim.fastq $OUTDIR
mv *_R2_001.unpaired.trim.fastq $OUTDIR


```

### Let's look at the trimmomatic output:

I navigated into the folder containing the log files produced by trimmomatic and I ran multiQc to get data summaries. 

``` bash

srun -p compute-hugemem -A merlab --nodes=1 --ntasks-per-node=1 --time=02:00:00 --mem=50G --pty /bin/bash

conda activate multiqc_env

cd ~/pollock_scripts

multiqc .

```

On average, we retained 88% of bases after trimming sequencing adapters and poor-quality bases. 





