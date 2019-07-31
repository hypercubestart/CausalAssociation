#!/bin/bash

#$ -N makefastq
#$ -o /proj/omics4tb2/alomana/scratch/messages.deconvolution.makefastq.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.deconvolution.makefastq.e.txt
#$ -pe smp 40
#$ -l hostname=baliga1
#$ -S /bin/bash

cd /users/alomana
source .bash_profile 

echo "Original paths given my bash_profile..."
echo $PATH
echo ""

echo "Adding extra paths for Cell Ranger..."
source /proj/omics4tb2/alomana/software/cellRanger/cellranger-3.0.2/sourceme.bash
echo $PATH
echo ""

echo "Adding path to bcl2fastq..."
export PATH="/usr/local/bcl2fastq-v2.20.0.422/bin:$PATH"
echo $PATH
echo ""

echo "where am i?"
uname -a
echo ""

echo "Running makefastq..." 
cd /proj/omics4tb2/alomana/projects/i18/results/deconvolution
time cellranger mkfastq --run=/proj/omics4tb2/alomana/projects/i18/data/deconvolution/original/190510_NS500720_0419_AHFGF7BGXB --samplesheet=/proj/omics4tb2/alomana/projects/i18/data/deconvolution/cellranger-bcl-samplesheet.csv --localcores=40 --localmem=90
