#!/bin/bash

#$ -N both
#$ -o /proj/omics4tb2/alomana/scratch/messages.cr.count.both.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.cr.count.both.e.txt
#$ -pe smp 40
#$ -l hostname=baliga2
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

echo "Running count..." 
cd /proj/omics4tb2/alomana/projects/i18/results/
time cellranger count --id=both_PSTAT_d4_a --transcriptome=/proj/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-hg19-and-mm10-3.0.0 --fastqs=/proj/omics4tb2/alomana/projects/i18/results/HTTLKBGX9/outs/fastq_path --sample=PSTAT_d4_a --expect-cells=1000 --localcores=40 --localmem=90
time cellranger count --id=both_PSTAT_d4_b --transcriptome=/proj/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-hg19-and-mm10-3.0.0 --fastqs=/proj/omics4tb2/alomana/projects/i18/results/HTTLKBGX9/outs/fastq_path --sample=PSTAT_d4_b --expect-cells=1000 --localcores=40 --localmem=90
echo "... run completed."
