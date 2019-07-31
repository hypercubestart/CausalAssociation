#!/bin/bash

#$ -N count
#$ -o /proj/omics4tb2/alomana/scratch/messages.deconvolution.mouse.count.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.deconvolution.mouse.count.e.txt
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
cd /proj/omics4tb2/alomana/projects/i18/results/deconvolution/mouse
time cellranger count --id=d2_dmso_mouse --transcriptome=/proj/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0 --fastqs=/proj/omics4tb2/alomana/projects/i18/results/deconvolution/HFGF7BGXB/outs/fastq_path --sample=d2_dmso --expect-cells=1000 --localcores=40 --localmem=90
time cellranger count --id=d4_dmso_mouse --transcriptome=/proj/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0 --fastqs=/proj/omics4tb2/alomana/projects/i18/results/deconvolution/HFGF7BGXB/outs/fastq_path --sample=d4_dmso --expect-cells=1000 --localcores=40 --localmem=90
time cellranger count --id=d4_pstat_r2_mouse --transcriptome=/proj/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0 --fastqs=/proj/omics4tb2/alomana/projects/i18/results/deconvolution/HFGF7BGXB/outs/fastq_path --sample=d4_pstat_r2 --expect-cells=1000 --localcores=40 --localmem=90
echo "... run completed."

# 1,control,untreated samples,NA
# 2,d2_pstat,2 day treatment with pitavastatin,NA
# 3,d2_dmso,2 day treatment with dmso vehicle,WT mice macrophage added (10% v:v)
# 4,d3_pstat,3 day treatment with pitavastatin,NA
# 5,d3_dmso,3 day treatment with dmso vehicle,NA
# 6,d4_pstat_r1,4 day treatment with pitavastatin technical rep1,NA
# 7,d4_dmso,4 day treatment with dmso vehicle,MAVK mice macrophage added (10% v:v)
# 8,d4_pstat_r2,4 day treatment with pitavastatin technical rep2,NLRP3 mice macrophage added (10% v:v)
