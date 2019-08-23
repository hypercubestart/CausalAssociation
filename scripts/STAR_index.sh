#!/bin/bash

#$ -N genomeIndex
#$ -o /proj/omics4tb2/aliu/scratch/messages.genomeIndex.o.txt
#$ -e /proj/omics4tb2/aliu/scratch/messages.genomeIndex.e.txt
#$ -pe smp 32
#$ -S /bin/bash

cd /users/aliu
source .bash_profile 

echo "Original paths given my bash_profile..."
echo $PATH
echo ""


echo "Adding STAR to path..."
export PATH="/proj/omics4tb2/aliu/tools/STAR-2.5.2b/bin/Linux_x86_64/:$PATH"
echo $PATH
echo ""

STAR --runThreadN 32 --runMode genomeGenerate --genomeDir /proj/omics4tb2/aliu/projects/causalAssociation/data/genomeIndex --genomeFastaFiles /proj/omics4tb2/aliu/projects/causalAssociation/data/referenceGenome/GCF_000001405.39_GRCh38.p13_genomic.fna --sjdbGTFfile /proj/omics4tb2/aliu/projects/causalAssociation/data/annotation/GCF_000001405.39_GRCh38.p13_genomic.gff --sjdbOverhang 100 --sjdbGTFtagExonParentTranscript Parent

