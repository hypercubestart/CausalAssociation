#!/bin/bash

#$ -N STARMapping
#$ -o /proj/omics4tb2/aliu/scratch/messages.STARMapping.o.txt
#$ -e /proj/omics4tb2/aliu/scratch/messages.STARMapping.e.txt
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

python /proj/omics4tb2/aliu/projects/causalAssociation/src/STAR_mapping.py
