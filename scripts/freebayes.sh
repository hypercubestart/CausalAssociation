#!/bin/bash

#$ -N freebayes
#$ -o /proj/omics4tb2/aliu/scratch/messages.freebayes.o.txt
#$ -e /proj/omics4tb2/aliu/scratch/messages.freebayes.e.txt
#$ -pe smp 32
#$ -S /bin/bash

cd /users/aliu
source .bash_profile 

python /proj/omics4tb2/aliu/projects/causalAssociation/src/freebayes.py
