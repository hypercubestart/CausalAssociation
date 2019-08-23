#!/bin/bash

#$ -N splitFastq
#$ -o /proj/omics4tb2/aliu/scratch/messages.splitFastq.o.txt
#$ -e /proj/omics4tb2/aliu/scratch/messages.splitFastq.e.txt
#$ -pe smp 32
#$ -S /bin/bash

cd /users/aliu
source .bash_profile 

python /proj/omics4tb2/aliu/projects/causalAssociation/src/splitFastq.py
