#!/bin/bash
#BSUB -J scrublet_batch             # LSF job name
#BSUB -o logs/scrublet_batch%J.out      # Name of the job output file 
#BSUB -e logs/scrublet_batch.%J.error    # Name of the job error file
#BSUB -n 5                 # Request cores
#BSUB -M 400000              # JMEM=MB of memory for this job 
#BSUB -R "rusage[mem=400000] span[hosts=1]"    # Make sure the job's assigned node has enough memory and that all cores are on the same worker node
#BSUB -notify done         # Whenever the job finishes, email
#BSUB -u chena23@seas.upenn.edu # User's email
# #BSUB -m node123              # Specify node

# :%s/JNAME/NEW_JNAME/g
# :%s/JMEM/NEW_JMEM/g

cd /project/hipaa_ycheng11lab/atlas/CAMR2024

mkdir -p 00_scrublet_batch

python scripts/scrublet_batch.py