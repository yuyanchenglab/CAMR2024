#!/bin/bash
#BSUB -J cellbender             # LSF job name
#BSUB -o cellbender%J.out      # Name of the job output file 
#BSUB -e cellbender.%J.error    # Name of the job error file
#BSUB -n 5                 # Request cores
#BSUB -M 200000              # JMEM=MB of memory for this job 
#BSUB -R "rusage[mem=200000] span[hosts=1]"    # Make sure the job's assigned node has enough memory and that all cores are on the same worker node
##BSUB -q gpu                 # Request GPU
#BSUB -notify done         # Whenever the job finishes, email
#BSUB -u chena23@seas.upenn.edu # User's email
# #BSUB -m node123              # Specify node

# :%s/JNAME/NEW_JNAME/g
# :%s/JMEM/NEW_JMEM/g

source py311env/bin/activate

cellbender remove-background \
    --input camr_scrublet_batch_filtered.h5ad \
    --output camr_scrublet_batch_filtered_bender.h5

