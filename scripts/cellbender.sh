#!/bin/bash
#BSUB -J cellbender             # LSF job name
#BSUB -o cellbender%J.out      # Name of the job output file 
#BSUB -e cellbender.%J.error    # Name of the job error file
## BSUB -n 2                 # Request cores
## BSUB -M 50000              # JMEM=MB of memory for this job 
#BSUB -q gpu                 # Request GPU

source py311env/bin/activate

module load cuda/11.8

cellbender remove-background \
    --cuda \
    --input camr_scrublet_batch_filtered.h5ad \
    --output camr_scrublet_batch_filtered_bender.h5

