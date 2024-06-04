 #!/bin/bash
 #BSUB -J my_test_job            # LSF job name
 #BSUB -o my_test_job.%J.out     # Name of the job output file 
 #BSUB -e my_test_job.%J.error   # Name of the job error file

 echo "this is a test"
 sleep 15
