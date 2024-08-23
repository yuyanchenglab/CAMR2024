#!/bin/bash

cd /project/hipaa_ycheng11lab/atlas/CAMR2024/

rm -f 01_QualityControl/mentioned_genes.txt

while read gene; do
  grep -o -m 1 $gene data/CAMR_paper_text.txt >> 01_QualityControl/mentioned_genes.txt
done < 01_QualityControl/1_genes.txt
