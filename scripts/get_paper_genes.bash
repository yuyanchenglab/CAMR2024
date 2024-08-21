#!/bin/bash

cd /project/hipaa_ycheng11lab/atlas/CAMR2024/

rm -f data/mentioned_genes.txt

while read gene; do
  grep -o -m 1 $gene data/CAMR_paper_text.txt >> data/mentioned_genes.txt
done < data/1_genes.txt
