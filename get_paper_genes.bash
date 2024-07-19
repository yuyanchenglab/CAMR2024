#!/bin/bash

echo "" > mentioned_genes.txt

while read gene; do
  grep -o -m 1 $gene CAMR_paper_text.txt >> mentioned_genes.txt
done < CAMR_genes.csv
