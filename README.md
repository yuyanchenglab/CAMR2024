Data Paper: https://www.biorxiv.org/content/10.1101/2024.01.24.577060v1.full

Data Download Link: https://cellxgene.cziscience.com/collections/a0c84e3f-a5ca-4481-b3a5-ccfda0a81ecc

Original Data Source: d0183df5-815d-48c2-bcfe-fbf9b716505c.h5ad
Published Data: 0ca84fa3-eaa3-456e-a24c-e13e225d7ba6.h5ad

Yuyan Org Notes: I need you to go back to organize your code and output. Rename your output files (designer notes excel sheet and dotplots) as V1-V3. Have a readme file on what each version means. I believe V1 : all genes queried and curated; no filtering; V2: queried and curated genes that are manually labeled as 'great'; V3: all genes that are labeled as 'great' on 09/16; and V4: Final set lableld by Yuyan (edited) 

I think we still need V3 as I'll need a shorter 'great' list and their dotplot to look at...thanks

Jeff Org Notes:

Everything up to script 4 should be ok.

Scripts 0.5, 05.1, and 10_Rename_Genes should be combined so that all additional tweaking of the adata object can be done in one spot.

Script 06 should be ok, but parts should be moved to to the adata update.

Scripts under 07 should have a good amount of notes written about it.

Scripts 08 and 08.1 should be combined and should be delineated from 09.

Scripts 11, 12, 13 & 14 should be combined such that each of them is the result of an if statement in the script itself.

I also need to delete portions that are no longer applicable or poorly named
