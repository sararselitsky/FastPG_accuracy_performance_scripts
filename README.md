# FastPG accuracy performance scripts
Scripts used in the manuscript "FastPG: Fast clustering of millions of single cells"
https://www.biorxiv.org/content/10.1101/2020.06.19.159749v1

Running fastpg_compare.R:
Need to download gold standard data and change appropriate directories.

To run: Rscript fastpg_compare.R cell_number iterations output_name

Example: Rscript fastpg_compare.R 10000 10 test.txt

Authors:
- FastPG vs. PhenoGraph accuracy: Sara Selitsky, sararselitsky@gmail.com
- FastPG vs. PhenoGraph time: Sara Selitsky and Tom Bodenheimer, bodenhei@email.unc.edu
- scRNA-seq analysis: Siyao Liu, siyao@email.unc.edu and Sara Selitsky

