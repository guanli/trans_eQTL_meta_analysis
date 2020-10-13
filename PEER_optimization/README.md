This pipeline finds the number of PEER factors that optimizes trans-QTL detection rate. <br />
  1. Use EPACT to run genome-wide associations between SNPs and each PEER factors to identify which PEER factors are under genetic control.  <br />
  2. Performs LD-pruning to the input genotype file <br />
  3. Use the pruned list of genetic variants, the input gene expression file and the input covariate file to conduct QTL mapping <br />
     using each of the PEER factors specified with the PEER factors that are under genetic control removed. <br />
  4. Generate QTL mapping results into ../peer_runs/factor_k/ folder. <br />
Therefore, by comparing the number of genes having eQTL or the number of significant SNP-gene associations, the number of PEER factors 
  that optimizes QTL detection rate can be determined.
  
Run the Snakemake pipeline with the following command. <br />
`snakemake --snakefile peer_opt_trans --configfile peer_opt_trans.json`

Output in the current folder  <br />
  -trans_eGene_peer_factors-qtltools_nominal-summary.tsv.gz: the number of genes with QTL from result using different number of PEER factors <br />
  -trans_eGene_peer_factors-qtltools_nominal-summary.pdf: a line plot displaying the above results
