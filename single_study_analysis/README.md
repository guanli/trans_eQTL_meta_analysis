This pipeline uses QTLtools to do the following two analyses.  <br />
  1. Nomimal trans-eQTL pass  <br />
  2. Approximate permutation pass to adjust for the number of tested SNPs and the number of tested genes.

Before running this pipeline <br />
  - Create a folder called data and place the input files for QTLtools (gene expession, genotype and covariate) there  <br />
  - Creat a folder called scripts/ and place the two scripts from QTLtools website (qtltools-check_beta_approx.R and qtltools-runFDR_cis.R) <http:https://qtltools.github.io/qtltools/>
  
Command to run the pipeline  <br />
  `snakemake --snakefile Snakefile --configfile config_analysis.json`. <br />
  The window size in the Snakefile is set to be 1000000, larger than the size of any chromosome in order to get inter-chromosomal SNP-gene pairs only. 

Options in the configuration file <br />
  - `TRANS_NOMINAL_CUTOFF` only SNP-gene pairs that pass this pvalue threshold will be saved. Please note that the result file will be extremely LARGE
  if saving all of the results regardless of pvalues.<br />
  - `PARAM_QTLTOOLS_APX_NULL_SMP_NUM` The number of randomly selected genes to run permutation tests for <br />
  - `PARAM_QTLTOOLS_ADJ_GENE_OUT` Output file name for the output after correcting for the number of tested genes 

Result files  <br />
  - Nominal pass: trans_qtltools_nominal.tsv.gz in the current folder.  <br />
  - Approximate permutation pass: the following files in qtltools_apx_null/  <br />
      - *best.txt.gz
      - *bins.txt.gz
      - *hits.txt.gz
      - A file whose name is specified in `PARAM_QTLTOOLS_ADJ_GENE_OUT`
  - QTL mapping results whose multiple testing correction are done based on results of the approximate permutation pass: the following files in trans_final_qtltools_nominal/
  
