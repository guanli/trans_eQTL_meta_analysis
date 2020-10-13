This pipeline 1) uses software METAL to perform QTL_meta-analysis, 2) extract results where the SNPs present in both studies and 
3) extract the snp-gene pairs that pass a pvalue threshold.

Command to run the pipeline <br /> 
`snakemake --snakefile meta final`

Parameters in the configuration file (meta.json). <br /> 
  -Required by METAL: `SCHEME`, `MARKERLABEL`, `PVALUELABEL`, `EFFECTLABEL`, `ALLELE`. <br /> 
  -Specify the folder where the within study eQTL mapping results are `STUDY_1_RES_DIR` and `STUDY_2_RES_DIR`.

Result files  <br />
  1. metal_res/ will have the output from METAL. The results of each gene are named meta_{gene_names}_1.tbl.gz and meta_{gene_names}_1.tbl.info.gz <br />
  2. metal_res_subset/ will have the results that are from SNPs present in both studies, which are obtained by subseting the resutls from 1. <br />
  3. metal_sig.txt will appear in the current folder, with the SNP-gene pairs that pass a pvalue threshold from 2. 
