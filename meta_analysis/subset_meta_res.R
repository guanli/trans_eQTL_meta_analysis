#!/usr/bin/env Rscript
.libPaths("/net/snowwhite/home/guanli/transFusion/packrat/lib/x86_64-pc-linux-gnu/3.4.4")
library(optparse)
library(data.table)

#example usage 
#Rscript subset_meta_res.R -i ~/fusion_gtex_trans_meta/metal_res/meta_ENSG00000187642.5_1.tbl.gz
#-o ~/fusion_gtex_trans_meta/metal_res_subset/metal_subset_
option_list <- list(
  make_option(c("-i", "--input_metal_res"), type="character", action="store", default= NULL, 
              metavar="character", help="input metal result file"),
  
   
  make_option(c("-o", "--outpfx"), type="character", default= NULL, 
              help="output prefix")
  
  )

opt = parse_args(OptionParser(option_list=option_list))

meta = data.frame(fread(paste("zcat", opt$input_metal_res), sep="\t", header=T, stringsAsFactors=F))
cat('Before subsetting results has', nrow(meta),' rows.\n')
meta = subset(meta, meta$Direction=='--' |  meta$Direction=='++' | meta$Direction=='+-' | meta$Direction=='-+')
cat('After subsetting results has', nrow(meta),' rows.\n')

gene_name =  strsplit(opt$input_metal_res, "meta_", fixed = T)[[1]][2]
gene_name =  strsplit(gene_name, "_1.tbl", fixed = T)[[1]][1]
    
gzfh = gzfile(paste0(opt$outpfx, gene_name, '.tab.gz'), 'w')
write.table(meta, gzfh, row.names = F, col.names = T, quote=F, sep="\t")
close(gzfh)