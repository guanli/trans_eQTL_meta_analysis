#!/usr/bin/env snakemake

"""
Snakemake QTL pipeline
----------------------

Snakemake pipeline for molecular trait QTL mapping
"""

def chunker(seq, size):
    """
    Chunks a big list (<seq>) into smaller groups, defined by <size>. 
    """
    res = []
    for el in seq:
        res.append(el)
        if len(res) == size:
            yield res
            res = []
    if res:
        yield res


rule trans_qtl: #now trans means > 1M
    input:
        # QTLtools nominal associations
        'trans_qtltools_nominal.tsv.gz',
        
        # QTLtools approximate permutation associations
        #expand('qtltools_apx_null/{out_put}', out_put=config['PARAM_QTLTOOLS_ADJ_GENE_OUT'])
        

# QTLtools: file prep ##########################################################
rule qtltools_prep__vcf:
    """
    QTLtools:   generate tabix index of vcf input files 
    """
    input:
        'data/genotypes.vcf.gz'
    output:
        'data/genotypes.vcf.gz.tbi'
    shell:
        'tabix {input}'


rule qtltools_prep__bed:
    """
    QTLtools:   generate tabix index of bed input files 
    """
    input:
        'data/moltraits_trans_qtltools.bed.gz'
    output:
        'data/moltraits_trans_qtltools.bed.gz.tbi'
    shell:
        'tabix {input}'
################################################################################


# QTLtools: nominal pass #######################################################
rule qtltools_proximal__nominal:
    """
    QTLtools:   nominal proximal mapping (association with all variants)
    """
    input:
        geno='data/genotypes.vcf.gz',
        geno_tbi='data/genotypes.vcf.gz.tbi',
        pheno='data/moltraits_trans_qtltools.bed.gz',
        pheno_tbi='data/moltraits_trans_qtltools.bed.gz.tbi'
    output:
        temp('trans_swarm_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.best.txt.gz'),
        temp('trans_swarm_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.bins.txt.gz'),
        temp('trans_swarm_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.hits.txt.gz')
    params:
        other_param=config['TRNAS_PARAM_QTLTOOLS'],
        nominal_cutoff=config['TRANS_NOMINAL_CUTOFF'],
        out_prefix = 'trans_swarm_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}'
    shell:
        'QTLtools trans --vcf {input.geno} --bed {input.pheno} '
            '--nominal --threshold {params.nominal_cutoff} ' 
            #'--normal ' # force the input phenotype to N(0, 1) 
            '--seed 15112011 ' # seed for random number generator
            '--window 1000000 ' # window size. default = 1Mb
            '--chunk {wildcards.j_cur} {wildcards.j_total} ' # for swarming
            '{params.other_param} ' # PARAM from config
            '--out {params.out_prefix}'

rule unzip_hits_gz:
    """
      unzip the swarm_qtltools_nominal/qtltools_nominal-{j_cur}_{j_total}.txt.hits.txt.gz file
      so that the next in-place zip can work
    """
    input:
        'trans_swarm_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.hits.txt.gz'
    output:
        'trans_swarm_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.hits.txt' 
    shell:
        'gzip -d {input}'

rule qtltools_proximal__nominal_concat:
    """
    QTLtools: swarm proximal mapping & combine the output
    
    Notes:
    - 1 based so first iteration is 1 not 0.
    - removes the temp input data
    """
    input:
        expand(
            'trans_swarm_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.hits.txt', 
            j_cur=range(1, config['TRANS_NJOBS_NOMINAL']+1), 
            j_total=config['TRANS_NJOBS_NOMINAL']
        )
    output:
        'trans_qtltools_nominal.tsv.gz',

    shell:
        # add the header to the output file
        'echo "pheno_id pheno_chr pheno_start '
            'var_id var_chr var_pos '
            'p_nominal dummy beta" | '
        'sed s/" "/"\t"/g | gzip -c > {output}; '
            
        # cat swarm files to make the output file
        'cat {input} | sed s/" "/"\t"/g | gzip -c >> {output}; '
        
        #rezip the hits.txt files
        'gzip {input}'
        
################################################################################
# QTLtools: approximate permutation pass #######################################################
rule approximate_null:
    """
    QTLtools:   1st step in trans approx pass -- build null distribution
    """
    input:
        geno='data/genotypes.vcf.gz',
        geno_tbi='data/genotypes.vcf.gz.tbi',
        pheno='data/moltraits.bed.gz',
        pheno_tbi='data/moltraits.bed.gz.tbi'
    output:
       best= expand('qtltools_apx_null/{outpfx}.best.txt.gz', outpfx = config['PARAM_QTLTOOLS_APX_NULL_OUT']),
       bins= expand('qtltools_apx_null/{outpfx}.bins.txt.gz', outpfx = config['PARAM_QTLTOOLS_APX_NULL_OUT']),
       hits= expand('qtltools_apx_null/{outpfx}.hits.txt.gz', outpfx = config['PARAM_QTLTOOLS_APX_NULL_OUT'])
    params:
        other_param=config['PARAM_QTLTOOLS'],
        smp_num=config['PARAM_QTLTOOLS_APX_NULL_SMP_NUM'],
        outpfx=config['PARAM_QTLTOOLS_APX_NULL_OUT']
    shell:
        'QTLtools trans --vcf {input.geno} --bed {input.pheno} '
            '--sample {params.smp_num} ' 
            '--normal ' # force the input phenotype to N(0, 1) 
            '--window 1000000 ' # window size. default = 1Mb
            '{params.other_param} ' # PARAM from config
            '--out qtltools_apx_null/{params.outpfx}'
            
rule adjust_variants:
    """
    QTLtools:   2nd step in trans approx pass -- Adjust the nominal P-values regarding the number of tested variants
    """
    input:
        geno='data/genotypes.vcf.gz',
        geno_tbi='data/genotypes.vcf.gz.tbi',
        pheno='data/moltraits.bed.gz',
        pheno_tbi='data/moltraits.bed.gz.tbi',
        adjust=rules.approximate_null.output.best 
    output:
       best= expand('qtltools_apx_null/{outpfx}.best.txt.gz', outpfx = config['PARAM_QTLTOOLS_ADJ_VAR_OUT']),
       bins= expand('qtltools_apx_null/{outpfx}.bins.txt.gz', outpfx = config['PARAM_QTLTOOLS_ADJ_VAR_OUT']),
       hits= expand('qtltools_apx_null/{outpfx}.hits.txt.gz', outpfx = config['PARAM_QTLTOOLS_ADJ_VAR_OUT'])
    params:
        adj_var_thresh=config['PARAM_QTLTOOLS_APX_NULL_ADJ_VAR_THRESH'],
        outpfx=config['PARAM_QTLTOOLS_ADJ_VAR_OUT']
    shell:
        'QTLtools trans --vcf {input.geno} --bed {input.pheno} '
            '--adjust {input.adjust} ' 
            '--normal ' # force the input phenotype to N(0, 1) 
            '--threshold {params.adj_var_thresh} ' # window size. default = 1Mb
            '--out qtltools_apx_null/{params.outpfx}'
            
rule adjust_genes:
    """
    QTLtools:   3rd step in trans approx pass -- Adjust the nominal P-values regarding the number of tested genes
    """
    input:
        best = rules.adjust_variants.output.best,
        hits = rules.adjust_variants.output.hits
    output:
        expand('qtltools_apx_null/{out_put}', out_put=config['PARAM_QTLTOOLS_ADJ_GENE_OUT'])
    params:
        QTLtools_adj_gene_script=config['PARAM_QTLTOOLS_ADJ_GENE_SCRIPT'],
        adj_gene_thresh=config['PARAM_QTLTOOLS_APX_NULL_ADJ_GENE_THRESH'],
        output=config['PARAM_QTLTOOLS_ADJ_GENE_OUT']
    shell:
        'Rscript {params.QTLtools_adj_gene_script} {input.best} {input.hits} {params.adj_gene_thresh} {output}'
    
            
################################################################################
# QTLtools: final nominal pass #######################################################        
rule qtltools_distal__nominal:
    """
    QTLtools:  final nominal distall mapping (after #peer factors is chosen)
    """
    input:
        geno='data/genotypes.vcf.gz',
        geno_tbi='data/genotypes.vcf.gz.tbi',
        pheno='data/moltraits_trans_qtltools.bed.gz',
        pheno_tbi='data/moltraits_trans_qtltools.bed.gz.tbi'
    output:
        'trans_final_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.best.txt.gz',
        'trans_final_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.bins.txt.gz',
        'trans_final_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.hits.txt.gz'
    params:
        other_param=config['TRNAS_PARAM_QTLTOOLS'],
        nominal_cutoff=config['TRANS_NOMINAL_CUTOFF'],
        out_prefix = 'trans_final_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}'
    shell:
        'QTLtools trans --vcf {input.geno} --bed {input.pheno} '
            '--nominal --threshold {params.nominal_cutoff} ' 
            #'--normal ' # force the input phenotype to N(0, 1) 
            '--seed 15112011 ' # seed for random number generator
            '--window 1000000 ' # window size. default = 1Mb
            '--chunk {wildcards.j_cur} {wildcards.j_total} ' # for swarming
            '{params.other_param} ' # PARAM from config
            '--out {params.out_prefix}'

rule qtltools_distal__nominal_final:
    input:
        expand('trans_final_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.best.txt.gz',j_cur=range(1, config['TRANS_NJOBS_NOMINAL_FINAL']+1),j_total=config['TRANS_NJOBS_NOMINAL_FINAL']),      
        expand('trans_final_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.bins.txt.gz',j_cur=range(1, config['TRANS_NJOBS_NOMINAL_FINAL']+1),j_total=config['TRANS_NJOBS_NOMINAL_FINAL']),
        expand('trans_final_qtltools_nominal/trans_qtltools_nominal-{j_cur}_{j_total}.hits.txt.gz',j_cur=range(1, config['TRANS_NJOBS_NOMINAL_FINAL']+1),j_total=config['TRANS_NJOBS_NOMINAL_FINAL'])



