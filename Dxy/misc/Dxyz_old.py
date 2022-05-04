import pandas as pd
import numpy as np
import sys

windows = [1000]
steps = [250]

rule all:
    input:
        expand('final/outlier_genes_{window}_{step}.bed', window=windows, step=steps)



rule run_dxyWindow_xy:
    input:
        'vrueppellii_ind66cov1pval5_majorminor5_snps.mafs',
        'vvNA_ind66cov1pval5_majorminor5_snps.mafs'
    output:
        'dxy_vRvvNA_{window}_{step}_snps'
    shell: '/home/joana/software/PopGenomicsTools/dxyWindow -winsize {wildcards.window} -stepsize {wildcards.step} {input[0]} {input[1]} > {output}'


rule run_dxyWindow_xz:
    input:
        'vrueppellii_ind66cov1pval5_majorminor5_snps.mafs',
        'vvEU_ind66cov1pval5_majorminor5_snps.mafs',
    output:
        'dxy_vRvvEU_{window}_{step}_snps'
    shell: '/home/joana/software/PopGenomicsTools/dxyWindow -winsize {wildcards.window} -stepsize {wildcards.step} {input[0]} {input[1]} > {output}'


rule run_makeDxyz_output:
    input:
        'dxy_vRvvNA_{window}_{step}_snps',
        'dxy_vRvvEU_{window}_{step}_snps'
    output:
        'dxyz_{window}_{step}_snps',
        'dxyz_{window}_{step}_snps.sorted',
        'final/outliers_{window}_{step}.bed'
    params:
        x=lambda wildcards: int(wildcards['window'])-1
    shell: """
    set +o pipefail;
    python /home/joana/scripts/mergeDxy2Dxyz.py {input[0]} {input[1]} {output[0]} &&
    cat {output[0]} | sort -V -k1,1 -k2,2 > {output[1]} &&
    cat {output[1]} | tail -n +2 | awk '($5 > {params.x})' | sort -g -k8,8 | head -n 1000   > {output[2]}
    """

rule intersect_outliers_and_genes:
    input:
        'final/outliers_{window}_{step}.bed',
        '/space/s1/joana/refgenomes/Functional_annotation/bedmap/genes.bed'
    output:
        temp('final/intersection_{window}_{step}.bed')
    shell:
        'bedtools intersect -a {input[0]} -b {input[1]} -loj > {output}'



rule filter_mapping:
    input: 'final/intersection_{window}_{step}.bed',
    output: temp('final/intersection_{window}_{step}.bed2'),
    run:
        in_df = pd.read_csv(input[0], '\t', header=None)
        out_df = in_df.groupby([0,1,2])[in_df.shape[1]-1].apply(lambda x: str(list(np.unique(list(x)))))
        out_df.to_csv(output[0], '\t', header=False)


rule merge_toFinal:
    input:
        'final/outliers_{window}_{step}.bed',
        'final/intersection_{window}_{step}.bed2'
    output:
        'final/outlier_genes_{window}_{step}.bed'
    run:
        df1 = pd.read_csv(input[0], '\t', header=None)
        df2 = pd.read_csv(input[1], '\t', header=None)
        result = pd.merge(df1,df2, how='inner', on=[0,1,2])
        result['3_y'] = result['3_y'].str.extractall('\((.*?)\)').groupby(level=0)[0].apply(lambda x: str(list(np.unique(list(x)))))
        result = result.fillna('[]')
        result.to_csv(output[0], '\t', header=False, index=False)


