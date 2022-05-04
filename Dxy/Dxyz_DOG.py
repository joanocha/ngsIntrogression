import pandas as pd
import numpy as np
import sys

windows = [5000]
steps = [1250]

rule all:
    input:
        expand('final/outlier_genes_{window}_{step}.bed', window=windows, step=steps)



rule run_dxyWindow_xy:
    input:
        'vzerda_intersection.maf.gz',
        'vrueppellii_intersection.maf.gz'
    output:
        'dxy_vZvR_{window}_{step}_snps'
    shell: '/home/joana/software/PopGenomicsTools/dxyWindow -winsize {wildcards.window} -stepsize {wildcards.step} {input[0]} {input[1]} > {output}'


rule run_dxyWindow_xz:
    input:
        'vzerda_intersection.maf.gz',
        'vvulpes_intersection.maf.gz',
    output:
        'dxy_vZvV_{window}_{step}_snps'
    shell: '/home/joana/software/PopGenomicsTools/dxyWindow -winsize {wildcards.window} -stepsize {wildcards.step} {input[0]} {input[1]} > {output}'


rule run_makeDxyz_output:
    input:
        'dxy_vZvR_{window}_{step}_snps',
        'dxy_vZvV_{window}_{step}_snps'
    output:
        'dxyz_{window}_{step}_snps'
    wildcard_constraints:
        step="\d+"
    shell: """
    python /home/joana/scripts/mergeDxy2Dxyz.py {input[0]} {input[1]} {output[0]}
    """

rule run_makeDxyz_output_pt2:
    input:
        'dxyz_{window}_{step}_snps'
    output:
        'final/outliers_{window}_{step}.bed'
    shell: """
    cat {input} | tail -n +2  > {output}
    """
rule fillna:
    input:
        'final/outliers_{window}_{step}.bed'
    output:
        'final/outliers_{window}_{step}.fillNA.bed'
    run:
        df = pd.read_csv(input[0], '\t', header=None, names = ['scaffold', 'start', 'end', ''])
        df = df.fillna('NA')
        df = df[]
        df.to_csv(output[0], '\t', header=False, index=False)



rule sort_bed_file:
    input: 'final/outliers_{window}_{step}.fillNA.bed',
    output: 'final/outliers_{window}_{step}.sorted.bed'
    shell:
         'bedtools sort -i {input} > {output}'


rule map_genes:
    input:
        'final/outliers_{window}_{step}.sorted.bed',
        '/space/s2/joana/refgenomes/Functional_annotation/bedmap/genesDog.sorted.bed'
    output:
        'final/outlier_genes_{window}_{step}.bed'
    shell:
        'bedtools map -a {input[0]} -b {input[1]} -c 4 -o distinct | sort -g -k8,8  > {output}'



#rule filter_mapping:
#    input: 'final/intersection_{window}_{step}.bed',
#    output: temp('final/intersection_{window}_{step}.bed2'),
#    run:
#        in_df = pd.read_csv(input[0], '\t', header=None)
#        out_df = in_df.groupby([0,1,2])[in_df.shape[1]-1].apply(lambda x: str(list(np.unique(list(x)))))
#        out_df.to_csv(output[0], '\t', header=False)
#
#
#rule merge_toFinal:
#    input:
#        'final/outliers_{window}_{step}.bed',
#        'final/intersection_{window}_{step}.bed2'
#    output:
#        'final/outlier_genes_{window}_{step}.bed'
#    run:
#        df1 = pd.read_csv(input[0], '\t', header=None)
#        df2 = pd.read_csv(input[1], '\t', header=None)
#        result = pd.merge(df1,df2, how='inner', on=[0,1,2])
#        result['3_y'] = result['3_y'].str.extractall('\((.*?)\)').groupby(level=0)[0].apply(lambda x: str(list(np.unique(list(x)))))
#        result = result.fillna('[]')
#        result.to_csv(output[0], '\t', header=False, index=False)

