import numpy as np
import pandas as pd
rule all:
    input: 'final/gene_list.txt'

rule get_outliers:
    input: 'output'
    output: 'final/outliers.bed'
    run:
        df = pd.read_csv(input[0], '\t')
        condition = (
            (df.Dstat > 0)
        )
        print(df.shape)
        df = df[condition]
        print(df.shape)
        df = df.sort_values(by="score", ascending=False).head(50)
        df.to_csv(output[0], '\t', header=False, index=False)

rule intersect_outliers_and_genes:
    input:
        'final/outliers.bed',
        '/space/s1/joana/refgenomes/Functional_annotation/bedmap/genes.bed'
    output:
        temp('final/intersection.bed')
    shell:
        'bedtools intersect -a {input[0]} -b {input[1]} -loj > {output}'

rule simplify:
    input: 'final/intersection.bed',
    output: temp('final/intersection.bed2'),
    shell: ' cat {input}| cut -f 1-11,22 > {output}'


rule filter_mapping:
    input: 'final/intersection.bed2',
    output: temp('final/intersection.bed3'),
    run:
        in_df = pd.read_csv(input[0], '\t', header=None)
        out_df = in_df.groupby([0,1,2])[11].apply(lambda x: str(list(np.unique(list(x)))))
        out_df.to_csv(output[0], '\t', header=False)


rule merge_toFinal:
    input:
        'final/outliers.bed',
        'final/intersection.bed3'
    output:
        'final/outlier_genes.bed'
    run:
        df1 = pd.read_csv(input[0], '\t', header=None)
        df2 = pd.read_csv(input[1], ',', header=None)
        result = pd.merge(df1,df2, how='inner', on=[0,1,2])
        result['3_y'] = result['3_y'].str.extractall('\((.*?)\)').groupby(level=0)[0].apply(lambda x: str(list(np.unique(list(x)))))
        result = result.fillna('[]')
        result.to_csv(output[0], '\t', header=False, index=False)

rule gene_list:
    input: 'final/outlier_genes.bed'
    output: 'final/gene_list.txt'
    run:
        df = pd.read_csv(input[0], '\t', header=None)
        x = np.concatenate([eval(l) for l in df[11].values.astype(str)])

        indexes = np.unique(x, return_index=True)[1]
        gene_list = [x[index] for index in sorted(indexes)]
        np.savetxt(output[0], gene_list, fmt='%s')

