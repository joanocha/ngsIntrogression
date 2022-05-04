import numpy as np
import pandas as pd

rule all:
    input:
        expand(
            'output/{abbababa}/{percentile}/outlier_genes.txt',
            abbababa=["f_hom"],
            percentile=[95.0, 97.0, 99.0]
        )

rule filter_gff:
    input:
        "/space/s2/joana/refgenomes/Functional_annotation/bedmap/vulpes_lagopus.liftoff.gff"
    output:
        temp("output/genes_filtered.bed")
    run:
        data = pd.read_csv(input[0], '\t', header=None)
        data["gene"] = data[8].str.split(";").str[0].str.split("=").str[1]
        data.to_csv(output[0], '\t', columns=[0, 3,4,"gene"], header=False, index=False)

rule sort_gff:
    input:
        'output/genes_filtered.bed'
    output:
        'output/genes_filtered.sorted.bed'
    shell:
        'bedtools sort -i {input} > {output}'

rule add_genes:
    input:
        config["input_file"],
        'output/genes_filtered.sorted.bed',
    output:
        "output/windows.bed",
    shell:
        'bedtools sort -i {input[0]} | bedtools map -a stdin -b {input[1]} -c 4 -o distinct > {output}'

rule get_outliers:
    input:
        "output/windows.bed",
    output:
        'output/{abbababa}/{percentile}/outliers.bed',
    run:
        abbababa = wildcards["abbababa"]
        percentile = float(wildcards["percentile"])
        df = pd.read_csv(input[0], '\t', header=None, names = ["CHR", "BLOCKstart", "BLOCKend", "Numer_P1P2P3O", "numSites_x", "Numer_P1P3P3O", "numSites_y", "f_hom"])
        threshold = df[abbababa].quantile(percentile / 100.0)
        df = df[df[abbababa] >= threshold]
        df = df.sort_values(by=abbababa, ascending=False)
        df.to_csv(output[0], '\t', index=False)

rule extract_genes:
    input:
        'output/{abbababa}/{percentile}/outliers.bed',
    output:
        'output/{abbababa}/{percentile}/outlier_genes.txt'
    run:
        df = pd.read_csv(input[0], '\t')
        print(df)
        df = df[df.genes != "."]
        df["genes"] =  df.genes.str.split(",")
        genes = pd.unique(np.concatenate(df.genes.values))
        pd.DataFrame(genes).to_csv(output[0], index=False, header=False)

