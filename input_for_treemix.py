import pandas as pd
import numpy as np
import sys


inputs = sys.argv[1:-1]
output = sys.argv[-1]

print(inputs)


def load_input(fname):
    df = pd.read_csv(fname, '\t', header=0,  usecols=["chromo", "position", "knownEM", 'nInd'])
    count_minor = round(2 * df.nInd * df.knownEM).map(int)
    count_major = (2 * df.nInd) - count_minor
    df[fname.split('_')[0]] = count_major.map(str) + ',' + count_minor.map(str)
    df = df.drop(columns=["knownEM", 'nInd'])
    return df


df = load_input(inputs[0])


print(df)

for fname in inputs[1:]:
    df2 = load_input(fname)
    df = pd.merge(df, df2, how='inner', on=['chromo', 'position'])
print(df)
df = df.drop(columns=["chromo", 'position'])
df.to_csv(output, sep=' ', index=False)
