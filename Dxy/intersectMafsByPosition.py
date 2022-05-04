import pandas as pd
import numpy as np
import sys


inputs = sys.argv[1:]
print(inputs)

index = pd.read_csv(inputs[0], '\t', index_col=[0, 1]).index

for fname in inputs[1:]:
    index = index.intersection(pd.read_csv(fname, '\t', index_col=[0,1]).index)

pd.DataFrame(index.values.tolist()).to_csv("intersection.txt", sep='\t', index=False)

#cat intersection.txt | sort -V -k1,1 -k2,2 > intersection.sorted.bed 
