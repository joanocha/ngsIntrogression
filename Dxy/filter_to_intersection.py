import pandas as pd
import numpy as np
import sys
#make sure intersection.sorted.txt with  intersection.txt | sort -V -k1,1 -k2,2 > intersection.sorted.txt
input_path = sys.argv[1]
intersection_path = sys.argv[2]
output_path = sys.argv[1] + ".filt_to_intersection.maf.gz"

index = pd.read_csv(intersection_path, '\t', index_col=[0, 1]).index.set_names(["chromo", "position"])
df = pd.read_csv(input_path, '\t', index_col=[0, 1])
print(df.shape)
df = df.loc[index]
print(df.shape)
df.reset_index().to_csv(output_path, "\t", index=False)
