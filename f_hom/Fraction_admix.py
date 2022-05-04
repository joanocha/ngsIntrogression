import pandas as pd
import numpy as np
import sys


numer1 = sys.argv[1]
numer2 = sys.argv[2]
fhom = sys.argv[3]


df = pd.read_csv(numer1, '\t', header=0,  usecols=['CHR','BLOCKstart','BLOCKend','Numer_P1P2P3', 'numSites'])
df2 = pd.read_csv(numer2, '\t', header=0,  usecols=['CHR','BLOCKstart','BLOCKend','Numer_P1P3P3', 'numSites'])
merged_df = pd.merge(df, df2, how='inner', on=['CHR', 'BLOCKstart', 'BLOCKend'])
print(merged_df)

merged_df['f_hom'] = merged_df['Numer_P1P2P3'] / merged_df['Numer_P1P3P3']

print(merged_df.head())

merged_df.to_csv(fhom, sep='\t', index=False)
