import pandas as pd
import numpy as np
import sys


numer1 = sys.argv[1]
numer2 = sys.argv[2]
numer3 = sys.argv[3]
fstat = sys.argv[4]


df = pd.read_csv(numer1, '\t', header=0,  usecols=['CHR','BLOCKstart','BLOCKend','Numer_P1P2P3', 'numSites'])
df2 = pd.read_csv(numer2, '\t', header=0,  usecols=['CHR','BLOCKstart','BLOCKend','Numer_P1P2P2', 'numSites'])
df3 = pd.read_csv(numer3, '\t', header=0,  usecols=['CHR','BLOCKstart','BLOCKend','Numer_P1P3P3', 'numSites'])

merged_df1 = pd.merge(df, df2, how='inner', on=['CHR', 'BLOCKstart', 'BLOCKend'])
merged_df = pd.merge(merged_df1, df3, how='inner', on=['CHR', 'BLOCKstart', 'BLOCKend'])

print(merged_df)

merged_df['Numer_D'] = merged_df[['Numer_P1P2P2', 'Numer_P1P3P3']].max(axis=1)
merged_df['f_hom'] = merged_df['Numer_P1P2P3'] / merged_df['Numer_P1P3P3']
merged_df['f_d'] = merged_df['Numer_P1P2P3'] / merged_df['Numer_D']
merged_df = merged_df.fillna()

print(merged_df.head())

merged_df.to_csv(fstat, sep='\t', index=False)
