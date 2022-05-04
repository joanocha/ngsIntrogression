import pandas as pd
import sys
import os
import numpy as np
from scipy.stats import chi2_contingency, chisquare

infile = sys.argv[1]
outfile = sys.argv[2]
#outlier_pos = sys.argv[3]
#outlier_neg = sys.argv[4]
tmp_outfile = outfile + '.unsorted'
data= pd.read_csv(infile, delim_whitespace=True)
simple_data =data.drop_duplicates(['CHR','BLOCKstart', 'BLOCKend'])
# simple_data = data.iloc[::3]
print(simple_data.head())


D = simple_data.Numer / simple_data.Denom

X = (simple_data.Denom + simple_data.Numer)/2.0
Y = (simple_data.Denom - simple_data.Numer)/2.0

global_x = X.sum()
global_y = Y.sum()

global_prop_x = global_x / (global_x + global_y)
global_prop_y = global_y / (global_x + global_y)

X = X.values
Y = Y.values

n_sites = X + Y
expected_X = n_sites * global_prop_x
expected_Y = n_sites * global_prop_y

chi2, pval = zip(*[list(chisquare([X[i], Y[i]], [expected_X[i], expected_Y[i]])) if n_sites[i] >= 5 else [0.0, 1.5] for i in range(len(X))])
score = -np.log10(pval)



#score  = [-np.log10(chi2_contingency(np.array([[X[i], Y[i]], [global_x, global_y]]))[1])
 #         if X[i] > 5 or Y[i] > 5 else np.nan
  #        for i in range(len(X))]
#chi2 = [chi2_contingency(np.array([[X[i], Y[i]], [global_x, global_y]]))[0]
 #       if X[i] > 5 or Y[i] > 5 else np.nan
  #      for i in range(len(X))]

simple_data['X'] = X
simple_data['Y'] = Y
simple_data['Dstat'] = D
simple_data['score'] = score
simple_data['chi2'] = chi2


simple_data.to_csv(tmp_outfile, '\t', index=False)
os.system("cat {} | sort -V -k1,1 -k2,2 > {}".format(tmp_outfile, outfile))
os.system("rm {}".format(tmp_outfile))

print(simple_data[D<0].sort_values('score', ascending=False).head(50))

#outliers_neg = simple_data[D<0].sort_values('score', ascending=False).head(50)
#outliers_neg.to_csv(outlier_neg, '\t', index=False)


#outliers_pos = simple_data[D>0].sort_values('score', ascending=False).head(50)
#outliers_pos.to_csv(outlier_pos, '\t', index=False)
