# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 19:10:19 2022

@author: xingjin1
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

homo = pd.read_csv('../BLAST/Mah1_blast_ncbi_nr_5k_chkd.csv')
homo = homo[homo['Per. ident'] >= 70]
bcts = set(homo['Scientific Name'])

pheno_id, abd_dir, name_col, abd_thr = 'acarbose', './Acarbose/', 'name', 0.01
rm_outliers = False

meta_table = pd.read_csv('all_runs_associated_with_%s.csv'%pheno_id, index_col=0)
meta_table = meta_table[~meta_table.index.duplicated()]
#meta_table = meta_table[meta_table['experiment type'] != 'Amplicon']

result = []
for fn in os.listdir(abd_dir):
    runid = fn.split('_')[0]
    if not runid in meta_table.index: continue
    df = pd.read_csv(abd_dir + fn)
    if pheno_id == 'acarbose':
        if 131567 not in df.tax_id.values: continue
        df['relative_abundance'] = 100 * df['self_count'] / df.loc[df['tax_id'] == 131567, 'total_count'].to_list()[0]
    #if df.shape[0] < 35: continue
    c = set.intersection(bcts, set(df[name_col]))
    if len(c) == 0:
        result.append([runid, df.shape[0], 0, 0])
    else:
        abdsum = df.loc[df[name_col].isin(c), 'relative_abundance'].sum()
        result.append([runid, df.shape[0], len(c), abdsum])
result = pd.DataFrame(data=result, columns=['RunID', 'SpeciesNum', 'OrgMatchNum', 'AbundanceSum'])
result.index = result.iloc[:, 0]
result = result.iloc[:, 1:]
result['Found'] = [v >= abd_thr for v in result['AbundanceSum']]
if rm_outliers:
    result = result[result['AbundanceSum'] < result['AbundanceSum'].quantile(0.95)]

print(result['OrgMatchNum'].value_counts())
plt.hist(result.loc[result['AbundanceSum'] >= abd_thr, 'AbundanceSum'], bins=20)
plt.show()

result = pd.concat([result, meta_table.reindex(result.index)], axis=1)
result.to_csv('Abundance_sum_%s.csv'%pheno_id)
'''
sns.boxplot(y='country', x='AbundanceSum', data=result)
#sns.swarmplot(x='country', y='AbundanceSum', data=result, color='grey')
plt.plot()

plt.figure(); plt.plot(result.SpeciesNum, result.AbundanceSum, 'o'); plt.show()
'''