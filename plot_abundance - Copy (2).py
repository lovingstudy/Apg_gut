# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 20:20:42 2022

@author: xingjin1
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import seaborn as sns
from scipy.stats import ranksums, fisher_exact

filenames = ['Abundance_sum_D003924.csv', 'Abundance_sum_D006262.csv']
pheno_list = ['T2D', 'Healthy']
result = []
for fn, pheno in zip(filenames, pheno_list):
    df = pd.read_csv(fn, index_col=0)
    df['Phenotype'] = pheno
    result.append(df)
result = pd.concat(result)
#result = result[result['host age in years'] >= 40]

#Clrs = {'T2D': ['royalblue', 'lightskyblue'], 'Healthy': ['darkorange', 'wheat']}
Clrs = {'T2D': ['cornflowerblue', 'lightskyblue'], 'Healthy': ['salmon', 'wheat']}
for pheno in pheno_list:
    FIG = plt.figure(figsize=(2.5, 2.5), dpi=300)
    labels = ['Mah+', 'Mah-']
    sizes = [result[(result['Phenotype'] == pheno) & (result['Found'] == B)].shape[0] \
             for B in [True, False]]
    explode = (0.1, 0)
    plt.pie(sizes, explode=explode, labels=labels, colors=Clrs[pheno], \
            autopct='%1.1f%%', startangle=270)
    plt.title(pheno + ' (' + str(sum(sizes)) + ')')
    plt.tight_layout()
    #FIG.savefig('Pie_%s.pdf'%pheno, transparent=True)
    FIG.savefig('Pie_%s.png'%pheno, transparent=True)

k = result[(result['Phenotype'] == 'T2D') & (result['Found'] == True)].shape[0]
K = result[result['Phenotype'] == 'T2D'].shape[0]
n = result[result['Found'] == True].shape[0]
N = result.shape[0]
print(fisher_exact([[k, K], [n, N]]))

res_found = result[result['Found'] == True]
FIG, axes = plt.subplots(len(pheno_list), 1, figsize=(8, 3 * len(pheno_list)), sharex=True)
for axi, pheno in enumerate(pheno_list):
    axes[axi].hist(res_found.loc[res_found['Phenotype'] == pheno, 'AbundanceSum'], \
                   color=Clrs[pheno][0], bins=20, label=pheno, log=False)
    axes[axi].legend()
    axes[axi].set_ylabel('Number of Samples')
axes[len(pheno_list) - 1].set_xlabel('Abundance %')
FIG.savefig('Histograms_%s.pdf'%('_'.join(pheno_list)), transparent=True)
'''
FIG = plt.figure(figsize=(8, 4))
for pheno in pheno_list:
    sns.histplot(res_found.loc[res_found['Phenotype'] == pheno, 'AbundanceSum'], \
                 stat='probability', bins=np.arange(0, 18, 0.5), element='step', \
                 fill=False, color=Clrs[pheno][0], label=pheno)
plt.xlabel('Abundance %')
plt.ylabel('Fraction')
plt.legend()
FIG.savefig('Histogram1_%s.pdf'%('_'.join(pheno_list)), transparent=True)
'''
FIG, axes = plt.subplots(2, 1, figsize=(6, 5), dpi=300, \
                         gridspec_kw={'height_ratios': [1, 3]})
for pheno in pheno_list:
    sns.histplot(res_found.loc[res_found['Phenotype'] == pheno, 'AbundanceSum'], \
                 stat='probability', bins=np.arange(0, 18, 0.5), element='step', \
                 fill=False, color=Clrs[pheno][0], label=pheno, ax=axes[0])
axes[0].set_ylim(0.1, 0.55)
axes[0].set_xlim(0, 18)
axes[0].set_ylabel(None)
axes[0].set_xlabel(None)
axes[0].spines.bottom.set_visible(False)
axes[0].set_xticks([])
axes[0].set_yticks(np.arange(0.15, 0.5, 0.2))
axes[0].legend()
for pheno in pheno_list:
    sns.histplot(res_found.loc[res_found['Phenotype'] == pheno, 'AbundanceSum'], \
                 stat='probability', bins=np.arange(0, 18, 0.5), element='step', \
                 fill=False, color=Clrs[pheno][0], label=pheno, ax=axes[1])
axes[1].set_ylim(0, 0.1)
axes[1].set_xlim(0, 18)
axes[1].set_ylabel(None)
axes[1].spines.top.set_visible(False)
axes[1].set_xlabel('Abundance %')
plt.tight_layout()
#FIG.savefig('Histogram1_%s.pdf'%('_'.join(pheno_list)), transparent=True)
FIG.savefig('Histogram1_%s.png'%('_'.join(pheno_list)), transparent=True)


A = res_found.loc[res_found['Phenotype'] == 'T2D', 'AbundanceSum']
B = res_found.loc[res_found['Phenotype'] == 'Healthy', 'AbundanceSum']
print(ranksums(A, B))

'''
A = res_age_c.loc[res_age_c['Phenotype'] == 'Healthy', 'AbundanceSum']
B = res_age_c.loc[res_age_c['Phenotype'] == 'T2D', 'AbundanceSum']
print(ranksums(A, B))
sns.kdeplot(A, label='Healthy')
sns.kdeplot(B, label='T2D')
plt.legend()
plt.xlim(0, 22)
plt.show()
sns.boxplot(x='Phenotype', y='AbundanceSum', data=res_age_c)
plt.show()
'''