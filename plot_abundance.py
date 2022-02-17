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
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False
import seaborn as sns
from scipy.stats import ranksums, fisher_exact


result = pd.read_csv('Abundance_sum_D006262.csv', index_col=0)
result['country'] = ['USA' if b else c for b, c in zip(result['country'].str.contains('United States of America'), result['country'])]
result['AgeGroup'] = ['%s~%s'%(int(age-age%10), int(age-age%10+9)) if not np.isnan(age) else np.nan for age in result['host age in years']]
bmi_thrs = [lambda x: (np.isnan(x), np.nan), lambda x: (x < 18.5, 'Underweight'), \
            lambda x: (18.5 <= x < 25, 'Normal'), lambda x: (25 <= x < 30, 'Overweight'), \
            lambda x: (x >= 30, 'Obesity')]
#result['BMIGroup'] = ['[%s, %s)'%(int(bmi-bmi%5), int(bmi-bmi%5+5)) if not np.isnan(bmi) else np.nan for bmi in result['BMI']]
tmp = [[f(i) for f in bmi_thrs] for i in result['BMI']]
result['BMIGroup'] = [[l for b, l in t if b][0] for t in tmp]

FIGSIZE = (4, 5)
pos_color, neg_color = ['cornflowerblue', 'silver']
p_thrs = [lambda x: (x < 0.001, '***'), lambda x: (0.001 <= x < 0.01, '**'), lambda x: (x >= 0.01, '*')]
bmi_compare_health = True

FIG = plt.figure(figsize=(2.5, 2.5), dpi=300)
labels = ['Apg+', 'Apg-']
sizes = [result[result['Found'] == B].shape[0] for B in [True, False]]
explode = (0.1, 0)
plt.pie(sizes, explode=explode, labels=labels, colors=[pos_color, neg_color], \
        autopct='%1.1f%%', startangle=0)
#plt.title('Prevalence' + ' (' + str(sum(sizes)) + ' in total)')
plt.tight_layout()
#FIG.savefig('Pie_%s.pdf'%pheno, transparent=True)
FIG.savefig('Pie_Healthy.png', transparent=True)


res_found = result[result['Found'] == True]
FIG, axes = plt.subplots(2, 1, figsize=(4, 4), dpi=300, \
                         gridspec_kw={'height_ratios': [1, 5]})
sns.histplot(res_found['AbundanceSum'], bins=np.arange(0, 14, 0.5), color=pos_color, ax=axes[0])
axes[0].set_ylim(1500, 1600)
axes[0].set_xlim(0, 14)
axes[0].set_ylabel(None)
axes[0].set_xlabel(None)
axes[0].spines.bottom.set_visible(False)
axes[0].set_xticks([])
axes[0].set_title('Total (%s)'%result.shape[0])

sns.histplot(res_found['AbundanceSum'], bins=np.arange(0, 14, 0.5), color=pos_color, ax=axes[1])
axes[1].set_ylim(0, 405)
axes[1].set_xlim(0, 14)
axes[1].set_ylabel('Sample Count')
axes[1].spines.top.set_visible(False)
axes[1].set_xlabel('Abundance %')
plt.tight_layout()
FIG.savefig('Histogram_Healthy.png', transparent=True)


cnt = result.country.value_counts()
res_country = result[result['country'].isin(cnt[cnt >= 100].index)]
prv = pd.DataFrame(res_country.groupby(by='country').mean()['Found'] * 100)
prv = prv.rename(columns={'Found': 'Apg+'})
prv.loc['All'] = [result.loc[~result['country'].isna(), 'Found'].mean() * 100]
prv['Apg-'] = 100 - prv['Apg+']
prv = pd.concat([prv[prv.index != 'All'].sort_values(by='Apg+'), prv.loc[['All']]])
new_index = []
for country in prv.index:
    sub = result[result['country'] == country]
    if country == 'All': sub = result[~result['country'].isna()]
    new_index.append(country + ' (%s)'%sub.shape[0])
    k = sub['Found'].sum()
    K = result.loc[~result['country'].isna(), 'Found'].sum()
    od, p = fisher_exact([[k, K], [sub.shape[0], result[~result['country'].isna()].shape[0]]])
    prv.loc[country, 'P'] = p
prv.index = new_index
print(prv)

FIG = plt.figure(figsize=FIGSIZE, dpi=300)
plt.barh(prv.index, prv['Apg+'], label='Apg+', color=pos_color)
plt.barh(prv.index, prv['Apg-'], label='Apg-', color=neg_color, left=prv['Apg+'])
for g in prv.index:
    p = prv.loc[g, 'P']
    if p >= 0.05: continue
    else:
        astr = [t(p) for t in p_thrs]
        astr = [txt for b, txt in astr if b][0]
        plt.text(101, g, astr, va='center')
plt.ylim(-1, prv.shape[0] + 1)
#plt.xlim(0, 130)
plt.legend(ncol=2, loc='upper left')
plt.xlabel('Prevalence %')
plt.title('Country')
plt.tight_layout()
FIG.savefig('Bars_Prevalence_countries.png', transparent=True)


cnt = result['AgeGroup'].value_counts()
res_age = result[result['AgeGroup'].isin(cnt[cnt >= 50].index)]
res_age = res_age.sort_values(by='AgeGroup', ascending=False)
prv = pd.DataFrame(res_age.groupby(by='AgeGroup').mean()['Found'] * 100)
prv = prv.rename(columns={'Found': 'Apg+'})
prv.loc['All'] = [result.loc[~result['host age in years'].isna(), 'Found'].mean() * 100]
prv['Apg-'] = 100 - prv['Apg+']
prv = pd.concat([prv[prv.index != 'All'].sort_index(ascending=False), prv.loc[['All']]])
new_index = []
for ag in prv.index:
    sub = result[result['AgeGroup'] == ag]
    if ag == 'All': sub = result[~result['host age in years'].isna()]
    new_index.append(ag + ' (%s)'%sub.shape[0])
    k = sub['Found'].sum()
    K = result.loc[~result['host age in years'].isna(), 'Found'].sum()
    od, p = fisher_exact([[k, K], [sub.shape[0], result[~result['host age in years'].isna()].shape[0]]])
    prv.loc[ag, 'P'] = p
prv.index = new_index
print(prv)
matplotlib.rcParams['font.size'] = 14
FIG = plt.figure(figsize=(5, 5), dpi=300)
plt.barh(prv.index, prv['Apg+'], label='Apg+', color=pos_color)
plt.barh(prv.index, prv['Apg-'], label='Apg-', color=neg_color, left=prv['Apg+'])
for g in prv.index:
    p = prv.loc[g, 'P']
    if p >= 0.05: continue
    else:
        astr = [t(p) for t in p_thrs]
        astr = [txt for b, txt in astr if b][0]
        plt.text(101, g, astr, va='center')
plt.ylim(-1, prv.shape[0] + 1)
#plt.xlim(0, 130)
plt.legend(ncol=2, loc='upper left')
plt.xlabel('Prevalence %')
plt.title('Age')
plt.tight_layout()
FIG.savefig('Bars_Prevalence_age.png', transparent=True)
matplotlib.rcParams['font.size'] = 10

cnt = result['BMIGroup'].value_counts()
res_bmi = result[(result['host age in years'] >= 20) & \
                 (result['BMIGroup'].isin(cnt[cnt >= 20].index)) & \
                 (~result['BMI'].isna())]
res_bmi = res_bmi.sort_values(by='BMIGroup', ascending=False)
res_bmi_h = res_bmi[res_bmi['BMIGroup'] == 'Normal']
prv = pd.DataFrame(res_bmi.groupby(by='BMIGroup').mean()['Found'] * 100)
prv = prv.rename(columns={'Found': 'Apg+'})
prv.loc['All'] = [res_bmi['Found'].mean() * 100]
prv['Apg-'] = 100 - prv['Apg+']
order = ['Obesity', 'Overweight', 'Normal', 'Underweight', 'All']
prv = pd.concat([prv[prv.index == i] for i in order])
new_index = []
for ag in prv.index:
    sub = res_bmi[res_bmi['BMIGroup'] == ag]
    if ag == 'All': sub = res_bmi.copy()
    new_index.append(ag + ' (%s)'%sub.shape[0])
    k = sub['Found'].sum()
    K = res_bmi_h['Found'].sum() if bmi_compare_health else res_bmi['Found'].sum()
    N = res_bmi_h.shape[0] if bmi_compare_health else res_bmi.shape[0]
    od, p = fisher_exact([[k, K], [sub.shape[0], N]])
    prv.loc[ag, 'P'] = p
prv.index = new_index
if bmi_compare_health: prv = prv.loc[~prv.index.str.startswith('All')]
print(prv)

FIG = plt.figure(figsize=FIGSIZE, dpi=300)
plt.barh(prv.index, prv['Apg+'], label='Apg+', color=pos_color)
plt.barh(prv.index, prv['Apg-'], label='Apg-', color=neg_color, left=prv['Apg+'])
for g in prv.index:
    p = prv.loc[g, 'P']
    if p >= 0.05: continue
    else:
        astr = [t(p) for t in p_thrs]
        astr = [txt for b, txt in astr if b][0]
        tmp = 'Normal' if bmi_compare_health else 'All'
        yt = [i for i in prv.index if i.startswith(tmp)][0]
        mid = (prv.index.get_loc(yt) + prv.index.get_loc(g)) / 2
        plt.plot([101, 105], [yt, yt], color='black', lw=0.5)
        plt.plot([105, 105], [yt, g], color='black', lw=0.5)
        plt.plot([101, 105], [g, g], color='black', lw=0.5)
        plt.text(106, mid, astr, va='center')
plt.ylim(-1, prv.shape[0])
#plt.xlim(0, 130)
plt.legend(ncol=2, loc='upper left')
plt.xlabel('Prevalence %')
plt.title('Fatness')
plt.tight_layout()
FIG.savefig('Bars_Prevalence_BMI.png', transparent=True)


res_gender = result[(result['host age in years'] >= 18) & (~result['gender'].isna())]
k = res_gender[(res_gender['gender'] == 'Male') & (res_gender['Found'] == True)].shape[0]
K = res_gender[res_gender['gender'] == 'Male'].shape[0]
n = res_gender['Found'].sum()
od, p = fisher_exact([[k, K], [n, res_gender.shape[0]]])
print('Gender:', '%.2f'%od, '%.2E'%p)
prv = pd.DataFrame(res_gender.groupby(by='gender').mean()['Found'] * 100)
prv = prv.rename(columns={'Found': 'Apg+'})
prv['Apg-'] = 100 - prv['Apg+']
new_index = []
for ag in prv.index:
    sub = result[result['gender'] == ag]
    new_index.append(ag + ' (%s)'%sub.shape[0])
prv.index = new_index
print(prv)

FIG = plt.figure(figsize=(3, 5), dpi=300)
plt.bar(prv.index, prv['Apg+'], label='Apg+', color=pos_color)
plt.bar(prv.index, prv['Apg-'], label='Apg-', color=neg_color, bottom=prv['Apg+'])
plt.ylim(0, 130)
plt.legend(ncol=2, loc='upper left')
plt.ylabel('Prevalence %')
plt.plot([0, 0], [101, 110], color='black', lw=0.5)
plt.plot([0, 1], [110, 110], color='black', lw=0.5)
plt.plot([1, 1], [101, 110], color='black', lw=0.5)
plt.text(x=0.5, y=111, s='**', ha='center', va='center')
plt.title('Gender')
plt.tight_layout()
FIG.savefig('Bars_Prevalence_gender.png', transparent=True)


matplotlib.rcParams['font.size'] = 14
FIG = plt.figure(dpi=300)
sns.boxenplot(x='AbundanceSum', y='country', data=res_country[res_country['Found'] == True])
plt.xlabel('Abundance %')
plt.ylabel(None)
plt.title('Country')
plt.tight_layout()
FIG.savefig('Boxplot_abundance_countries.png', transparent=True)

FIG = plt.figure(dpi=300)
sns.boxenplot(x='AbundanceSum', y='AgeGroup', data=res_age[res_age['Found'] == True])
plt.xlabel('Abundance %')
plt.ylabel(None)
plt.title('Age')
plt.tight_layout()
FIG.savefig('Boxplot_abundance_age.png', transparent=True)

FIG = plt.figure(dpi=300)
sns.boxenplot(x='AbundanceSum', y='BMIGroup', data=res_bmi[res_bmi['Found'] == True])
plt.xlabel('Abundance %')
plt.ylabel(None)
plt.title('Fatness')
plt.tight_layout()
FIG.savefig('Boxplot_abundance_BMI.png', transparent=True)

FIG = plt.figure(dpi=300)
sns.boxenplot(x='AbundanceSum', y='gender', data=res_gender[res_gender['Found'] == True])
plt.xlabel('Abundance %')
plt.ylabel(None)
plt.title('Gender')
plt.tight_layout()
FIG.savefig('Boxplot_abundance_gender.png', transparent=True)




