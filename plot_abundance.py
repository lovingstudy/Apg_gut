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


result = pd.read_csv('Abundance_sum_D006262_idt85.csv', index_col=0)
result['country'] = ['USA' if b else c for b, c in zip(result['country'].str.contains('United States of America'), result['country'])]
result['AgeGroup'] = ['%s~%s'%(int(age-age%10), int(age-age%10+9)) if not np.isnan(age) else np.nan for age in result['host age in years']]
bmi_thrs = [lambda x: (np.isnan(x), np.nan), lambda x: (x < 18.5, 'Underweight'), \
            lambda x: (18.5 <= x < 25, 'Normal'), lambda x: (25 <= x < 30, 'Overweight'), \
            lambda x: (x >= 30, 'Obesity')]
tmp = [[f(i) for f in bmi_thrs] for i in result['BMI']]
result['BMIGroup'] = [[l for b, l in t if b][0] for t in tmp]

FIGSIZE = (4, 5)
pos_color, neg_color = ['cornflowerblue', 'silver']
p_thrs = [lambda x: (x < 0.001, '***'), lambda x: (0.001 <= x < 0.01, '**'), \
          lambda x: (0.01 <= x < 0.05, '*'), lambda x: (x >= 0.05, '')]
bmi_compare_health = True

cnt = result.country.value_counts()
res_country = result[result['country'].isin(cnt[cnt >= 100].index)]
prv = pd.DataFrame(res_country.groupby(by='country').mean()['Found'] * 100)
prv = prv.rename(columns={'Found': 'apg+'})
prv.loc['All'] = [result.loc[~result['country'].isna(), 'Found'].mean() * 100]
prv['apg-'] = 100 - prv['apg+']
prv = pd.concat([prv[prv.index != 'All'].sort_values(by='apg+'), prv.loc[['All']]])
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
plt.barh(prv.index, prv['apg+'], label='apg+', color=pos_color)
plt.barh(prv.index, prv['apg-'], label='apg-', color=neg_color, left=prv['apg+'])
for g in prv.index:
    p = prv.loc[g, 'P']
    if p >= 1: continue
    else:
        astr = [t(p) for t in p_thrs]
        astr = [txt for b, txt in astr if b][0]
        plt.text(101, g, astr, va='center')
        plt.text(98, g, 'p = %.2E'%p, va='center', ha='right')
plt.ylim(-1, prv.shape[0] + 1)
#plt.xlim(0, 130)
plt.legend(ncol=2, loc='upper left')
plt.xlabel('Prevalence %')
plt.title('Country')
plt.tight_layout()
FIG.savefig('Bars_Prevalence_countries_idt85.pdf', transparent=True)


cnt = result['AgeGroup'].value_counts()
res_age = result[result['AgeGroup'].isin(cnt[cnt >= 50].index)]
res_age = res_age.sort_values(by='AgeGroup', ascending=False)
prv = pd.DataFrame(res_age.groupby(by='AgeGroup').mean()['Found'] * 100)
prv = prv.rename(columns={'Found': 'apg+'})
prv.loc['All'] = [result.loc[~result['host age in years'].isna(), 'Found'].mean() * 100]
prv['apg-'] = 100 - prv['apg+']
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
plt.barh(prv.index, prv['apg+'], label='apg+', color=pos_color)
plt.barh(prv.index, prv['apg-'], label='apg-', color=neg_color, left=prv['apg+'])
for g in prv.index:
    p = prv.loc[g, 'P']
    if p >= 1: continue
    else:
        astr = [t(p) for t in p_thrs]
        astr = [txt for b, txt in astr if b][0]
        plt.text(101, g, astr, va='center')
        plt.text(98, g, 'p = %.2E'%p, va='center', ha='right')
plt.ylim(-1, prv.shape[0] + 1)
#plt.xlim(0, 130)
plt.legend(ncol=2, loc='upper left')
plt.xlabel('Prevalence %')
plt.title('Age')
plt.tight_layout()
FIG.savefig('Bars_Prevalence_age_idt85.pdf', transparent=True)
matplotlib.rcParams['font.size'] = 10

cnt = result['BMIGroup'].value_counts()
res_bmi = result[(result['host age in years'] >= 20) & \
                 (result['BMIGroup'].isin(cnt[cnt >= 20].index)) & \
                 (~result['BMI'].isna())]
res_bmi = res_bmi.sort_values(by='BMIGroup', ascending=False)
res_bmi_h = res_bmi[res_bmi['BMIGroup'] == 'Normal']
prv = pd.DataFrame(res_bmi.groupby(by='BMIGroup').mean()['Found'] * 100)
prv = prv.rename(columns={'Found': 'apg+'})
prv.loc['All'] = [res_bmi['Found'].mean() * 100]
prv['apg-'] = 100 - prv['apg+']
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
plt.barh(prv.index, prv['apg+'], label='apg+', color=pos_color)
plt.barh(prv.index, prv['apg-'], label='apg-', color=neg_color, left=prv['apg+'])
for g in prv.index:
    p = prv.loc[g, 'P']
    if p >= 1: continue
    else:
        astr = [t(p) for t in p_thrs]
        astr = [txt for b, txt in astr if b][0]
        plt.text(101, g, astr, va='center')
        plt.text(98, g, 'p = %.2E'%p, va='center', ha='right')
        '''
        tmp = 'Normal' if bmi_compare_health else 'All'
        yt = [i for i in prv.index if i.startswith(tmp)][0]
        mid = (prv.index.get_loc(yt) + prv.index.get_loc(g)) / 2
        plt.plot([101, 105], [yt, yt], color='black', lw=0.5)
        plt.plot([105, 105], [yt, g], color='black', lw=0.5)
        plt.plot([101, 105], [g, g], color='black', lw=0.5)
        plt.text(106, mid, astr, va='center')
        '''
plt.ylim(-1, prv.shape[0])
#plt.xlim(0, 130)
plt.legend(ncol=2, loc='upper left')
plt.xlabel('Prevalence %')
plt.title('Body Mass')
plt.tight_layout()
FIG.savefig('Bars_Prevalence_BMI_idt85.pdf', transparent=True)


res_gender = result[(result['host age in years'] >= 18) & (~result['gender'].isna())]
k = res_gender[(res_gender['gender'] == 'Male') & (res_gender['Found'] == True)].shape[0]
K = res_gender[res_gender['gender'] == 'Male'].shape[0]
n = res_gender['Found'].sum()
od, p = fisher_exact([[k, K], [n, res_gender.shape[0]]])
print('Gender:', '%.2f'%od, '%.2E'%p)
prv = pd.DataFrame(res_gender.groupby(by='gender').mean()['Found'] * 100)
prv = prv.rename(columns={'Found': 'apg+'})
prv['apg-'] = 100 - prv['apg+']
new_index = []
for ag in prv.index:
    sub = result[result['gender'] == ag]
    new_index.append(ag + ' (%s)'%sub.shape[0])
prv.index = new_index
print(prv)

FIG = plt.figure(figsize=(3, 5), dpi=300)
plt.bar(prv.index, prv['apg+'], label='apg+', color=pos_color)
plt.bar(prv.index, prv['apg-'], label='apg-', color=neg_color, bottom=prv['apg+'])
plt.ylim(0, 130)
plt.legend(ncol=2, loc='upper left')
plt.ylabel('Prevalence %')

plt.plot([0, 0], [101, 110], color='black', lw=0.5)
plt.plot([0, 1], [110, 110], color='black', lw=0.5)
plt.plot([1, 1], [101, 110], color='black', lw=0.5)
plt.text(x=0.5, y=111, s='p = %.2f'%p, ha='center', va='bottom')

plt.title('Gender')
plt.tight_layout()
FIG.savefig('Bars_Prevalence_gender_idt85.pdf', transparent=True)


matplotlib.rcParams['font.size'] = 14
tmp = res_country[res_country['Found'] == True]
tmp['country'] = [i + ' (%s)'%tmp['country'].value_counts().loc[i] for i in tmp['country']]
FIG = plt.figure(dpi=300)
sns.boxenplot(x='AbundanceSum', y='country', data=tmp)
#plt.xlim(-1, 65)
plt.xlabel('Abundance %')
plt.ylabel(None)
plt.title('Country')
plt.tight_layout()
FIG.savefig('Boxplot_abundance_countries_idt85.pdf', transparent=True)

FIG = plt.figure(dpi=300)
tmp = res_age[res_age['Found'] == True]
tmp['AgeGroup'] = [i + ' (%s)'%tmp['AgeGroup'].value_counts().loc[i] for i in tmp['AgeGroup']]
sns.boxenplot(x='AbundanceSum', y='AgeGroup', data=tmp)
#plt.xlim(-1, 20)
plt.xlabel('Abundance %')
plt.ylabel(None)
plt.title('Age')
plt.tight_layout()
FIG.savefig('Boxplot_abundance_age_idt85.pdf', transparent=True)

FIG = plt.figure(dpi=300)
tmp = res_bmi[res_bmi['Found'] == True]
tmp['BMIGroup'] = [i + ' (%s)'%tmp['BMIGroup'].value_counts().loc[i] for i in tmp['BMIGroup']]
sns.boxenplot(x='AbundanceSum', y='BMIGroup', data=tmp)
#plt.xlim(-1, 20)
plt.xlabel('Abundance %')
plt.ylabel(None)
plt.title('Body Mass')
plt.tight_layout()
FIG.savefig('Boxplot_abundance_BMI_idt85.pdf', transparent=True)

FIG = plt.figure(dpi=300)
tmp = res_gender[res_gender['Found'] == True]
tmp['gender'] = [i + ' (%s)'%tmp['gender'].value_counts().loc[i] for i in tmp['gender']]
sns.boxenplot(x='AbundanceSum', y='gender', data=tmp)
#plt.xlim(-1, 20)
plt.xlabel('Abundance %')
plt.ylabel(None)
plt.title('Gender')
plt.tight_layout()
FIG.savefig('Boxplot_abundance_gender_idt85.pdf', transparent=True)




