# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 20:06:48 2022

@author: xingjin1
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['font.family'] = 'Arial'
import seaborn as sns


acar_neg = pd.read_csv('Abundance_sum_D006262_idt85.csv', index_col=0)
acar_pos = pd.read_csv('apg_acar_1.csv')
acar_pos['Found'] = [True if i > 0.0001 else False for i in acar_pos['apg_abundance']]

prev_neg = acar_neg['Found'].mean() * 100
prev_pos = acar_pos['Found'].mean() * 100
FIG = plt.figure(figsize=(2.5, 3), dpi=300)
plt.bar(['Healthy\n(%s)'%acar_neg.shape[0], 'T2DM\n(%s)'%acar_pos.shape[0]], \
        [prev_neg, prev_pos], color=['cornflowerblue', 'salmon'], label='Apg+')
plt.text(0, prev_neg + 1, '%.1f'%prev_neg + '%', color='black', ha='center')
plt.text(1, prev_pos + 1, '%.1f'%prev_pos + '%', color='black', ha='center')
'''
plt.bar(['T2DM\n(%s)'%acar_pos.shape[0]], [prev_pos], color='salmon', label='Apg+')
plt.bar(['Healthy\n(%s)'%acar_neg.shape[0], 'T2DM\n(%s)'%acar_pos.shape[0]], \
        [100 - prev_neg, 100 - prev_pos], color='silver', label='Apg-', \
        bottom=[prev_neg, prev_pos])
'''
plt.ylabel('apg+ samples / Total samples')
#plt.xlim(-0.5, 3)
#plt.legend(loc='upper right')
plt.tight_layout()
FIG.savefig('Bars_prevl_healthy_acarbose_idt85.pdf', transparent=True)

neg = acar_neg['AbundanceSum']
pos = acar_pos['apg_abundance'] * 100

FIG, axes = plt.subplots(2, 1, figsize=(4, 3.5), dpi=300, \
                         gridspec_kw={'height_ratios': [1, 5]})
sns.histplot(pos, stat='probability', bins=np.arange(0.0, 10, 0.1), element='step', \
             color='salmon', fill=False, ax=axes[0], label='T2DM')
sns.histplot(neg, stat='probability', bins=np.arange(0.0, 10, 0.1), element='step', \
             color='cornflowerblue', fill=False, ax=axes[0], label='Healthy')
axes[0].set_ylim(0.9, 1)
axes[0].set_xlim(0, 10)
axes[0].set_ylabel(None)
axes[0].set_xlabel(None)
axes[0].spines["bottom"].set_visible(False)
axes[0].set_xticks([])
axes[0].legend(loc='upper center')
sns.histplot(pos, stat='probability', bins=np.arange(0.0, 10, 0.1), element='step', \
             color='salmon', fill=False, ax=axes[1])
sns.histplot(neg, stat='probability', bins=np.arange(0.0, 10, 0.1), element='step', \
             color='cornflowerblue', fill=False, ax=axes[1], label='Healthy')
axes[1].set_ylim(0, 0.25)
axes[1].set_xlim(0, 10)
axes[1].set_ylabel('apg+ samples / Total samples')
axes[1].set_xlabel('The abundance of apg+ bacteria (%)')
#plt.xlabel('Abundance %')
plt.tight_layout()
FIG.savefig('Histogram_Healthy_Acarbose_idt85.pdf', transparent=True)
