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
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False
import seaborn as sns


acar_neg = pd.read_csv('Abundance_sum_D006262.csv', index_col=0)
acar_pos = pd.read_csv('apg_acar_1.csv')
acar_pos['Found'] = [True if i > 0.0001 else False for i in acar_pos['apg_abundance']]

prev_neg = acar_neg['Found'].mean() * 100
prev_pos = acar_pos['Found'].mean() * 100
FIG = plt.figure(figsize=(2.5, 3), dpi=300)
plt.bar(['Healthy\n(%s)'%acar_neg.shape[0], 'T2DM\n(%s)'%acar_pos.shape[0]], \
        [prev_neg, prev_pos], color=['cornflowerblue', 'salmon'], label='Apg+')
'''
plt.bar(['T2DM\n(%s)'%acar_pos.shape[0]], [prev_pos], color='salmon', label='Apg+')
plt.bar(['Healthy\n(%s)'%acar_neg.shape[0], 'T2DM\n(%s)'%acar_pos.shape[0]], \
        [100 - prev_neg, 100 - prev_pos], color='silver', label='Apg-', \
        bottom=[prev_neg, prev_pos])
'''
#plt.ylabel('Prevalence %')
#plt.xlim(-0.5, 3)
#plt.legend(loc='upper right')
plt.tight_layout()
FIG.savefig('Bars_prevl_healthy_acarbose.png', transparent=True)

neg = acar_neg['AbundanceSum']
pos = acar_pos['apg_abundance'] * 100
FIG, axes = plt.subplots(2, 1, figsize=(5, 4), dpi=300, \
                         gridspec_kw={'height_ratios': [1, 5]})
sns.histplot(pos, stat='probability', bins=np.arange(0.0, 70, 0.5), element='step', \
             color='salmon', fill=False, ax=axes[0])
sns.histplot(neg, stat='probability', bins=np.arange(0.0, 20, 0.5), element='step', \
             color='cornflowerblue', fill=False, ax=axes[0])
axes[0].set_ylim(0.55, 0.75)
axes[0].set_ylabel(None)
axes[0].set_xlabel(None)
axes[0].spines.bottom.set_visible(False)
axes[0].set_xticks([])
sns.histplot(pos, stat='probability', bins=np.arange(0.0, 70, 0.5), element='step', \
             color='salmon', fill=False, ax=axes[1])
sns.histplot(neg, stat='probability', bins=np.arange(0.0, 20, 0.5), element='step', \
             color='cornflowerblue', fill=False, ax=axes[1])
axes[1].set_ylim(0, 0.15)
plt.xlabel('Abundance %')
plt.tight_layout()
FIG.savefig('Histogram_Healthy_Acarbose.png', transparent=True)

FIG, axes = plt.subplots(2, 2, figsize=(5, 4), dpi=300, \
                         gridspec_kw={'height_ratios': [1, 5], 'width_ratios': [3, 1]})
sns.histplot(pos, stat='probability', bins=np.arange(0.0, 70, 0.5), element='step', \
             color='salmon', fill=False, ax=axes[0][0], label='T2DM')
sns.histplot(neg, stat='probability', bins=np.arange(0.0, 15, 0.5), element='step', \
             color='cornflowerblue', fill=False, ax=axes[0][0], label='Healthy')
axes[0][0].set_ylim(0.55, 0.75)
axes[0][0].set_xlim(0, 15)
axes[0][0].set_ylabel(None)
axes[0][0].set_xlabel(None)
axes[0][0].spines.bottom.set_visible(False)
axes[0][0].set_xticks([])
axes[0][0].legend(loc='upper center')
sns.histplot(pos, stat='probability', bins=np.arange(0.0, 70, 0.5), element='step', \
             color='salmon', fill=False, ax=axes[1][0])
sns.histplot(neg, stat='probability', bins=np.arange(0.0, 15, 0.5), element='step', \
             color='cornflowerblue', fill=False, ax=axes[1][0])
axes[1][0].set_ylim(0, 0.15)
axes[1][0].set_xlim(0, 15)
axes[1][0].set_xlabel(None)
axes[1][0].set_ylabel(None)
sns.histplot(pos, stat='probability', bins=np.arange(0.0, 70, 2), element='step', \
             color='salmon', fill=False, ax=axes[1][1])
sns.histplot(neg, stat='probability', bins=np.arange(0.0, 15, 2), element='step', \
             color='cornflowerblue', fill=False, ax=axes[1][1])
axes[1][1].set_ylim(0, 0.15)
axes[1][1].set_xlim(15, 72)
axes[1][1].set_xlabel(None)
axes[1][1].set_ylabel(None)
axes[1][1].set_yticks([])
axes[1][1].spines.left.set_visible(False)

axes[0][1].set_axis_off()
plt.tight_layout()
FIG.savefig('Histogram_Healthy_Acarbose-1.png', transparent=True)
