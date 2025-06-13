# ##### Human Adipose Depot DIA (Figure 8)
# ##### Yue (Winnie) Wen, Alex Zelter, Mike Riffle, Michael J. MacCoss, Nina Isoherranen
# ##### Department of Pharmaceutics, Department of Genome Science, University of Washington-Seattle
# ##### 06/12/2025

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn
from matplotlib.patches import Patch
from matplotlib.ticker import ScalarFormatter

# figure a

cleandata = pd.read_csv("figure-6-7-8-data.csv")
paired_participants = [3501,3502,3503,3504,3505,3506,3507,
                               3508,3509,3510,3512,3514,3515,
                               3516,3517,3518,3519,3520,3521,3522,
                               3523,3524,3525,3527,3528,3529,
                               3531,3532]
cleandata = cleandata[cleandata.id.isin(paired_participants)]

OM_SC_paired_tests = pd.read_csv("figure-8a-data.csv")
protein_top10_OM = OM_SC_paired_tests[OM_SC_paired_tests.Fold_Change<0.5].sort_values(by='Adjusted_P_Value', ascending=True).Protein[:10]
protein_top10_SC = OM_SC_paired_tests[OM_SC_paired_tests.Fold_Change>2].sort_values(by='Adjusted_P_Value', ascending=True).Protein[:10]

cleandata_top10_OM = pd.merge(pd.DataFrame({'protein':protein_top10_OM}), cleandata, left_on='protein', right_on='protein',how='left')
cleandata_top10_SC = pd.merge(pd.DataFrame({'protein':protein_top10_SC}), cleandata, left_on='protein', right_on='protein',how='left')

cleandata_top10_OM_pivot = cleandata_top10_OM.pivot(index=['protein','id'], columns='tissue_type', values='value').reset_index()
cleandata_top10_OM_pivot['value'] = np.log2(cleandata_top10_OM_pivot.SQ/cleandata_top10_OM_pivot.OM)
cleandata_top10_OM_pivot['protein_cleaned'] = cleandata_top10_OM_pivot['protein'].str.split('|').str[2].str.split('_').str[0]

cleandata_top10_SC_pivot = cleandata_top10_SC.pivot(index=['protein','id'], columns='tissue_type', values='value').reset_index()
cleandata_top10_SC_pivot['value'] = np.log2(cleandata_top10_SC_pivot.SQ/cleandata_top10_SC_pivot.OM)
cleandata_top10_SC_pivot['protein_cleaned'] = cleandata_top10_SC_pivot['protein'].str.split('|').str[2].str.split('_').str[0]

order_lst = OM_SC_paired_tests[OM_SC_paired_tests.Protein.isin(cleandata_top10_OM.protein)].\
            sort_values(by='Adjusted_P_Value').\
            Protein.str.split('|').str[2].str.split('_').str[0].tolist() + \
                OM_SC_paired_tests[OM_SC_paired_tests.Protein.isin(cleandata_top10_SC.protein)]. \
            sort_values(by='Adjusted_P_Value'). \
            Protein.str.split('|').str[2].str.split('_').str[0].tolist()
om_col = [153/255,0/255,0/255]
sc_col = [0/255,76/255,153/255]
palette = 10*[om_col] + 10*[sc_col]

adj_p_val_data = OM_SC_paired_tests[OM_SC_paired_tests.Protein.isin(cleandata_top10_OM.protein.tolist()+cleandata_top10_SC.protein.tolist())]
adj_p_val_data['protein_cleaned'] = adj_p_val_data.Protein.str.split('|').str[2].str.split('_').str[0]

figs, axes = plt.subplots(nrows=2, ncols=1, figsize= (12,9), gridspec_kw={'height_ratios': [1, 4]})
sbn.barplot(ax=axes[0], data=adj_p_val_data, x='protein_cleaned', y='Adjusted_P_Value', order=order_lst, palette=palette)
sbn.boxplot(ax=axes[1], data=pd.concat([cleandata_top10_OM_pivot, cleandata_top10_SC_pivot]), x='protein_cleaned', y='value',fill=False, color=[0/255,0/255,0/255],linewidth=2,
            order=order_lst)
sbn.stripplot(ax=axes[1], data=pd.concat([cleandata_top10_OM_pivot, cleandata_top10_SC_pivot]), x='protein_cleaned', y='value', palette=palette,
              s=10,alpha = 0.2, 
              order=order_lst)
sbn.despine(left=True)
axes[1].tick_params(bottom=False)
sbn.despine(top=True, bottom=True)
axes[1].axhline(y=0, color='black')
axes[0].set_title("Top 10 Proteins Enriched in OM and SC", size=20)
axes[0].set_ylabel('Adjusted p value', fontsize=18)
axes[1].set_ylabel('log$_{2}$ fold change (SC/OM)', fontsize=18)
axes[1].set_ylim(-7, 7)
axes[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axes[1].tick_params(axis='x', rotation=90)
axes[0].yaxis.get_offset_text().set_fontsize(18)

for ax in axes:
    ax.set_xticks(ax.get_xticks())         # Needed before rotating xticks
    ax.set_xlabel('', fontsize=18)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)

legend_elements = [
    Patch(facecolor=om_col, label='Enriched in OM'),
    Patch(facecolor=sc_col, label='Enriched in SC')
]

ax.legend(handles=legend_elements, title='', fontsize=16, title_fontsize=13, loc='lower right')

plt.tight_layout()
plt.savefig('figure-8a.pdf', dpi=1200, bbox_inches='tight')
plt.clf()

#figure b

ALDH_family = ['sp|P00352|AL1A1_HUMAN','sp|O94788|AL1A2_HUMAN','sp|Q06278|AOXA_HUMAN']
cleandata_ALDH_family = pd.merge(pd.DataFrame({'protein':ALDH_family}), cleandata, left_on='protein', right_on='protein',how='left')

cleandata_ALDH_family['tissue_type'] = cleandata_ALDH_family.tissue_type.str.replace('SQ','SC')
cleandata_ALDH_family['tissue_bmi'] = cleandata_ALDH_family['tissue_type']+'_'+cleandata_ALDH_family['bmi']
cleandata_ALDH_family['protein_cleaned'] = cleandata_ALDH_family['protein'].str.split('|').str[2].str.split('_').str[0]
cleandata_ALDH_family['log10_value'] = np.log10(cleandata_ALDH_family['value'])

cleandata_ALDH_family['protein_tissue_type'] = cleandata_ALDH_family['protein_cleaned'] + ' - ' + cleandata_ALDH_family['tissue_type']

Guo_data = pd.read_csv("figure-8b-data.csv")

fig, ax = plt.subplots(figsize=(8,6))
ax = sbn.boxplot(data=Guo_data, x='protein_tissue_type', y='value', 
    order = ['AL1A1 - OM', 'AL1A1 - SC', 'AL1A2 - OM', 'AL1A2 - SC', 'AOXA - OM', 'AOXA - SC'], 
    fill=False, linewidth=3,palette={
      'AL1A1 - OM': [153/255,0,0], 
      'AL1A1 - SC':[0,76/255,153/255], 
      'AL1A2 - OM': [153/255,0,0], 
      'AL1A2 - SC': [0,76/255,153/255], 
      'AOXA - OM': [153/255,0,0], 
      'AOXA - SC': [0,76/255,153/255]
    })
sbn.stripplot(data=Guo_data,  x='protein_tissue_type', y='value', hue="tissue_type", jitter=False,
              order = ['AL1A1 - OM', 'AL1A1 - SC', 'AL1A2 - OM', 'AL1A2 - SC', 'AOXA - OM', 'AOXA - SC'],
            linewidth=1.5, marker='o', s= 8, palette={
      'AL1A1 - OM': [153/255,0,0], 
      'AL1A1 - SC':[0,76/255,153/255], 
      'AL1A2 - OM': [153/255,0,0], 
      'AL1A2 - SC': [0,76/255,153/255], 
      'AOXA - OM': [153/255,0,0], 
      'AOXA - SC': [0,76/255,153/255],
      'OM': [153/255,0,0],
      'SC':[0,76/255,153/255]
    },alpha = 0.3)
plt.xticks(rotation=90)
plt.legend(title = None)
plt.setp(ax.get_legend().get_texts(), fontsize='16')
ax.set_ylabel('Protein expression (pmol/gm adipose tissue)',fontsize = 15)

formatter = ScalarFormatter()
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 10))
ax.yaxis.set_major_formatter(formatter)

plt.ticklabel_format(axis= 'y', style='sci',scilimits=(0,0))

ax.set_xlabel('',fontsize = 18)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(axis='x', labelsize=14)
plt.title("Targeted MRM: atRA-synthesizing enzymes", size = 16)
ax.yaxis.offsetText.set_fontsize(16)

plt.text(x=4.02, y=270, s='OM: n=14', fontsize=16, color=[153/255,0,0])
plt.text(x=4.10, y=295, s='SC: n=15', fontsize=16, color=[0,76/255,153/255])

plt.savefig('figure-8b.pdf',dpi=1200, bbox_inches='tight')
plt.clf()

# figure c

fig, ax = plt.subplots(figsize=(8,6))
ax = sbn.boxplot(data=cleandata_ALDH_family, x='protein_tissue_type', y='value', 
    order = ['AL1A1 - OM', 'AL1A1 - SC', 'AL1A2 - OM', 'AL1A2 - SC', 'AOXA - OM', 'AOXA - SC'], 
    fill=False, linewidth=3,palette={
      'AL1A1 - OM': [153/255,0,0], 
      'AL1A1 - SC':[0,76/255,153/255], 
      'AL1A2 - OM': [153/255,0,0], 
      'AL1A2 - SC': [0,76/255,153/255], 
      'AOXA - OM': [153/255,0,0], 
      'AOXA - SC': [0,76/255,153/255]
    })
sbn.stripplot(data=cleandata_ALDH_family,  x='protein_tissue_type', y='value', hue="tissue_type",
              order = ['AL1A1 - OM', 'AL1A1 - SC', 'AL1A2 - OM', 'AL1A2 - SC', 'AOXA - OM', 'AOXA - SC'],
            linewidth=1.5, marker='o', s= 8, palette={
      'AL1A1 - OM': [153/255,0,0], 
      'AL1A1 - SC':[0,76/255,153/255], 
      'AL1A2 - OM': [153/255,0,0], 
      'AL1A2 - SC': [0,76/255,153/255], 
      'AOXA - OM': [153/255,0,0], 
      'AOXA - SC': [0,76/255,153/255],
      'OM': [153/255,0,0],
      'SC':[0,76/255,153/255]
    },alpha = 0.3)
plt.xticks(rotation=90)
plt.legend(title = None)
plt.setp(ax.get_legend().get_texts(), fontsize='16')
ax.set_ylabel('Protein peak area',fontsize = 18)

formatter = ScalarFormatter()
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 10))
ax.yaxis.set_major_formatter(formatter)

plt.ticklabel_format(axis= 'y', style='sci',scilimits=(0,0))

ax.set_xlabel('',fontsize = 18)
ax.tick_params(axis='y', labelsize=16)
ax.tick_params(axis='x', labelsize=16)
plt.title("DIA: atRA-synthesizing enzymes", size = 18)
ax.yaxis.offsetText.set_fontsize(16)

plt.text(x=4.08, y=8400000000, s='OM: n=28', fontsize=16, color=[153/255,0,0])
plt.text(x=4.13, y=9200000000, s='SC: n=31', fontsize=16, color=[0,76/255,153/255])

plt.savefig('figure-8c.pdf',dpi=1200, bbox_inches='tight')
