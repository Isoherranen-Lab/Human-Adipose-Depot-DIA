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

cleandata = pd.read_csv("figure-6-7-8-9-data.csv")
paired_participants = [3501,3502,3503,3504,3505,3506,3507,
                               3508,3509,3510,3512,3514,3515,
                               3516,3517,3518,3519,3520,3521,3522,
                               3523,3524,3525,3527,3528,3529,
                               3531,3532]
cleandata = cleandata[cleandata.id.isin(paired_participants)]

OM_SC_paired_tests = pd.read_csv("figure-8a-data.csv")
OM_SC_paired_tests['abs_Fold_Change'] = abs(OM_SC_paired_tests['Fold_Change'])

protein_top10_OM = OM_SC_paired_tests[OM_SC_paired_tests.Fold_Change<0.5].sort_values(by=['Adjusted_P_Value','abs_Fold_Change'], ascending=[True,True]).Protein[:10]
protein_top10_SC = OM_SC_paired_tests[OM_SC_paired_tests.Fold_Change>2].sort_values(by=['Adjusted_P_Value','abs_Fold_Change'], ascending=[True,False]).Protein[:10]

cleandata_top10_OM = pd.merge(pd.DataFrame({'protein':protein_top10_OM}), cleandata, left_on='protein', right_on='protein',how='left')
cleandata_top10_SC = pd.merge(pd.DataFrame({'protein':protein_top10_SC}), cleandata, left_on='protein', right_on='protein',how='left')

cleandata_top10_OM_pivot = cleandata_top10_OM.pivot(index=['protein','id'], columns='tissue_type', values='value').reset_index()
cleandata_top10_OM_pivot['value'] = np.log2(cleandata_top10_OM_pivot.SQ/cleandata_top10_OM_pivot.OM)
cleandata_top10_OM_pivot['protein_cleaned'] = cleandata_top10_OM_pivot['protein'].str.split('|').str[2].str.split('_').str[0]

cleandata_top10_SC_pivot = cleandata_top10_SC.pivot(index=['protein','id'], columns='tissue_type', values='value').reset_index()
cleandata_top10_SC_pivot['value'] = np.log2(cleandata_top10_SC_pivot.SQ/cleandata_top10_SC_pivot.OM)
cleandata_top10_SC_pivot['protein_cleaned'] = cleandata_top10_SC_pivot['protein'].str.split('|').str[2].str.split('_').str[0]
 
order_lst = OM_SC_paired_tests[OM_SC_paired_tests.Protein.isin(cleandata_top10_OM.protein)].\
            sort_values(by=['Adjusted_P_Value','abs_Fold_Change'],ascending=[True,True]).\
            Protein.str.split('|').str[2].str.split('_').str[0].tolist() + \
                OM_SC_paired_tests[OM_SC_paired_tests.Protein.isin(cleandata_top10_SC.protein)]. \
            sort_values(by=['Adjusted_P_Value', 'abs_Fold_Change'], ascending=[True,False]). \
            Protein.str.split('|').str[2].str.split('_').str[0].tolist()
om_col = [153/255,0/255,0/255]
sc_col = [0/255,76/255,153/255]
palette = 10*[om_col] + 10*[sc_col]

adj_p_val_data = OM_SC_paired_tests[OM_SC_paired_tests.Protein.isin(cleandata_top10_OM.protein.tolist()+cleandata_top10_SC.protein.tolist())]
adj_p_val_data['protein_cleaned'] = adj_p_val_data.Protein.str.split('|').str[2].str.split('_').str[0]
adj_p_val_data['neg_log_10_value'] = -np.log10(adj_p_val_data['Adjusted_P_Value'])

figs, axes = plt.subplots(nrows=2, ncols=1, figsize= (12,9), gridspec_kw={'height_ratios': [1, 4]})
sbn.barplot(ax=axes[0], data=adj_p_val_data, x='protein_cleaned', y='Adjusted_P_Value', 
            order=order_lst, palette=palette,legend =False,width=0.7)

sbn.boxplot(ax=axes[1], data=pd.concat([cleandata_top10_OM_pivot, cleandata_top10_SC_pivot]), x='protein_cleaned', y='value',fill=False, color=[0/255,0/255,0/255],linewidth=2,
            order=order_lst,legend =False,showfliers = False)
sbn.stripplot(ax=axes[1], data=pd.concat([cleandata_top10_OM_pivot, cleandata_top10_SC_pivot]), x='protein_cleaned', y='value', palette=palette,
              s=10,alpha = 0.2, 
              order=order_lst,legend =False)
sbn.despine(left=True)
axes[1].tick_params(bottom=False)
sbn.despine(top=True, bottom=True)
axes[1].axhline(y=0, color='black')
axes[0].set_title("Top 10 Proteins Enriched in OM and SQ", size=20)
axes[0].set_ylabel('Adj. p value', fontsize=17)

axes[1].set_ylabel('log$_{2}$ fold change (SQ/OM)', fontsize=18)
axes[1].set_ylim(-7, 7)
axes[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axes[1].tick_params(axis='x', rotation=90)
axes[0].yaxis.get_offset_text().set_fontsize(18)
axes[0].ticklabel_format(style='plain', axis='y')  

for ax in axes:
    ax.set_xticks(ax.get_xticks())         
    ax.legend(title=None)
    ax.set_xlabel('', fontsize=18)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.legend_.remove()
    
legend = ax.get_legend()
if legend:
    legend.remove()
    
legend_elements = [
    Patch(facecolor=om_col, label='Enriched in OM'),
    Patch(facecolor=sc_col, label='Enriched in SQ')
]

ax.legend(handles=legend_elements, title='', fontsize=16, title_fontsize=13, loc='lower right')

plt.tight_layout()
plt.savefig('figure-8a.pdf', dpi=1200, bbox_inches='tight')
plt.savefig('figure-8a.png', dpi=1200, bbox_inches='tight')
plt.clf()

#figure b

ALDH_family = ['sp|P00352|AL1A1_HUMAN','sp|O94788|AL1A2_HUMAN','sp|Q06278|AOXA_HUMAN']
cleandata_ALDH_family = pd.merge(pd.DataFrame({'protein':ALDH_family}), cleandata, left_on='protein', right_on='protein',how='left')

cleandata_ALDH_family['tissue_bmi'] = cleandata_ALDH_family['tissue_type']+'_'+cleandata_ALDH_family['bmi']
cleandata_ALDH_family['protein_cleaned'] = cleandata_ALDH_family['protein'].str.split('|').str[2].str.split('_').str[0]
cleandata_ALDH_family['log10_value'] = np.log10(cleandata_ALDH_family['value'])

cleandata_ALDH_family['protein_tissue_type'] = cleandata_ALDH_family['protein_cleaned'] + ' - ' + cleandata_ALDH_family['tissue_type']

Guo_data = pd.read_csv("figure-8b-data.csv")
Guo_data.loc[Guo_data.tissue_type=='SC','tissue_type'] = 'SQ'
Guo_data['protein_tissue_type'] = Guo_data['protein_tissue_type'].str.replace('SC','SQ')

fig, ax = plt.subplots(figsize=(6.5,4.8))
ax = sbn.boxplot(data=cleandata_ALDH_family, x='protein_tissue_type', y=cleandata_ALDH_family['value']/10000000000, showfliers = False, 
    order = ['AL1A1 - OM', 'AL1A1 - SQ', 'AL1A2 - OM', 'AL1A2 - SQ', 'AOXA - OM', 'AOXA - SQ'], 
    fill=False, linewidth=3,palette={
      'AL1A1 - OM': [153/255,0,0], 
      'AL1A1 - SQ':[0,76/255,153/255], 
      'AL1A2 - OM': [153/255,0,0], 
      'AL1A2 - SQ': [0,76/255,153/255], 
      'AOXA - OM': [153/255,0,0], 
      'AOXA - SQ': [0,76/255,153/255]
    })
sbn.stripplot(data=cleandata_ALDH_family,  x='protein_tissue_type', y=cleandata_ALDH_family['value']/10000000000, hue="tissue_type",
              order = ['AL1A1 - OM', 'AL1A1 - SQ', 'AL1A2 - OM', 'AL1A2 - SQ', 'AOXA - OM', 'AOXA - SQ'],
            linewidth=1.5, marker='o', s= 8, palette={
      'AL1A1 - OM': [153/255,0,0], 
      'AL1A1 - SQ':[0,76/255,153/255], 
      'AL1A2 - OM': [153/255,0,0], 
      'AL1A2 - SQ': [0,76/255,153/255], 
      'AOXA - OM': [153/255,0,0], 
      'AOXA - SQ': [0,76/255,153/255],
      'OM': [153/255,0,0],
      'SQ':[0,76/255,153/255]
    },alpha = 0.3)
plt.xticks(rotation=90, fontsize = 14)
handles, labels = ax.get_legend_handles_labels()
new_labels = ["SQ (n=31)", "OM (n=28)"]
ax.legend(handles=handles[::-1], labels=new_labels[::-1])
plt.setp(ax.get_legend().get_texts(), fontsize='16')
ax.set_ylabel('Protein peak area (x10$^{10}$)',fontsize = 16)

formatter = ScalarFormatter()
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 10))
ax.yaxis.set_major_formatter(formatter)
ax.set_xlabel('',fontsize = 14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(axis='x', labelsize=14)
plt.title("DIA: atRA-synthesizing enzymes", size = 16)
ax.yaxis.offsetText.set_fontsize(16)

plt.savefig('figure-8b.pdf',dpi=1200, bbox_inches='tight')
plt.savefig('figure-8b.png',dpi=1200, bbox_inches='tight')
plt.clf()


# panel c
fig, ax = plt.subplots(figsize=(6.5,4.8))
sbn.boxplot(data=Guo_data,  x='protein_tissue_type', y='value', hue="tissue_type", showfliers=False,
              order = ['AL1A1 - OM', 'AL1A1 - SQ', 'AL1A2 - OM', 'AL1A2 - SQ', 'AOXA - OM', 'AOXA - SQ'],
            fill=False, linewidth=3, palette={
      'AL1A1 - OM': [153/255,0,0], 
      'AL1A1 - SQ':[0,76/255,153/255], 
      'AL1A2 - OM': [153/255,0,0], 
      'AL1A2 - SQ': [0,76/255,153/255], 
      'AOXA - OM': [153/255,0,0], 
      'AOXA - SQ': [0,76/255,153/255],
      'OM': [153/255,0,0],
      'SQ':[0,76/255,153/255]
    })
sbn.stripplot(data=Guo_data,  x='protein_tissue_type', y='value', hue="tissue_type", 
              order = ['AL1A1 - OM', 'AL1A1 - SQ', 'AL1A2 - OM', 'AL1A2 - SQ', 'AOXA - OM', 'AOXA - SQ'],
            linewidth=1.5, marker='o', s= 8, palette={
      'AL1A1 - OM': [153/255,0,0], 
      'AL1A1 - SQ':[0,76/255,153/255], 
      'AL1A2 - OM': [153/255,0,0], 
      'AL1A2 - SQ': [0,76/255,153/255], 
      'AOXA - OM': [153/255,0,0], 
      'AOXA - SQ': [0,76/255,153/255],
      'OM': [153/255,0,0],
      'SQ':[0,76/255,153/255]
    },alpha = 0.3)
plt.xticks(rotation=90)
handles, labels = ax.get_legend_handles_labels()
new_labels = ["SQ (n=15)", "OM (n=14)"]
ax.legend(handles=handles[::-1], labels=new_labels[::-1])
plt.setp(ax.get_legend().get_texts(), fontsize='16')
ax.set_ylabel('Protein expression (pmol/gm adipose tissue)',fontsize = 14)

ax.ticklabel_format(style='plain', axis='y') 

ax.set_xlabel('',fontsize = 18)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(axis='x', labelsize=14)
plt.title("Targeted MRM: atRA-synthesizing enzymes", size = 16)
ax.yaxis.offsetText.set_fontsize(16)

plt.savefig('figure-8c.pdf',dpi=1200, bbox_inches='tight')
plt.savefig('figure-8c.png',dpi=1200, bbox_inches='tight')
