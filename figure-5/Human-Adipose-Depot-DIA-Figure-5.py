# ##### Human Adipose Depot DIA (Figure 5)
# ##### Yue (Winnie) Wen, Alex Zelter, Mike Riffle, Michael J. MacCoss, Nina Isoherranen
# ##### Department of Pharmaceutics, Department of Genome Science, University of Washington-Seattle
# ##### 06/12/2025

# load packages and function

import pandas as pd
from scipy import stats
from scipy.stats import shapiro
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn
import pyqtgraph as pg
import statsmodels.api as sm
import pingouin as pg
def getCleanColumnName(input):
    split_lst = input.split('.')
    return split_lst[1] + "_" + split_lst[2] + "_" + split_lst[4].split('_')[0]

# panel a

unnormalized_peptidedata = pd.read_csv('figure-5ab-data.tsv', sep='\t')
for i in unnormalized_peptidedata.columns[3:]:
    unnormalized_peptidedata=unnormalized_peptidedata.rename(columns={i:getCleanColumnName(i)})
unnormalized_peptidedata = pd.melt(unnormalized_peptidedata,id_vars="modifiedSequence", value_vars=unnormalized_peptidedata.columns[3:])
unnormalized_peptidedata['log2_value'] = unnormalized_peptidedata['value']
unnormalized_peptidedata['value'] = 2**(unnormalized_peptidedata['log2_value'])-1
unnormalized_peptidedata['tissue_type'] = unnormalized_peptidedata['variable'].str.split('_').str[1]
unnormalized_peptidedata['bmi'] = unnormalized_peptidedata['variable'].str.split('_').str[2]
unnormalized_peptidedata = unnormalized_peptidedata.rename(columns={'variable':'id'})
unnormalized_peptidedata['id'] = unnormalized_peptidedata['id'].str.split('_').str[0]
unnormalized_peptidedata.loc[unnormalized_peptidedata.tissue_type == "SQ", 'tissue_type'] = 'SC'
unnormalized_peptidedata['id_tissue_type'] = unnormalized_peptidedata['id'] + "-" + unnormalized_peptidedata['tissue_type']
plt.figure(figsize=(15,5))
ax=sbn.boxplot(x='id_tissue_type',y='value', data=unnormalized_peptidedata, hue='tissue_type', palette=[[0/255,76/255,153/255],[153/255,0,0]],
            order=unnormalized_peptidedata.sort_values(['tissue_type','id'])['id_tissue_type'])
for patch in ax.patches:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.6))
plt.xticks(rotation=90)
plt.yticks(fontsize=14)
plt.xticks(fontsize=10)
plt.xlabel("")
plt.ylabel("Peptide peak area", size = 16)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('figure-5a-boxplot.pdf',dpi=1200, bbox_inches='tight')
plt.clf()
plt.figure(figsize=(5,5))
sbn.kdeplot(y='value', data=unnormalized_peptidedata, fill=True, color='grey')
plt.ylabel("Peptide peak area", size = 16)
plt.savefig('figure-5a-kde-density.pdf',dpi=1200, bbox_inches='tight')
plt.clf()

# panel b
plt.figure(figsize=(15,5))
ax = sbn.boxplot(x='id_tissue_type',y='log2_value', data=unnormalized_peptidedata, hue='tissue_type', palette=[[0/255,76/255,153/255],[153/255,0,0]],
            order=unnormalized_peptidedata.sort_values(['tissue_type','id'])['id_tissue_type'])
for patch in ax.patches:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.6))
plt.xticks(rotation=90)
plt.yticks(fontsize=14)
plt.xticks(fontsize=10)
plt.xlabel("")
plt.ylabel("log$_{2}$(peptide peak area)", size = 16)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('figure-5b-boxplot.pdf',dpi=1200, bbox_inches='tight')
plt.clf()
plt.figure(figsize=(5,5))
sbn.kdeplot(y='log2_value', data=unnormalized_peptidedata, fill=True, color='grey')
plt.ylabel("log$_{2}$(peptide peak area)", size = 16)
plt.savefig('figure-5b-kde-density.pdf',dpi=1200, bbox_inches='tight')
plt.clf()

# panel c

rawpeptidedata = pd.read_csv("figure-5c-data.tsv",sep='\t')
for i in rawpeptidedata.columns[3:]:
    rawpeptidedata=rawpeptidedata.rename(columns={i:getCleanColumnName(i)})

rawpeptidedata = pd.melt(rawpeptidedata,id_vars="modifiedSequence", value_vars=rawpeptidedata.columns[3:])

rawpeptidedata['tissue_type'] = rawpeptidedata['variable'].str.split('_').str[1]
rawpeptidedata['bmi'] = rawpeptidedata['variable'].str.split('_').str[2]
rawpeptidedata = rawpeptidedata.rename(columns={'variable':'id'})
rawpeptidedata['id'] = rawpeptidedata['id'].str.split('_').str[0]
rawpeptidedata.loc[rawpeptidedata.tissue_type == "SQ", 'tissue_type'] = 'SC'
rawpeptidedata['id_tissue_type'] = rawpeptidedata['id'] + "-" + rawpeptidedata['tissue_type']
plt.figure(figsize=(15,5))
ax = sbn.boxplot(x='id_tissue_type',y='value', data=rawpeptidedata, hue='tissue_type', palette=[[0/255,76/255,153/255],[153/255,0,0]],
            order=rawpeptidedata.sort_values(['tissue_type','id'])['id_tissue_type'])
for patch in ax.patches:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.6))
plt.xticks(rotation=90)
plt.yticks(fontsize=14)
plt.xticks(fontsize=10)
plt.xlabel("")
plt.ylabel("log$_{2}$(peptide peak area)", size = 16)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('figure-5c-boxplot.pdf',dpi=1200, bbox_inches='tight')
plt.clf()
plt.figure(figsize=(5,5))
sbn.kdeplot(y='value', data=rawpeptidedata, fill=True, color='grey')
plt.ylabel("log$_{2}$(peptide peak area)", size = 16)
plt.savefig('figure-5c-kde-density.pdf',dpi=1200, bbox_inches='tight')