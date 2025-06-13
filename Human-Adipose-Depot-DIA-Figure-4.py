# ##### Human Adipose Depot DIA Code (Figure 4)
# ##### Yue (Winnie) Wen, Alex Zelter, Mike Riffle, Michael J. MacCoss, Nina Isoherranen
# ##### Department of Pharmaceutics, Department of Genome Science, University of Washington-Seattle
# ##### 06/12/2025

import pandas as pd
from scipy import stats
from scipy.stats import shapiro
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn
import pyqtgraph as pg
import statsmodels.api as sm
import pingouin as pg

# panel a

data = pd.read_csv("figure-4-data.csv")
data['tissue'] = data['Sample Name'].str.split('-').str[2]
data['bmi'] = data['Sample Name'].str.split('-').str[-1]
boxplot_data1 = data[['Peptide Detections', 'Protein Detections', 'tissue']]
boxplot_data1.columns = ['Peptide Detections', 'Protein Detections', 'category']
boxplot_data2 = data[['Peptide Detections', 'Protein Detections', 'bmi']]
boxplot_data2.columns = ['Peptide Detections', 'Protein Detections', 'category']
boxplot_data = pd.concat(
    [boxplot_data1, boxplot_data2]
)
boxplot_data.loc[boxplot_data.category == 'SQ', 'category'] = 'SC'
fig, ax = plt.subplots(figsize=(6,6))
sbn.boxplot(x=boxplot_data['category'], y=boxplot_data['Peptide Detections'], order=['OM','SC','obese','lean'],fill=False,linewidth=2,color='black')
sbn.stripplot(x=boxplot_data['category'], y=boxplot_data['Peptide Detections'], order=['OM','SC','obese','lean'],s=10,alpha = 0.3,palette=[[153/255,0,0], [0,76/255,153/255], [153/255,0,153/255], [0,153/255,76/255]])
plt.title("Peptide Detections", size = 20)
plt.ylabel('Peptide Detections',fontsize = 18)
plt.xlabel('',fontsize = 18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.savefig('figure-4a.pdf',dpi=1200, bbox_inches='tight')

# panel b

plt.clf()
fig, ax = plt.subplots(figsize=(6,6))

sbn.boxplot(x=boxplot_data['category'], y=boxplot_data['Protein Detections'], order=['OM','SC','obese','lean'],fill=False,linewidth=2,color='black')
sbn.stripplot(x=boxplot_data['category'], y=boxplot_data['Protein Detections'], order=['OM','SC','obese','lean'],s=10,alpha = 0.3,palette=[[153/255,0,0], [0,76/255,153/255], [153/255,0,153/255], [0,153/255,76/255]])

plt.title("Protein Detections", size = 20)
plt.ylabel('Protein Detections',fontsize = 18)
plt.xlabel('',fontsize = 18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.savefig('figure-4b.pdf',dpi=1200, bbox_inches='tight')