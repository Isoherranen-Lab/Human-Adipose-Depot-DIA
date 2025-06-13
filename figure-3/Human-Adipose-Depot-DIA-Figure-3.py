# ##### Human Adipose Depot DIA Code (Figure 3)
# ##### Yue (Winnie) Wen, Alex Zelter, Mike Riffle, Nina Isoherranen
# ##### Department of Pharmaceutics, Department of Genome Science, University of Washington-Seattle
# ##### 06/12/2025

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn
from matplotlib.ticker import ScalarFormatter

# figure a

peptide_CV_all = pd.read_csv('figure-3a-data.csv')

fig, ax = plt.subplots(figsize=(4.5,4.5))
ax = sbn.boxplot(data=peptide_CV_all, x='Protein', y='All CV', 
    fill=False, linewidth=3, color = [0,0,0],
    showmeans=True,
    meanline=True,
    showfliers=False,
    meanprops={'color': 'black', 'linestyle': '-', 'linewidth': 2},
    medianprops={'visible': False},
    whiskerprops={'visible': False},
    showcaps=False,
    boxprops={'visible': False})
plt.ylim(0,40)
sbn.swarmplot(data=peptide_CV_all, x='Protein', y='All CV', 
            linewidth=1.5, marker='o', s= 10, alpha = 0.9, color = [76/255,0,153/255])
plt.xticks(rotation=90)
plt.title("Process Controls\n(spiked into samples)",fontsize=14)
ax.set_ylabel('Coefficient of Variation (%)',fontsize = 12)

formatter = ScalarFormatter()
formatter.set_scientific(False)
formatter.set_powerlimits((-1, 10))
ax.yaxis.set_major_formatter(formatter)

ax.set_xlabel('',fontsize = 16)
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='x', labelsize=12)
ax.yaxis.offsetText.set_fontsize(16)

plt.savefig('figure-3a.pdf',dpi=1200, bbox_inches='tight')
plt.clf()

# panel b

System_Suitability_CV = pd.read_csv('figure-3b-data.csv')
System_Suitability_CV.head(3)

fig, ax = plt.subplots(figsize=(3,4.5))
ax = sbn.boxplot(data=System_Suitability_CV, x='Protein', y='All CV', 
    fill=False, linewidth=3, color = [0,0,0],
    showmeans=True,
    meanline=True,
    showfliers=False,
    meanprops={'color': 'black', 'linestyle': '-', 'linewidth': 2},
    medianprops={'visible': False},
    whiskerprops={'visible': False},
    showcaps=False,
    boxprops={'visible': False})
plt.ylim(0,40)
sbn.swarmplot(data=System_Suitability_CV, x='Protein', y='All CV', 
            linewidth=1.5, marker='o', s= 10, alpha = 0.9, color = [76/255,0,153/255])
plt.xticks(rotation=90)
plt.title("System Suitability\n(sample-independent)",fontsize=14)
ax.set_ylabel('Coefficient of Variation (%)',fontsize = 12)

formatter = ScalarFormatter()
formatter.set_scientific(False)
formatter.set_powerlimits((-1, 10))
ax.yaxis.set_major_formatter(formatter)

ax.set_xlabel('',fontsize = 16)
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='x', labelsize=12)
ax.yaxis.offsetText.set_fontsize(16)

plt.savefig('figure-3b.pdf',dpi=1200, bbox_inches='tight')