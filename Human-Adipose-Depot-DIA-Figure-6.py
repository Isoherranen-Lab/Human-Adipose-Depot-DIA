# ##### Human Adipose Depot DIA (Figure 6)
# ##### Yue (Winnie) Wen, Alex Zelter, Mike Riffle, Michael J. MacCoss, Nina Isoherranen
# ##### Department of Pharmaceutics, Department of Genome Science, University of Washington-Seattle
# ##### 06/12/2025

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn
from matplotlib.ticker import ScalarFormatter
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import matplotlib.patches as mpatches

# figure a

def plotEllipse(x, y, ax, n=3.0, color='none', **kwargs):
    if x.size != y.size:
        raise ValueError("x and y must be the same size")
    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=color, **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n
    mean_x = np.mean(x)
    scale_y = np.sqrt(cov[1, 1]) * n
    mean_y = np.mean(y)
    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

# Proteomics dataset
proteomics_data_cleaned = pd.read_csv('figure-6-7-8-data.csv')
proteomics_data_cleaned.loc[proteomics_data_cleaned.tissue_type == 'SQ','tissue_type'] = 'SC'
protein_by_id_tt = proteomics_data_cleaned.pivot(index=['id','tissue_type','bmi'],columns='protein',values='raw_value').reset_index()

# Initialize the StandardScaler and Perform PCA
protein_col = protein_by_id_tt.columns[3:]
for i in protein_col:
    protein_by_id_tt[i] = StandardScaler().fit_transform(protein_by_id_tt[[i]])
pca = PCA(n_components=2)  # Choose the number of components
principal_components = pca.fit_transform(protein_by_id_tt[protein_col])
principal_components = pd.DataFrame(data=principal_components,columns=['PC1','PC2'])
principal_components_with_raw = pd.merge(principal_components, protein_by_id_tt, left_index=True, right_index=True)

# Plot
fig, axs = plt.subplots(1, 1, figsize=(8, 8))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
axs.scatter(data=principal_components_with_raw, x='PC1', y='PC2', c=principal_components_with_raw['tissue_type'].map({'OM':[153/255,0/255,0/255], 'SC':[0/255,76/255,153/255]}))
axs.set_xlabel("Principal Component 1 (20.0%)", size = 18)
axs.set_ylabel("Principal Component 2 (13.1%)", size = 18)
axs.legend(handles=[
    mpatches.Patch(color=[0/255,76/255,153/255] , label='SC'),
    mpatches.Patch(color=[153/255,0/255,0/255] , label='OM'),
],fontsize = 16)
plotEllipse(x=principal_components_with_raw[principal_components_with_raw.tissue_type == 'OM']['PC1'], 
    y = principal_components_with_raw[principal_components_with_raw.tissue_type == 'OM']['PC2'], 
    ax=axs, n=1.96, edgecolor=[153/255,0/255,0/255],color = [153/255,0/255,0/255], alpha = 0.1)
plotEllipse(x=principal_components_with_raw[principal_components_with_raw.tissue_type == 'SC']['PC1'], 
    y = principal_components_with_raw[principal_components_with_raw.tissue_type == 'SC']['PC2'], 
    ax=axs, n=1.96, edgecolor=[0/255,76/255,153/255],color = [0/255,76/255,153/255], alpha = 0.1)

plt.savefig('figure-6a.pdf', dpi=1200, bbox_inches='tight')
plt.clf()

# figure b

# read data file
cleandata = pd.read_csv("figure-6-7-8-data.csv")

protein_list = []
OM_mean_list =[]
SQ_mean_list = []
for i in cleandata['protein'].unique():
    temp = pd.pivot_table(cleandata[cleandata.protein==i], values='raw_value', index='id', columns='tissue_type')
    protein_list.append(i)
    OM_mean_list.append(np.mean(temp["OM"].dropna()))
    SQ_mean_list.append(np.mean(temp["SQ"].dropna()))
dynamic_range_data = pd.DataFrame({"Protein":protein_list, "OM Mean log2 PA":OM_mean_list, "SQ Mean log2 PA":SQ_mean_list})

dynamic_range_data['OM Mean log10 PA'] = np.log10(2**dynamic_range_data["OM Mean log2 PA"])
dynamic_range_data['SQ Mean log10 PA'] = np.log10(2**dynamic_range_data["SQ Mean log2 PA"])

dynamic_range_data['OM_rank'] = dynamic_range_data['OM Mean log10 PA'].rank(ascending=False).astype(int)
dynamic_range_data['SQ_rank'] = dynamic_range_data['SQ Mean log10 PA'].rank(ascending=False).astype(int)


# figure b

fig, axs = plt.subplots(1, 1, figsize=(6,6))
sbn.scatterplot(data=dynamic_range_data, x='OM_rank', y='OM Mean log10 PA',color=[153/255,0/255,0/255],alpha = 0.5,s=200)
plt.ylim(0,12)
temp1 = dynamic_range_data[dynamic_range_data.Protein.isin(['sp|P15090|FABP4_HUMAN',
                                                            'sp|O60240|PLIN1_HUMAN',
                                                            'sp|P18206|VINC_HUMAN',
                                                            'sp|P61026|RAB10_HUMAN',
                                                            'sp|Q15848|ADIPO_HUMAN',
                                                            'sp|P14672|GLUT4_HUMAN',
                                                            'sp|P41159|LEP_HUMAN / sp|LEP_HUMAN|P41159'])].reset_index(drop=True)

marker_protein_names = ['PLIN1','GLUT4','FABP4','VINC','LEP','RAB10','ADIPO']
for i in range(7):
    if i != 0:
        plt.annotate(marker_protein_names[i], (temp1['OM_rank'][i],temp1['OM Mean log10 PA'][i]),
                 xytext = (temp1['OM_rank'][i]+250,temp1['OM Mean log10 PA'][i]+0.15), arrowprops = dict(arrowstyle='-'))
    else:
        plt.annotate(marker_protein_names[i], (temp1['OM_rank'][i],temp1['OM Mean log10 PA'][i]),
                 xytext = (temp1['OM_rank'][i]+400,temp1['OM Mean log10 PA'][i]-0.18), arrowprops = dict(arrowstyle='-'))
sbn.scatterplot(data=temp1, x='OM_rank', y='OM Mean log10 PA', color=[0/255,153/255,76/255],s=200, alpha=0.8,label = "Adipose-Specific Markers")

plt.xlabel("Rank",fontsize = 14)
plt.ylabel('Log$_{10}$ protein peak area',fontsize = 14)
plt.title('Omental Adipose Tissue',fontsize = 14)
plt.legend(loc='upper right',fontsize = 12)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.text(x=3700, y=10, s='4799 proteins', fontsize=12, color=[153/255,0,0])

plt.savefig('figure-6b.pdf',dpi=1200, bbox_inches='tight')


# panel c

fig, axs = plt.subplots(1, 1, figsize=(6, 6))
sbn.scatterplot(data=dynamic_range_data, x='SQ_rank', y='SQ Mean log10 PA',color=[0/255,76/255,153/255],alpha = 0.2,s=200)
plt.ylim(0,12)
temp1 = dynamic_range_data[dynamic_range_data.Protein.isin(['sp|P15090|FABP4_HUMAN',
                                                            'sp|O60240|PLIN1_HUMAN',
                                                            'sp|P18206|VINC_HUMAN',
                                                            'sp|P61026|RAB10_HUMAN',
                                                            'sp|Q15848|ADIPO_HUMAN',
                                                            'sp|P14672|GLUT4_HUMAN',
                                                            'sp|P41159|LEP_HUMAN / sp|LEP_HUMAN|P41159'])].reset_index(drop=True)
# Based on order in the dataset
marker_protein_names = ['PLIN1','GLUT4','FABP4','VINC','LEP','RAB10','ADIPO']
for i in range(7):
    if i != 0 and i != 6:
        plt.annotate(marker_protein_names[i], (temp1['SQ_rank'][i],temp1['SQ Mean log10 PA'][i]),
                 xytext = (temp1['SQ_rank'][i]+300,temp1['SQ Mean log10 PA'][i]+0.1), arrowprops = dict(arrowstyle='-'))
    else:
        plt.annotate(marker_protein_names[i], (temp1['SQ_rank'][i],temp1['SQ Mean log10 PA'][i]),
                 xytext = (temp1['SQ_rank'][i]+400,temp1['SQ Mean log10 PA'][i]-0.35*(i/3+1)), arrowprops = dict(arrowstyle='-'))
sbn.scatterplot(data=temp1, x='SQ_rank', y='SQ Mean log10 PA', color=[0/255,153/255,76/255],s=200, alpha=0.8, label = "Adipose-Specific Markers" )

plt.xlabel("Rank",fontsize = 14)
plt.ylabel('Log$_{10}$ protein peak area',fontsize = 14)
plt.title('Subcutaneous Adipose Tissue',fontsize = 14)
plt.legend(loc='upper right',fontsize = 12)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.text(x=3700, y=10, s='4799 proteins', fontsize=12, color=[0,76/255,153/255])

plt.savefig('figure-6c.pdf',dpi=1200, bbox_inches='tight')
plt.clf()