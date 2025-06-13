# ##### Human Adipose Depot DIA (Figure 7)
# ##### Yue (Winnie) Wen, Alex Zelter, Mike Riffle, Michael J. MacCoss, Nina Isoherranen
# ##### Department of Pharmaceutics, Department of Genome Science, University of Washington-Seattle
# ##### 06/12/2025

import pandas as pd
from scipy import stats
from scipy.stats import shapiro
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn
import statsmodels.api as sm
from adjustText import adjust_text
from sklearn.preprocessing import StandardScaler
import warnings
from matplotlib.patches import Patch
warnings.filterwarnings('ignore')

# figure b

def getProteinName(input):
    return input.split('|')[2].split('_')[0]

paired_participants = [3501,3502,3503,3504,3505,3506,3507,
                        3508,3509,3510,3512,3514,3515,
                        3516,3517,3518,3519,3520,3521,3522,
                        3523,3524,3525,3527,3528,3529,
                        3531,3532]
cleandata = pd.read_csv('figure-6-7-8-data.csv')
cleandata = cleandata[cleandata.id.isin(paired_participants)]

tested_protein_list = []
test_statistics_list = []
OM_mean_PA = []
SC_mean_PA = []
p_value_list = []
fold_change_list = []
for i in cleandata['protein'].unique():
    temp = pd.pivot_table(cleandata[cleandata.protein==i], values='value', index='id', columns='tissue_type')
    stat1,p1 = stats.wilcoxon(temp.OM, temp.SQ, nan_policy='omit', zero_method='pratt')
    OM_mean_PA.append(np.mean(temp.OM))
    SC_mean_PA.append(np.mean(temp.SQ))
    test_statistics_list.append(stat1)
    p_value_list.append(p1)
    tested_protein_list.append(i)
    fold_change_list.append(np.mean(np.log2(temp.SQ/temp.OM)))
OM_SC_paired_tests = pd.DataFrame({"Protein":tested_protein_list, "Mean OM Protein Abundance":OM_mean_PA, "Mean SC Protein Abundance": SC_mean_PA, "log2(Fold Change)":fold_change_list, "P Value":p_value_list})
OM_SC_paired_tests["Adjusted P Value"] = stats.false_discovery_control(OM_SC_paired_tests["P Value"])
OM_SC_paired_tests = OM_SC_paired_tests.sort_values(by="Adjusted P Value", ascending=True)
OM_SC_paired_tests["Neg_log10_Adjusted_P_Value"] = -np.log10(OM_SC_paired_tests["Adjusted P Value"])

threshold = -np.log10(0.05)

plt.scatter(data=OM_SC_paired_tests, x='log2(Fold Change)', y='Neg_log10_Adjusted_P_Value', label=None, color=[160/255,160/255,160/255],alpha=0.2,s=20)
plt.axhline(threshold,color="black",linestyle="--")
plt.axvline(1,color="black",linestyle="--")
plt.axvline(-1,color="black",linestyle="--")
plt.xlim(-10,10)
plt.ylim(-0.1,6.2)
plt.xlabel("log$_{2}$(fold change)", size = 18)
plt.ylabel("-log$_{10}$(adj. p-value)", size = 18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.locator_params(axis='y', nbins=4)
plt.locator_params(axis='x', nbins=4)
OM_is_greater = OM_SC_paired_tests[(OM_SC_paired_tests['Neg_log10_Adjusted_P_Value']>threshold)&(OM_SC_paired_tests['log2(Fold Change)']<-1)]
SC_is_greater = OM_SC_paired_tests[(OM_SC_paired_tests['Neg_log10_Adjusted_P_Value']>threshold)&(OM_SC_paired_tests['log2(Fold Change)']>1)]

OM_is_greater_points = OM_is_greater.sort_values(by=['Neg_log10_Adjusted_P_Value','log2(Fold Change)'], ascending=[False,True])[:10].reset_index(drop=True)
SC_is_greater_points = SC_is_greater.sort_values(by=['Neg_log10_Adjusted_P_Value','log2(Fold Change)'], ascending=[False,True])[:10].reset_index(drop=True)
OM_high_fold_change_points = OM_is_greater[(OM_is_greater["log2(Fold Change)"] < -4)&(OM_is_greater.Neg_log10_Adjusted_P_Value < 4)].reset_index(drop=True)
SC_high_fold_change_points = SC_is_greater[(SC_is_greater["log2(Fold Change)"] > 3)&(SC_is_greater.Neg_log10_Adjusted_P_Value < 4)].reset_index(drop=True)

OM_is_smaller = OM_SC_paired_tests[(OM_SC_paired_tests['Neg_log10_Adjusted_P_Value']>threshold)&(OM_SC_paired_tests['log2(Fold Change)']>-1)&(OM_SC_paired_tests['log2(Fold Change)']<0)]
SC_is_smaller = OM_SC_paired_tests[(OM_SC_paired_tests['Neg_log10_Adjusted_P_Value']>threshold)&(OM_SC_paired_tests['log2(Fold Change)']<1)&(OM_SC_paired_tests['log2(Fold Change)']>0)]

plt.scatter(x=OM_is_greater['log2(Fold Change)'],y=OM_is_greater['Neg_log10_Adjusted_P_Value'],edgecolor='none',label="Enriched in OM",color=[153/255,0/255,0/255],s=30)
plt.scatter(x=SC_is_greater['log2(Fold Change)'],y=SC_is_greater['Neg_log10_Adjusted_P_Value'],edgecolor='none',label="Enriched in SC",color=[0/255,76/255,153/255],s=30)
plt.scatter(x=OM_is_smaller['log2(Fold Change)'],y=OM_is_smaller['Neg_log10_Adjusted_P_Value'],edgecolor='none',color=[255/255,204/255,204/255],s=30)
plt.scatter(x=SC_is_smaller['log2(Fold Change)'],y=SC_is_smaller['Neg_log10_Adjusted_P_Value'],edgecolor='none',color=[153/255,204/255,255/255],s=30)

legend = plt.legend()
plt.setp(legend.get_texts(), fontsize=16)

OM_x_axis = [-6.3,-6.7,-7.2,-6.9,-5.8,-3.5,-2.9,-6,-3.7,-4.7]
OM_y_axis = [5.25,4.7,4.9,5.5,5.7,5.7,5.9,4.5,4.2,4.45]
SC_x_axis = [1.2, 2.1, 3.5, 3.2, 5.5, 5.0,5.2,5.0,4.5,2.1]
SC_y_axis = [5.9, 5.7, 5.5, 5.25, 5.15, 4.9,4.6,4.3,4.0,4.1]

for i in range(10):
    plt.annotate(getProteinName(OM_is_greater_points['Protein'][i]), 
                 (OM_is_greater_points['log2(Fold Change)'][i], OM_is_greater_points['Neg_log10_Adjusted_P_Value'][i]), 
                 xytext = (OM_x_axis[i], OM_y_axis[i]), arrowprops = dict(arrowstyle='-'))

for i in range(10):
    plt.annotate(getProteinName(SC_is_greater_points['Protein'][i]), 
                 (SC_is_greater_points['log2(Fold Change)'][i], SC_is_greater_points['Neg_log10_Adjusted_P_Value'][i]), 
                 xytext = (SC_x_axis[i], SC_y_axis[i]), arrowprops = dict(arrowstyle='-'))
adjust_text([plt.text(SC_high_fold_change_points['log2(Fold Change)'][i], 
                      SC_high_fold_change_points['Neg_log10_Adjusted_P_Value'][i],
                      getProteinName(SC_high_fold_change_points['Protein'][i])) for i in range(len(SC_high_fold_change_points))], only_move={'points':'y','texts':'y'}, arrowprops = dict(arrowstyle='-'))
adjust_text([plt.text(OM_high_fold_change_points['log2(Fold Change)'][i], 
                      OM_high_fold_change_points['Neg_log10_Adjusted_P_Value'][i],
                      getProteinName(OM_high_fold_change_points['Protein'][i])) for i in range(len(OM_high_fold_change_points))], only_move={'points':'y','texts':'y'}, arrowprops = dict(arrowstyle='-'))

plt.legend(loc='lower left', prop={'size':10})

plt.savefig('figure-7b.pdf',dpi=1200, bbox_inches='tight')
plt.clf()

# figure c

def getProteinName(input):
    return input.split('|')[2].split('_')[0]


tissue_type_feature = pd.read_csv('figure-7c-data.csv')
proteomics_data_cleaned = pd.read_csv('figure-6-7-8-data.csv')
tissue_type_wsrtest = proteomics_data_cleaned[proteomics_data_cleaned.protein.isin(tissue_type_feature.Protein.unique())]
tissue_type_wsrtest['Protein_ID'] = tissue_type_wsrtest.protein.str.split("|").str[1]
tissue_type_wsrtest = tissue_type_wsrtest[(tissue_type_wsrtest.id != 3504) | (tissue_type_wsrtest.tissue_type != "SQ")]

for i in tissue_type_wsrtest['protein'].unique():
    scalar = StandardScaler()
    tissue_type_wsrtest.loc[tissue_type_wsrtest['protein']==i, 'raw_value'] = scalar.fit_transform(tissue_type_wsrtest[tissue_type_wsrtest['protein']==i][['raw_value']])
tissue_type_wsrtest.loc[tissue_type_wsrtest.tissue_type=="SQ","tissue_type"]="SC"

metafile = pd.read_csv('figure-7c-metadata.csv')
metafile['ID'] = metafile['ID'].str[3:].astype('int')
metafile['age group'] = 'age 20-40'
metafile.loc[np.logical_and(metafile.age >= 40, metafile.age < 60), 'age group'] = 'age 40-60'
metafile.loc[metafile.age >= 60, 'age group'] = 'age≥60'
metafile['bmi group'] = 'bmi<25'
metafile.loc[np.logical_and(metafile.bmi >= 25, metafile.bmi < 30), 'bmi group'] = 'bmi 25-30'
metafile.loc[metafile.bmi >= 30, 'bmi group'] = 'bmi≥30'

tissue_type_wsrtest_data = pd.merge(tissue_type_wsrtest, metafile[['ID','sex','age group','bmi group']], left_on='id', right_on='ID', how='left')

tissue_type_wsrtest_data['protein_name'] = tissue_type_wsrtest_data['protein'].str.split('|').str[2].str.split('_').str[0]
tissue_type_wsrtest_data_protein = tissue_type_wsrtest_data.pivot(index = ['id','tissue_type'], columns = ['protein_name'], values= 'raw_value').reset_index()
tissue_type_wsrtest_data_protein = tissue_type_wsrtest_data_protein.set_index(['id','tissue_type'])
tissue_type_wsrtest_data_protein.head(20)

tissue_type_wsrtest_data_final = pd.merge(tissue_type_wsrtest_data_protein,
                                          tissue_type_wsrtest_data[['sex','id','age group','bmi group','tissue_type']].drop_duplicates().set_index(['id','tissue_type'], drop=False), 
    left_index=True, right_index=True, how='left').drop(columns=['id'])
tissue_type_wsrtest_data_final.head(3)

tissue_type_wsrtest_data_final = tissue_type_wsrtest_data_final.rename(columns={'tissue_type': 'tissue type'})

def clustermap_wrap_rotated(data_input,label_features, colors):
    data = data_input.assign(blank=0)

    data_t = data.drop(label_features+['blank'], axis=1)
    data_t = data_t.reset_index()
    data_t['id_tissuetype'] = data_t['id'].astype('str') + '-' + data_t['tissue_type']
    data_t.drop(columns=['id','tissue_type'], inplace=True)
    data_t = data_t.set_index('id_tissuetype').transpose()

    labels = pd.DataFrame()
    legends = []
    for l, c in zip(label_features, colors):
        lut = dict(zip(data[l].unique(), c))
        if l == label_features[0]:
            labels = data[l].map(lut)
            labels.index = data_t.columns
        else:
            temp_label = data[l].map(lut)
            temp_label.index = data_t.columns
            labels = pd.concat([labels, temp_label], axis=1)
        for i,j in zip(data[l].unique(), c):
            legends.append(Patch(color=j, label=i))
    labels = pd.concat([labels, data['blank'].map(dict(zip(data['blank'],['white']))).rename('')], axis=1)  

    g = sbn.clustermap(data_t, col_colors=labels, col_cluster=True, colors_ratio=0.022,
                       cmap = sbn.color_palette("mako_r", as_cmap=True),yticklabels=False,xticklabels=1,figsize=(12,14))
    
    g.fig.subplots_adjust(right=0.7)
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.set_xlabel("")
    g.tick_params(bottom=False,right=False)
    g.ax_cbar.set_position((0.73, 0.45, .03, .2))
    plt.legend(handles=legends, bbox_to_anchor = (3.5, -0.3),fontsize = 12)
    plt.text(x=0,y=5,s='Z Score', fontsize=12)
    plt.savefig('figure-7c.pdf',dpi=1200, bbox_inches='tight')


clustermap_wrap_rotated(tissue_type_wsrtest_data_final, ['tissue type','sex','age group','bmi group'],
    [
        [(153/255,0,0),(0,76/255,153/255)],
        [(204/255,153/255,1),(51/255,0,102/255)],
        [(153/255,76/255,0),(1,128/255,0),(1,178/255,102/255)],
        [(51/255,102/255,0),(102/255,102/255,0),(178/255,1,102/255)]
    ])



