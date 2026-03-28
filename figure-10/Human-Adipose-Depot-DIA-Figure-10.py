##### Human Adipose Depot DIA (Figure 10)
##### Yue (Winnie) Wen, Alex Zelter, Mike Riffle, Michael J. MacCoss, Nina Isoherranen
##### Department of Pharmaceutics, Department of Genome Science, University of Washington-Seattle
##### 09/2025

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn

overlapped_between_OM_and_SQ = pd.read_csv(r'figure-6-7-8-9-10-data.csv')
meta_data = pd.read_csv(r"figure-7c-9-10-metadata.csv")

# panel b
temp = pd.pivot_table(overlapped_between_OM_and_SQ[overlapped_between_OM_and_SQ.protein=="sp|Q13610|PWP1_HUMAN"], values='value', index='id', columns='tissue_type')
temp_set_SQ = pd.concat([temp.reset_index().SQ,meta_data.SQ_log_area],axis=1).dropna()
temp_set_OM = pd.concat([temp.reset_index().OM,meta_data.OM_log_area],axis=1).dropna()
plt.figure()
g = sbn.regplot(data=temp_set_SQ, x='SQ_log_area', y=temp_set_SQ['SQ']/(10**7),color=[0,76/255,153/255], label="SQ (adj. p=0.988)")
h = sbn.regplot(data=temp_set_OM, x='OM_log_area', y=temp_set_OM['OM']/(10**7),color=[153/255,0,0], label="OM (adj. p=0.00099)")
plt.title("sp|Q13610|PWP1_HUMAN", size = 16)
plt.xlabel("log(adipocyte size)", size = 16)
plt.ylim(0,2.4)
plt.ylabel("Protein peak area (x $10^{7}$)", size = 16)
plt.legend(fontsize=13, loc='upper left')
plt.savefig('figure-10b.png',dpi=1200, bbox_inches='tight')
plt.savefig('figure-10b.pdf',dpi=1200, bbox_inches='tight')
plt.close()

# panel c
temp = pd.pivot_table(overlapped_between_OM_and_SQ[overlapped_between_OM_and_SQ.protein=="sp|Q8NFP9|NBEA_HUMAN"], values='value', index='id', columns='tissue_type')
temp_set_SQ = pd.concat([temp.reset_index().SQ,meta_data.SQ_log_area],axis=1).dropna()
temp_set_OM = pd.concat([temp.reset_index().OM,meta_data.OM_log_area],axis=1).dropna()
plt.figure()
g = sbn.regplot(data=temp_set_SQ, x='SQ_log_area', y=temp_set_SQ['SQ']/(10**7),color=[0,76/255,153/255], label="SQ (adj. p=0.795)")
h = sbn.regplot(data=temp_set_OM, x='OM_log_area', y=temp_set_OM['OM']/(10**7),color=[153/255,0,0], label="OM (adj. p=0.00099)")
plt.title("sp|Q8NFP9|NBEA_HUMAN", size = 16)
plt.xlabel("log(adipocyte size)", size = 16)
plt.ylabel("Protein peak area (x $10^{7}$)", size = 16)
plt.legend(fontsize=13, loc='upper left')
plt.savefig('figure-10c.png',dpi=1200, bbox_inches='tight')
plt.savefig('figure-10c.pdf',dpi=1200, bbox_inches='tight')
plt.close()

# panel d
temp = pd.pivot_table(overlapped_between_OM_and_SQ[overlapped_between_OM_and_SQ.protein=="sp|P43351|RAD52_HUMAN"], values='value', index='id', columns='tissue_type')
temp_set_SQ = pd.concat([temp.reset_index().SQ,meta_data.SQ_log_area],axis=1).dropna()
temp_set_OM = pd.concat([temp.reset_index().OM,meta_data.OM_log_area],axis=1).dropna()
plt.figure()
g = sbn.regplot(data=temp_set_SQ, x='SQ_log_area', y=temp_set_SQ['SQ']/(10**7),color=[0,76/255,153/255], label="SQ (adj. p=0.829)")
h = sbn.regplot(data=temp_set_OM, x='OM_log_area', y=temp_set_OM['OM']/(10**7),color=[153/255,0,0], label="OM (adj. p=0.0014)")
plt.title("sp|P43351|RAD52_HUMAN", size = 16)
plt.xlabel("log(adipocyte size)", size = 16)
plt.ylabel("Protein peak area (x $10^{7}$)", size = 16)
plt.legend(fontsize=13, loc='upper left')
plt.savefig('figure-10d.png',dpi=1200, bbox_inches='tight')
plt.savefig('figure-10d.pdf',dpi=1200, bbox_inches='tight')
plt.close()

# panel e
temp = pd.pivot_table(overlapped_between_OM_and_SQ[overlapped_between_OM_and_SQ.protein=="sp|P23141|EST1_HUMAN"], values='value', index='id', columns='tissue_type')
temp_set_SQ = pd.concat([temp.reset_index().SQ,meta_data.SQ_log_area],axis=1).dropna()
temp_set_OM = pd.concat([temp.reset_index().OM,meta_data.OM_log_area],axis=1).dropna()
plt.figure()
g = sbn.regplot(data=temp_set_SQ, x='SQ_log_area', y=temp_set_SQ['SQ']/(10**9),color=[0,76/255,153/255], label="SQ (adj. p=0.779)")
h = sbn.regplot(data=temp_set_OM, x='OM_log_area', y=temp_set_OM['OM']/(10**9),color=[153/255,0,0], label="OM (adj. p=0.0030)")
plt.title("sp|P23141|EST1_HUMAN", size = 16)
plt.xlabel("log(adipocyte size)", size = 16)
plt.ylabel("Protein peak area (x $10^{9}$)", size = 16)
plt.legend(fontsize=13, loc='upper left')
plt.savefig('figure-10e.png',dpi=1200, bbox_inches='tight')
plt.savefig('figure-10e.pdf',dpi=1200, bbox_inches='tight')
plt.close()