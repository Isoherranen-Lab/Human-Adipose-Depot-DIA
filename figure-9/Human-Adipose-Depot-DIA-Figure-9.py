##### Human Adipose Depot DIA (Figure 9)
##### Yue (Winnie) Wen, Alex Zelter, Mike Riffle, Michael J. MacCoss, Nina Isoherranen
##### Department of Pharmaceutics, Department of Genome Science, University of Washington-Seattle
##### 09/2025

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn

overlapped_between_OM_and_SQ = pd.read_csv(r'figure-6-7-8-9-10-data.csv')
meta_data = pd.read_csv(r"figure-7c-9-10-metadata.csv")

# panel b
temp = pd.pivot_table(overlapped_between_OM_and_SQ[overlapped_between_OM_and_SQ.protein=="sp|Q15811|ITSN1_HUMAN"], values='value', index='id', columns='tissue_type')
temp_set_SQ = pd.concat([temp.reset_index().SQ,meta_data.bmi],axis=1).dropna()
temp_set_OM = pd.concat([temp.reset_index().OM,meta_data.bmi],axis=1).dropna()
plt.figure()
g = sbn.regplot(data=temp_set_SQ, x='bmi', y=temp_set_SQ['SQ']/(10**8),color=[0,76/255,153/255], label="SQ (adj. p=0.034)")
h = sbn.regplot(data=temp_set_OM, x='bmi', y=temp_set_OM['OM']/(10**8),color=[153/255,0,0], label="OM (adj. p=0.0037)")
plt.title("sp|Q15811|ITSN1_HUMAN", size = 16)
plt.xlabel("Body Mass Index (kg/m$^{2}$)", size = 16)
plt.ylabel("Protein Peak Area (x $10^{8}$)", size = 16)
plt.legend(fontsize=13, loc='upper left')
plt.savefig('figure-9b.png',dpi=1200, bbox_inches='tight')
plt.savefig('figure-9b.pdf',dpi=1200, bbox_inches='tight')
plt.close()

# panel c
temp = pd.pivot_table(overlapped_between_OM_and_SQ[overlapped_between_OM_and_SQ.protein=="sp|P52594|AGFG1_HUMAN"], values='value', index='id', columns='tissue_type')
temp_set_SQ = pd.concat([temp.reset_index().SQ,meta_data.bmi],axis=1).dropna()
temp_set_OM = pd.concat([temp.reset_index().OM,meta_data.bmi],axis=1).dropna()
plt.figure()
g = sbn.regplot(data=temp_set_SQ, x='bmi', y=temp_set_SQ['SQ']/(10**7),color=[0,76/255,153/255], label="SQ (adj. p=0.017)")
h = sbn.regplot(data=temp_set_OM, x='bmi', y=temp_set_OM['OM']/(10**7),color=[153/255,0,0], label="OM (adj. p=0.011)")
plt.title("sp|P52594|AGFG1_HUMAN", size = 16)
plt.xlabel("Body Mass Index (kg/m$^{2}$)", size = 16)
plt.ylabel("Protein Peak Area (x $10^{7}$)", size = 16)
plt.legend(fontsize=13, loc='upper left')
plt.savefig('figure-9c.png',dpi=1200, bbox_inches='tight')
plt.savefig('figure-9c.pdf',dpi=1200, bbox_inches='tight')
plt.close()

# panel d
temp = pd.pivot_table(overlapped_between_OM_and_SQ[overlapped_between_OM_and_SQ.protein=="sp|Q86WU2|LDHD_HUMAN"], values='value', index='id', columns='tissue_type')
temp_set_SQ = pd.concat([temp.reset_index().SQ,meta_data.bmi],axis=1).dropna()
temp_set_OM = pd.concat([temp.reset_index().OM,meta_data.bmi],axis=1).dropna()
plt.figure()
h = sbn.regplot(data=temp_set_OM, x='bmi', y=temp_set_OM['OM']/(10**8),color=[153/255,0,0], label="OM (adj. p=0.011)")
g = sbn.regplot(data=temp_set_SQ, x='bmi', y=temp_set_SQ['SQ']/(10**8),color=[0,76/255,153/255], label="SQ (adj. p=0.034)")
plt.title("sp|Q86WU2|LDHD_HUMAN", size = 16)
plt.xlabel("Body Mass Index (kg/m$^{2}$)", size = 16)
plt.ylabel("Protein Peak Area (x $10^{8}$)", size = 16)
plt.legend(fontsize=13, loc='upper right')
plt.savefig('figure-9d.png',dpi=1200, bbox_inches='tight')
plt.savefig('figure-9d.pdf',dpi=1200, bbox_inches='tight')
plt.close()

# panel e
temp = pd.pivot_table(overlapped_between_OM_and_SQ[overlapped_between_OM_and_SQ.protein=="sp|P53007|TXTP_HUMAN"], values='value', index='id', columns='tissue_type')
temp_set_SQ = pd.concat([temp.reset_index().SQ,meta_data.bmi],axis=1).dropna()
temp_set_OM = pd.concat([temp.reset_index().OM,meta_data.bmi],axis=1).dropna()
plt.figure()
h = sbn.regplot(data=temp_set_OM, x='bmi', y=temp_set_OM['OM']/(10**8),color=[153/255,0,0], label="OM (adj. p=0.011)")
g = sbn.regplot(data=temp_set_SQ, x='bmi', y=temp_set_SQ['SQ']/(10**8),color=[0,76/255,153/255], label="SQ (adj. p=0.036)")
plt.title("sp|P53007|TXTP_HUMAN", size = 16)
plt.xlabel("Body Mass Index (kg/m$^{2}$)", size = 16)
plt.ylabel("Protein Peak Area (x $10^{8}$)", size = 16)
plt.legend(fontsize=13, loc='upper right')
plt.savefig('figure-9e.png',dpi=1200, bbox_inches='tight')
plt.savefig('figure-9e.pdf',dpi=1200, bbox_inches='tight')
plt.close()