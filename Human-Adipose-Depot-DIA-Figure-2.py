# ##### Human Adipose Depot DIA Code (Figure 2)
# ##### Yue (Winnie) Wen, Alex Zelter, Mike Riffle, Michael J. MacCoss, Nina Isoherranen
# ##### Department of Pharmaceutics, Department of Genome Science, University of Washington-Seattle
# ##### 06/12/2025

# load packages
import pandas as pd
from scipy import stats
from scipy.stats import shapiro
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn
import pyqtgraph as pg
import statsmodels.api as sm
import pingouin as pg

# panel b

adipocyte_data = pd.read_csv("figure-2-data.csv")
y1 = adipocyte_data['OMLogArea']
w1 = adipocyte_data['OM Weight']
y2 = adipocyte_data['SCLogArea']
w2 = adipocyte_data['SC Weight']
X = sm.add_constant(adipocyte_data['BMI'])
wls1 = sm.WLS(y1,X,weights=w1,missing="drop")
wls2 = sm.WLS(y2,X,weights=w2,missing="drop")
result_om = wls1.fit()
result_sc = wls2.fit()
plt.scatter(data=adipocyte_data, x='BMI', y='OMLogArea', s=adipocyte_data['OM Weight']*1000, color=[153/255,0/255,0/255],label="OM")
plt.scatter(data=adipocyte_data, x='BMI', y='SCLogArea', s=adipocyte_data['SC Weight']*1000, color=[0/255,76/255,153/255],label="SC")
plt.plot(np.unique(X), result_om.params.const + result_om.params.BMI*np.unique(X), color=[153/255,0/255,0/255])
plt.plot(np.unique(X), result_sc.params.const + result_sc.params.BMI*np.unique(X), color=[0/255,76/255,153/255])
plt.xlim(17,58)
plt.ylim(3.2,4)
plt.legend(loc='upper left', prop={'size':14})
plt.xlabel("BMI (kg/$m^2$)", size = 18)
plt.ylabel("Log adipocyte size", size = 18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.text(x=18, y=3.8, s='OM: p-value = 0.023', fontsize=13, color=[153/255,0,0])
plt.text(x=18.2, y=3.75, s='SC: p-value = 0.017', fontsize=13, color=[0,76/255,153/255])

plt.savefig('figure-2b.pdf',dpi=1200, bbox_inches='tight')
plt.clf()

# panel c

Depot_adipocyte_size_difference = 10**(adipocyte_data.SCLogArea) - 10**(adipocyte_data.OMLogArea)
sbn.histplot(x=Depot_adipocyte_size_difference,color = [76/255,0,153/255], alpha=0.3,bins="sturges", stat='count', kde=True, linewidth = 0)
plt.axvline(x=0, color = [0,0,0], linewidth = 2, linestyle ="--")
plt.axvline(x=1065, color = [76/255,0,153/255], linewidth = 2)
plt.xlabel("SC Adipocyte Size - OM Adipocyte Size (µm${^2}$)",size =16)
plt.ylabel("Frequency",size =16)
plt.xticks(fontsize =14)
plt.yticks(fontsize =14)
plt.text(x=1350, y=7.5, s='μ = 1065', fontsize=11, color=[76/255,0,153/255])
plt.ylim(0,8)
plt.savefig('figure-2c.pdf',dpi=1200, bbox_inches='tight')