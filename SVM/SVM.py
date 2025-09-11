#In[]: Jr-Fong Dang
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import svm
import matplotlib.pyplot as plt
import numpy as np
df = pd.read_csv('/home/guest002/Diskarray/tung_quanzi_1_2_umap.csv')
df = df.rename({'Unnamed: 0':'Cell'},axis=1)

for i in range(0,len(df)):
    df.loc[i,'Data'] = df.loc[i,'Cell'].split('_')[0]
    df.loc[i,'Barcode'] = df.loc[i,'Cell'].split('_')[1]
    if df.loc[i,'Data']=="Tung":
        df.loc[i,'Class'] = 0
    else:
        df.loc[i,'Class'] = 1

X = df[["umap_1","umap_2"]]
y = df[["Class"]]

clf=svm.SVC(kernel='linear')
clf.fit(X,y)
clf.predict(X)
clf.support_vectors_
clf.intercept_
clf.coef_

slope = -1*(clf.coef_[0][0]/clf.coef_[0][1])

intercept = -1*(clf.intercept_[0]/clf.coef_[0][1])

umap_1 = np.linspace(-6, 2, 10)

umap_2 = slope * umap_1 + intercept
plt.ylim(-6, 6)
plt.plot(umap_1, umap_2,color='r',linestyle='--')