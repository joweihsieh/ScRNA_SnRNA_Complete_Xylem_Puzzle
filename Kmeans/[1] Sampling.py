#Jr-Fong Dang
import pandas as pd
import matplotlib.pyplot as plt

umap = pd.read_csv('/home/guest002/Diskarray/tung_quanzi_2_3_umap.csv')
umap = umap.rename({'Unnamed: 0':'Cell'},axis=1)
for i in range(0,len(umap)):
    umap.loc[i,'Data'] = umap.loc[i,'Cell'].split('_')[0]
    umap.loc[i,'Barcode'] = umap.loc[i,'Cell'].split('_')[1]

umap_quanzi = umap[(umap['Data']=="Quanzi2")|(umap['Data']=="Quanzi3")]
area_1 = umap_quanzi [(umap_quanzi ['umap_1']<=-2.1)&(umap_quanzi ['umap_2']>=2)&(umap_quanzi ['umap_2']<=5.2)].sample(100)
area_2 = umap_quanzi [(umap_quanzi ['umap_1']>=-2.1)&(umap_quanzi ['umap_1']<=0)&(umap_quanzi ['umap_2']>=2)&(umap_quanzi ['umap_2']<=5.2)].sample(49)
area_3 = umap_quanzi [(umap_quanzi ['umap_1']>=0)&(umap_quanzi ['umap_1']<=2.1)&(umap_quanzi ['umap_2']>=2)&(umap_quanzi ['umap_2']<=5.2)].sample(54)
area_4 = umap_quanzi [(umap_quanzi ['umap_1']>=2.1)&(umap_quanzi ['umap_1']<=4.2)&(umap_quanzi ['umap_2']>=2)&(umap_quanzi ['umap_2']<=5.2)].sample(27)
area_5 = umap_quanzi [(umap_quanzi ['umap_1']<=-2.1)&(umap_quanzi ['umap_2']>=-2)&(umap_quanzi ['umap_2']<=2)].sample(100)
area_6 = umap_quanzi [(umap_quanzi ['umap_1']>=-2.1)&(umap_quanzi ['umap_1']<=0)&(umap_quanzi ['umap_2']>=-2)&(umap_quanzi ['umap_2']<=2)].sample(67)
area_7 = umap_quanzi [(umap_quanzi ['umap_1']>=0)&(umap_quanzi ['umap_1']<=2.1)&(umap_quanzi ['umap_2']>=-2)&(umap_quanzi ['umap_2']<=2)].sample(67)
area_8 = umap_quanzi [(umap_quanzi ['umap_1']>=2.1)&(umap_quanzi ['umap_1']<=4.2)&(umap_quanzi ['umap_2']>=-2)&(umap_quanzi ['umap_2']<=2)].sample(60)
area_9 = umap_quanzi [(umap_quanzi ['umap_1']<=-2.1)&(umap_quanzi ['umap_2']>=-5.2)&(umap_quanzi ['umap_2']<=-2)].sample(27)
area_10 = umap_quanzi [(umap_quanzi ['umap_1']>=-2.1)&(umap_quanzi ['umap_1']<=0)&(umap_quanzi ['umap_2']>=-5.2)&(umap_quanzi ['umap_2']<=-2)].sample(27)
area_11 = umap_quanzi [(umap_quanzi ['umap_1']>=0)&(umap_quanzi ['umap_1']<=2.1)&(umap_quanzi ['umap_2']>=-5.2)&(umap_quanzi ['umap_2']<=-2)].sample(54)
area_12 = umap_quanzi [(umap_quanzi ['umap_1']>=2.1)&(umap_quanzi ['umap_1']<=4.2)&(umap_quanzi ['umap_2']>=-5.2)&(umap_quanzi ['umap_2']<=-2)].sample(44)
target = umap_quanzi [(umap_quanzi ['umap_1']<=-2.1)&(umap_quanzi ['umap_2']>=-2)&(umap_quanzi ['umap_2']<=2)]
frames = [area_1,area_2,area_3,area_4,area_5,area_6,area_7,area_8,area_9,area_10,area_11,area_12]
result = pd.concat(frames,ignore_index=True)