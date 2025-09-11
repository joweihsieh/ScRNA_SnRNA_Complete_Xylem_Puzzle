#Jr-Fong Dang
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('/tung_quanzi_2_3_umap.csv')
df = df.rename({'Unnamed: 0':'Cell'},axis=1)
df_color_ref = pd.read_csv('/all_plotting_tables/plotting_TenX_Ptr_color.csv')

for i in range(0,len(df)):
    df.loc[i,'Data'] = df.loc[i,'Cell'].split('_')[0]
    df.loc[i,'Barcode'] = df.loc[i,'Cell'].split('_')[1]
    if df.loc[i,'Data'] =='Tung':
        for j in range(0,len(df_color_ref)):
            if df.loc[i,'Barcode'] == df_color_ref.loc[j,'Barcode']:
                df.loc[i,'Color'] = df_color_ref.loc[j,'Color']
        df.loc[i,'Alpha'] = 1
     elif df.loc[i,'Data'] =='Quanzi2':
        df.loc[i,'Color'] = '#C59738'
        df.loc[i,'Alpha'] = 1
    elif df.loc[i,'Data'] =='Quanzi3':
        #df.loc[i,'Color'] = '#808080'
        df.loc[i,'Color'] = '#000000'
        df.loc[i,'Alpha'] = 0

plt.scatter(df['umap_1'],df['umap_2'],label = "Wood",s=0.4,c=df['Color'],alpha=df['Alpha'])

# Layer1 Umap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

layer_tung = pd.read_csv('/tung_bio_laser_layer_correlation_final.csv')
layer_quanzi_2 = pd.read_csv('/quanzi_bio_2_laser_layer_correlation_final.csv')
layer_quanzi_3 = pd.read_csv('/quanzi_bio_3_laser_layer_correlation_final.csv')

layer_tung = layer_tung[layer_tung['Layer']=='L1']
layer_quanzi_2 = layer_quanzi_2[layer_quanzi_2['Layer']=='L1']
layer_quanzi_3 = layer_quanzi_3[layer_quanzi_3['Layer']=='L1']

frames = [layer_tung,layer_quanzi_2,layer_quanzi_3]
result = pd.concat(frames,ignore_index=True)

umap = pd.read_csv('/home/guest002/Diskarray/tung_quanzi_2_3_umap.csv')
umap = umap.rename({'Unnamed: 0':'Barcode'},axis=1)

show = result.merge(umap, left_on='Barcode', right_on='Barcode')

for i in range(0,len(show)):
    show.loc[i,'Type']=show.loc[i,'Barcode'].split('_')[0]
    if show.loc[i,'Pearson']>=0.40:
        show.loc[i,'Color'] = '#cf1111'
        show.loc[i,'Alpha'] = 1
    elif (show.loc[i,'Pearson']<0.4)&(show.loc[i,'Pearson']>=0.3):
        show.loc[i,'Color'] = '#cf1111'
        show.loc[i,'Alpha'] = 0.2
    else:
        show.loc[i,'Color'] = '#bab5b5'
        show.loc[i,'Alpha'] = 0.1

plt.scatter(show['umap_1'],show['umap_2'],label = "Wood",s=0.4,c=show['Color'],alpha=show['Alpha'])

df = show

for i in range(0,len(df)):
    if df.loc[i,'Type']=='Quanzi2':
    #if df.loc[i,'Type']=='Quanzi3':
        df.loc[i,'Alpha'] = 0

plt.scatter(df['umap_1'],df['umap_2'],label = "Wood",s=0.4,c=df['Color'],alpha=df['Alpha'])


#LCM Umap
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

umap = pd.read_csv('/home/guest002/Diskarray/tung_quanzi_2_3_umap.csv')
umap = umap.rename({'Unnamed: 0':'Barcode'},axis=1)

tung = pd.read_csv('/home/guest002/Diskarray/tung_LCM_vessel_correlation.csv')
tung = tung.rename({'Barcode':'Cell'},axis=1)

for i in range(0,len(tung)):
    tung.loc[i,'Barcode'] = tung.loc[i,'Type']+"_"+tung.loc[i,'Cell']
tung = tung.drop(['Type'],axis=1)
tung = tung.drop(['Cell'],axis=1)
tung = tung.merge(umap, left_on='Barcode', right_on='Barcode')

quanzi_2 = pd.read_csv('/home/guest002/Diskarray/quanzi_bio_2_LCM_vessel_correlation.csv')
quanzi_2 = quanzi_2.rename({'Barcode':'Cell'},axis=1)

for i in range(0,len(quanzi_2)):
    quanzi_2.loc[i,'Barcode'] = quanzi_2.loc[i,'Type']+"2_"+quanzi_2.loc[i,'Cell']
quanzi_2 = quanzi_2.drop(['Type'],axis=1)
quanzi_2 = quanzi_2.drop(['Cell'],axis=1)
quanzi_2 = quanzi_2.merge(umap, left_on='Barcode', right_on='Barcode')

quanzi_3 = pd.read_csv('/home/guest002/Diskarray/quanzi_bio_3_LCM_vessel_correlation.csv')
quanzi_3 = quanzi_3.rename({'Barcode':'Cell'},axis=1)

for i in range(0,len(quanzi_3)):
    quanzi_3.loc[i,'Barcode'] = quanzi_3.loc[i,'Type']+"3_"+quanzi_3.loc[i,'Cell']
quanzi_3 = quanzi_3.drop(['Type'],axis=1)
quanzi_3 = quanzi_3.drop(['Cell'],axis=1)
quanzi_3 = quanzi_3.merge(umap, left_on='Barcode', right_on='Barcode')

frames = [tung,quanzi_2,quanzi_3]
show = pd.concat(frames,ignore_index=True)

for i in range(0,len(show)):
    show.loc[i,'Type']=show.loc[i,'Barcode'].split('_')[0]
    #if show.loc[i,'Pearson']>0.425:
    if show.loc[i,'Pearson']>0.345:
        show.loc[i,'Color'] = '#cf1111'
        show.loc[i,'Alpha'] = 1
    #elif (show.loc[i,'Pearson']<=0.425)&(show.loc[i,'Pearson']>=0.32):
    elif (show.loc[i,'Pearson']<=0.345)&(show.loc[i,'Pearson']>=0.26):
        show.loc[i,'Color'] = '#cf1111'
        show.loc[i,'Alpha'] = 0.1
    else:
        show.loc[i,'Color'] = '#bab5b5'
        show.loc[i,'Alpha'] = 0.2

plt.scatter(show['umap_1'],show['umap_2'],label = "Wood",s=0.4,c=show['Color'],alpha=show['Alpha'])

df = show

for i in range(0,len(df)):
    if df.loc[i,'Type']=='Quanzi2':
    #if df.loc[i,'Type']=='Quanzi3':
        df.loc[i,'Alpha'] = 0


plt.scatter(df['umap_1'],df['umap_2'],label = "Wood",s=0.4,c=df['Color'],alpha=df['Alpha'])

#In[]: Umap Specific GENID
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

tung = pd.read_csv('/home/guest002/Diskarray/tung_bio_raw_umi_500.csv')
quanzi_2 = pd.read_csv('/home/guest002/Diskarray/quanzi_bio_2_raw_umi_500.csv')
quanzi_3 = pd.read_csv('/home/guest002/Diskarray/quanzi_bio_3_raw_umi_500.csv')

target_col_early =['Barcode','Potri.004G002500.v4.1','Potri.014G123400.v4.1','Potri.014G161200.v4.1',
'Potri.014G022200.v4.1',
'Potri.T011801.v4.1',
'Potri.003G100500.v4.1',
'Potri.010G141600.v4.1',
'Potri.002G029100.v4.1']
target_col_late = ['Barcode','Potri.006G222200.v4.1',
'Potri.009G012200.v4.1',
'Potri.011G047700.v4.1',
'Potri.006G129900.v4.1',
'Potri.007G120500.v4.1',
'Potri.005G256000.v4.1']

tung_rev = tung.filter(target_col_early)
quanzi_2_rev = quanzi_2.filter(target_col_early)
quanzi_3_rev = quanzi_3.filter(target_col_early)

frames = [tung_rev,quanzi_2_rev,quanzi_3_rev]
result = pd.concat(frames,ignore_index=True)

umap = pd.read_csv('/home/guest002/Diskarray/tung_quanzi_2_3_umap.csv')
umap = umap.rename({'Unnamed: 0':'Barcode'},axis=1)

show = result.merge(umap, left_on='Barcode', right_on='Barcode')

sel = target_col_early[7]
plt.hist(show[sel], 25)

for i in range(0,len(show)):
    if show.loc[i,sel]>15:
        show.loc[i,'Color'] = '#cf1111'
        show.loc[i,'Alpha'] = 1
    elif (show.loc[i,sel]<=15)&(show.loc[i,sel]>14):
        show.loc[i,'Color'] = '#cf1111'
        show.loc[i,'Alpha'] = 0.4
    else:
        show.loc[i,'Color'] = '#bab5b5'
        show.loc[i,'Alpha'] = 0.01

plt.scatter(show['umap_1'],show['umap_2'],label = "Wood",s=0.4,c=show['Color'],alpha=show['Alpha'])
