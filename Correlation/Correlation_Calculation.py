#In[]: Correlation Jr-Fong Dang
import pandas as pd
from scipy.stats import spearmanr
from scipy import stats

laser = pd.read_csv('/home/guest002/Diskarray/laser_normal.csv')
quanzi = pd.read_csv('/home/guest002/Diskarray/tung_bio_raw_umi_500.csv')

quanzi_col = quanzi.columns
laser_col = laser.columns
intersect = quanzi_col & laser_col
laser_bar = pd.unique(laser['Barcode'])
quanzi_bar= pd.unique(quanzi['Barcode'])

for i in range(0,len(quanzi_bar)):
    qua_i = quanzi[quanzi['Barcode']==quanzi_bar[i]].reset_index()
    cor_table = [ [ "" for i in range(4) ] for j in range(len(laser_bar))]
    for j in range(0,len(laser_bar)):
        las_j = laser[laser['Barcode']==laser_bar[j]].reset_index()
        target = qua_i.filter(intersect[1:len(intersect)])
        compare = las_j.filter(intersect[1:len(intersect)])
        cor_table[j][0] = quanzi_bar[i]
        cor_table[j][1] = laser_bar[j]
        cor_table[j][2] = spearmanr(target.iloc[0,:], compare.iloc[0,:])[0]
        cor_table[j][3] = stats.pearsonr(target.iloc[0,:], compare.iloc[0,:])[0]
    cor_table = pd.DataFrame(cor_table,columns=['Barcode','Laser','Spearman','Pearson'])
    if i==0:
        cor_table.to_csv('/home/guest002/Diskarray/tung_bio_laser_correlation_final.csv',index = False)
    else:
        cor_table.to_csv('/home/guest002/Diskarray/tung_bio_laser_correlation_final.csv',index = False,header = False, mode = 'a')
    print(i)


df = pd.read_csv('/home/guest002/Diskarray/tung_bio_laser_correlation_final.csv')

for i in range(0,len(df)):
    df.loc[i,'Bio']=df.loc[i,'Laser'].split('_')[1]
    df.loc[i,'Layer']=df.loc[i,'Laser'].split('_')[2]

bar = pd.unique(df['Barcode'])

for i in range(0,len(bar)):
    bar_i = df[df['Barcode']==bar[i]]
    ly = pd.unique(bar_i['Layer'])
    cor_table = [ [ "" for i in range(4) ] for j in range(len(ly))]
    for j in range(0,len(ly)):
        ly_j=bar_i[bar_i['Layer']==ly[j]]
        cor_table[j][0] = bar[i]
        cor_table[j][1] = ly[j]
        cor_table[j][2] = ly_j['Spearman'].mean()
        cor_table[j][3] = ly_j['Pearson'].mean()
    cor_table = pd.DataFrame(cor_table,columns=['Barcode','Layer','Spearman','Pearson'])
    if i==0:
        cor_table.to_csv('/home/guest002/Diskarray/quanzi_bio_laser_layer_correlation_final.csv',index = False)
    else:
        cor_table.to_csv('/home/guest002/Diskarray/quanzi_bio_laser_layer_correlation_final.csv',index = False,header = False, mode = 'a')
    print(i)