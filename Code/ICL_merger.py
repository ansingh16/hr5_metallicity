import configparser
import h5py
import pandas as pd 
import os 

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('./params.ini')

# Reading the data locations
Timedat = pd.read_csv('./Time_data.csv')
output = parser.get('Paths','outdir')
clus_half = pd.read_csv(f'{output}/halfmass.csv')
# LBT_sort = clus_half.sort_values(by=['LBT']).groupby(['LBT'])


clusfile = parser.get('Paths','clusfile')
Fofd = parser.get('Paths','Fofdir')
clusters = pd.read_csv(clusfile)



i=0
for clus in clusters['HostHaloID'].unique():

                Red = clus_half.loc[clus_half['HostHaloID']==clus,'Redshift'].values[0]
                half_snap = Timedat.loc[Timedat['Redshift']==Red,'Snapshot'].values[0]

                MAHfile = pd.read_csv(f"{output}/{clus}_MAH.csv")

                beyond_half = MAHfile[MAHfile['snap']>half_snap]

                for snap,haloid in zip(beyond_half['snap'].values,beyond_half['HostHaloID'].values):
                        print(snap,haloid)
                        if os.path.exists(f"./{snap}.dat"):
                                with open(f"./{snap}.dat",'a') as f:
                                        data = f.write(f"{haloid}\n")
                                        f.close()
                        else:
                                with open(f"./{snap}.dat",'w') as f:
                                        f.write("HostHaloID\n")
                                        f.write(f"{haloid}\n")
                                        f.close()
                
                i=i+1

print(i)

# check if file exists append the data else create a new file
# if os.path.exists(f"{output}/Total_data.csv"):
#     Total_data.to_csv(f"{output}/Total_data.csv",mode='a',header=False)
# else: 
#     Total_data.to_csv(f"{output}/Total_data.csv",header=True)

