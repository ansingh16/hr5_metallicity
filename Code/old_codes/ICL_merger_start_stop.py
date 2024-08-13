# This program reads the data from New_mergerID.csv 
# the file contains the HostHaloID, start and stop snapshots
# we are modifying the file ICL_merger.py to read the data for the mergers
# between two snapshots start and stop for a cluster

import configparser
import h5py
import pandas as pd 
import os 

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('../params.ini')

# Reading the data locations
Timedat = pd.read_csv('../Time_data.dat')
output = parser.get('Paths','outdir')
clusmer = pd.read_csv('../New_mergerID.csv')
# LBT_sort = clus_half.sort_values(by=['LBT']).groupby(['LBT'])




i=0
for clus,start,stop in zip(clusmer['HostHaloID'],clusmer['start'],clusmer['end']):

                
                MAHfile = pd.read_csv(f"{output}/{clus}_MAH.csv")

                merger = MAHfile.loc[(MAHfile['snap']>=start)&(MAHfile['snap']<=stop)]

                for snap,haloid in zip(merger['snap'].values,merger['HostHaloID'].values):
                        #print(snap,haloid)

                        
                        if os.path.exists(f"./{snap}.dat"):
                                
                                datfile = pd.read_csv(f'./{snap}.dat')
                                if haloid not in datfile.values:
                                        with open(f"./{snap}.dat",'a') as f:
                                                data = f.write(f"{haloid}\n")
                                                f.close()
                                
                        else:
                                with open(f"./{snap}.dat",'w') as f:
                                        f.write("HostHaloID\n")
                                        f.write(f"{haloid}\n")
                                        f.close()
                
                print(i)
                i=i+1

print(i)

# check if file exists append the data else create a new file
# if os.path.exists(f"{output}/Total_data.csv"):
#     Total_data.to_csv(f"{output}/Total_data.csv",mode='a',header=False)
# else: 
#     Total_data.to_csv(f"{output}/Total_data.csv",header=True)

