import configparser
import h5py
import pandas as pd 
import os 
from parallelbar import progress_map

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('../params.ini')

# Reading the data locations
Timedat = parser.get('Paths','hr5outs')
output = parser.get('Paths','outdir')
clus_half = pd.read_csv(f'{output}/halfmass.csv')
# LBT_sort = clus_half.sort_values(by=['LBT']).groupby(['LBT'])


clusfile = parser.get('Paths','clusfile')
Fofd = parser.get('Paths','Fofdir')
clusters = pd.read_csv(clusfile)


def Tree(clus):

        
        MAHfile = pd.read_csv(f"{output}/{clus}_MAH.csv")
        
        MAHfile = MAHfile.loc[MAHfile['ClusMass(Msun)']>=1.0e12]

        for snap,haloid in zip(MAHfile['snap'].values,MAHfile['HostHaloID'].values):
                
                if os.path.exists(f"../snapfiles/{snap}.dat"):
                        with open(f"../snapfiles/{snap}.dat",'a') as f:
                                data = f.write(f"{haloid},{clus}\n")
                                f.close()
                else:
                        with open(f"../snapfiles/{snap}.dat",'w') as f:
                                f.write("HostHaloID,LastClus\n")
                                f.write(f"{haloid},{clus}\n")
                                f.close()


cluslist = clusters['HostHaloID'].unique()


progress_map(Tree,cluslist, chunk_size=1,n_cpu=4)

# read all snapfiles using pandas in the snapfiles directory and sort with respect to HostHaloID

snaplist = os.listdir('../snapfiles')
snaplist = [int(snap.split('.')[0]) for snap in snaplist]
snaplist.sort()

for snap in snaplist:
        # read the snapfiles and sort with respect to HostHaloID
        snapfile = pd.read_csv(f'../snapfiles/{snap}.dat')
        snapfile = snapfile.sort_values(by=['HostHaloID'])
        snapfile.to_csv(f'../snapfiles/{snap}.dat',index=False)

        

        

