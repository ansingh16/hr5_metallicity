import HR5_cluster as hr5
import pandas as pd
from multiprocessing import Pool, cpu_count
import tqdm
import numpy as np



def process_cluster(clusid):
    # print(f"Processing cluster {clusid}")
    clusdict = []
    clus = hr5.Cluster(snapshot, clusid)
    for galid in clus.get_galids():
            
           
            mor = mrp.loc[mrp['ID']==int(galid),'sersicn'].values[0]
            
            clusdict.append({'clusID':clusid, 'sersicn':mor, 'ID':galid})
    
           

    return clusdict

# Get all the IDs of clusters present at the given snapshot
snapshot = 296

# Define the instance of the class Cluster with the snapshot
clus296 = pd.read_csv('../Data/groups5e13.csv')

snap296 = pd.read_parquet('/scratch/ankitsingh/Galaxy_catalogs/galaxy_catalogue_296.parquet')


# Define the dtype for the columns in the morph file
dtype = {
        'ID': np.int32,
        'mstar': np.float32,
        'galstarmass': np.float32,
        'rmssize': np.float32,
        'asym': np.float32,
        'einasn': np.float16,
        'coni': np.float32,
        'rhalf': np.float32,
        'sersicn': np.float16,
        'vrot': np.float32,
        'vsig': np.float32,
        'sfr': np.float32
    }

keys = list(dtype.keys())
values = list(dtype.values())

# Load the morph data
mrp = np.genfromtxt('../Data/galmorf.out_296.txt', dtype=values)

# Convert morph to a DataFrame
mrp = pd.DataFrame(mrp)
mrp.columns = keys


allgal=[]
for clusid in tqdm.tqdm(clus296.HostHaloID.values):
         
    allgal.append(process_cluster(clusid))
        
allgal = [item for sublist in allgal for item in sublist]

morpho = pd.DataFrame(allgal)

morpho.to_csv('../Data/morphology.csv',index=False)
