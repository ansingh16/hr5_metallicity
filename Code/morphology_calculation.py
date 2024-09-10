import HR5_cluster as hr5
import pandas as pd
from multiprocessing import Pool, cpu_count
import tqdm
import numpy as np
import os
import configparser
from tqdm.contrib.concurrent import process_map  # For multiprocessing

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('../params.ini')

outdir = parser.get('Paths','outdir')
galcats = parser.get('Paths','galaxycats')
morphs = parser.get('Paths','morphs')
snapdir = parser.get('Paths','snapfiles')

def process_cluster(clusid,snap):
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

    # print(f"Processing cluster {clusid}")
     # Load the morph data
    mrp = np.genfromtxt(f'{morphs}galmorf.out{snap}', dtype=values)

    # Convert morph to a DataFrame
    mrp = pd.DataFrame(mrp)
    mrp.columns = keys


    clusdict = []
    clus = hr5.Cluster(int(snap), clusid)
    for galid in clus.get_galids():
            
            try:
                mor = mrp.loc[mrp['ID']==int(galid),'sersicn'].values[0]
                clusdict.append({'clusID':clusid, 'sersicn':mor, 'ID':galid})
            except Exception as e:
                clusdict.append({'clusID':clusid, 'sersicn':-9.0, 'ID':galid})
                # print(f"Error processing {clusid} {galid} with error {e}")
                pass
                

                
            
            
                
    return clusdict
def process_snap(snap):

    print(f"Processing snapshot {snap}")
    
    try:

        clusters = pd.read_csv(f'{snapdir}/{int(snap)}.dat')


        allgal=[]
        for clusid in tqdm.tqdm(clusters.HostHaloID.values):
                
            allgal.append(process_cluster(clusid,snap))
                
        allgal = [item for sublist in allgal for item in sublist]

        morpho = pd.DataFrame(allgal)

        morpho.to_csv(f'{morphs}/morphology_{snap}.csv',index=False)

        print(f"Done processing snapshot {snap}")

    except Exception as e:
         # get error 
         print(f"Error processing snapshot {snap} with error {e}")
    

# main function
if __name__ == '__main__':

    # get snaps
    snaps = [ os.path.splitext(file)[0].zfill(3) for file in os.listdir(snapdir) if  int(os.path.splitext(file)[0])>50]

    # for snap in snaps:
    #     process_snap(snap)

    with Pool(6) as p:
        p.map(process_snap, snaps)

    # process_snap('185')

