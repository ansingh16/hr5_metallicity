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

         
        
def process_snap(snap):

    print(f"Processing snapshot {snap}")
    
    # snapi = pd.read_parquet(f'{galcats}/galaxy_catalogue_{snap}.parquet')
    clusters = pd.read_csv(f'{snapdir}/{snap}.dat')

    snapsave = pd.DataFrame(columns=['ID','HostHaloID'])

    
    for clus in tqdm.tqdm(clusters['HostHaloID']):
        clusi = hr5.Cluster(snap,clus)

        galids = clusi.get_galids()
        galids = [int(gal) for gal in galids ]
        
        snapsave = pd.concat([snapsave,pd.DataFrame({'ID':galids,'HostHaloID':[clusi.clusID for _ in galids]})])


    snapsave.to_csv(f'{outdir}/priya_galcat/{snap}.csv',index=False)

if __name__ == '__main__':

    # get snaps
    snaps = [ os.path.splitext(file)[0] for file in os.listdir(snapdir) if  int(os.path.splitext(file)[0])>50]

    with Pool(6) as p:
        p.map(process_snap, snaps)
