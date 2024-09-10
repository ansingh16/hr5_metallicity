import pandas as pd
import HR5_cluster as hr5
import configparser
from tqdm.contrib.concurrent import process_map  # For multiprocessing
import os 
from multiprocessing import Pool
import numpy as np 

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('../params.ini')

outdir = parser.get('Paths','outdir')
galcats = parser.get('Paths','galaxycats')
morphs = parser.get('Paths','morphs')
snapdir = parser.get('Paths','snapfiles')


# get files in snapfiles directory

snaps = [ int(os.path.splitext(file)[0]) for file in os.listdir(snapdir) if  int(os.path.splitext(file)[0])>50]

def process_snap(snap):

    snapm = str(snap).zfill(3)
    # get only elliptical galaxies
    morpho = pd.read_csv(f'{morphs}/morphology_{snapm}.csv')

    e_galaxies = morpho[morpho['sersicn']>2.5]

    e_galaxies.set_index('ID',inplace=True)

    # clusid = 1664541

    Ana = hr5.Analysis(snap)

    Ana.get_slope_data(galids=e_galaxies.index,clusids=e_galaxies.clusID,rmax=4,rbin_width=0.3,var='feh',dump='True')


    slope_median = Ana.slope_df['slope_feh'].median()


    return snap,slope_median


# main function
if __name__ == '__main__':

    # time_df = pd.read_csv('../Data/Time_data.csv')
    
    # get snaps
    snaps = [ int(os.path.splitext(file)[0]) for file in os.listdir(snapdir) if  int(os.path.splitext(file)[0])>50]


    # for snap in snaps:
    #     process_snap(snap)
    # snaps=[101]
    # process_snap(101)
    with Pool(6) as p:
        res = p.map(process_snap, snaps)

    # convert to dataframe
    df = pd.DataFrame(res, columns=['snapshot', 'slope'])
    df['snapshot'] = df['snapshot'].astype(np.int64)
    df.set_index('snapshot',inplace=True)
    df['slope'] = df['slope'].astype(np.float16)

    # get snapshot data
    time_df = pd.read_csv('../Data/Time_data.csv')
    time_df.columns = ['snapshot','redshift','LBT','dx']
    time_df.set_index('snapshot',inplace=True)

    # get final data
    df = df.join(time_df,how='left')


    df.to_csv(f'{outdir}/slope_redshift.csv',index=False)
    

