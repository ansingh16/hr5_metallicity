import numpy as np
import pandas as pd 
import h5py
# plt.style.use('/home/ankitsingh/HR5_AGN/paper_style.mplstyle')
import tqdm 
import math
import configparser
import os
import multiprocessing

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('/home/ankitsingh/hr5_metalicity/params.ini')

# clusfile = parser.get('Paths','clusfile')
Fofd = parser.get('Paths','Fofdir')
outdir = parser.get('Paths','outdir')
snapfiles = parser.get('Paths','snapfiles')

h0=float(parser.get('Setting','h0'))

# get snapshots which are to be analysed

import glob

files =  [filename for filename in os.listdir(snapfiles) if os.path.isfile(os.path.join(snapfiles, filename))]




snaps = [os.path.splitext(file)[0] for file in files]


snaps.sort(reverse=True)

snaps.remove('99')
snaps.remove('34')
snaps.remove('35')
print(snaps)


def Correct(snapno):

    print("Processing: ",snapno)
    
    snp = str(snapno).zfill(3)

    snapdat = pd.read_parquet(f"/scratch/ankitsingh/Galaxy_catalogs/galaxy_catalogue_{snp}.parquet")


    with h5py.File(f"{outdir}clusters{snapno}.hdf5", "r") as f:
        
        cluslist = list(f.keys())
        cluslist.remove('status')
        for clus in cluslist:
            
            
                        
            grouped = snapdat.groupby('HostHaloID')

            group = grouped.get_group(int(clus))

            actmass = format(group['HostMtot(Msun)'].values[0],'.3e')

        

            clusmtot = format(f[f'/{clus}/'].attrs['mtot'],'.3e')

            if actmass != clusmtot:
                print(snapno)

            #     cluster = f[f'/{clus}/']
            
            #     for par in ['mtot','mdm','mgas','msink','mstar','pos']:

            #             cluster.attrs[par] = cluster.attrs[par]/(h0**2)
                    
            #     if clus !='ICL':
            #         for gal in list(f[f'/{clus}'].keys()):

            #             #print("galaxy ",gal)
            #             # write positions of stars in subhalo

            #             pars = ['posstar','posdm','posgas','massstar','massdm','massgas']

            #             for par in pars:
                            
            #                 dat = np.array(f[f'{clus}/{gal}/{par}'])/(h0**2)

            #                 del f[f'{clus}/{gal}/{par}']

            #                 f[f'{clus}/{gal}/{par}'] = dat

                        

            #             galaxy = f[f'/{clus}/{gal}']
                        
            #             for par in ['mtot','mdm','mgas','msink','mstar','pos']:
            #                 galaxy.attrs[par] = galaxy.attrs[par]/(h0**2)

            # else:
            #     print(snapno,clus,"already done")


Correct(35)

