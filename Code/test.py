import h5py 



import glob


import configparser
import os
import tqdm 

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('/home/ankitsingh/hr5_metalicity/params.ini')

# clusfile = parser.get('Paths','clusfile')
Fofd = parser.get('Paths','Fofdir')
outdir = parser.get('Paths','outdir')
snapfiles = parser.get('Paths','snapfiles')


files =  [filename for filename in os.listdir(snapfiles) if os.path.isfile(os.path.join(snapfiles, filename))]



snaps = [os.path.splitext(file)[0] for file in files]

# snaps.remove('296')

snaps.sort(reverse=True)

print(snaps)


for snapno in [34,35]:

    with h5py.File(f"{outdir}clusters{snapno}.hdf5", "r")as f:

        cluslist=list(f.keys())
        cluslist.remove('status')
        for clus in cluslist:
            
            for gal in f[f'/{clus}'].keys():
                
                print(gal)
                for var in list(f[f'/{clus}/{gal}/'].keys()):

                    print(f[f'/{clus}/{gal}/{var}'])