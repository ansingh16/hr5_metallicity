from telnetlib import GA
import numpy as np
from scipy import interpolate,linalg
import numpy as np
from struct import unpack
import pandas as pd 
import h5py
from scipy.io import FortranFile
import matplotlib.pyplot as plt 
# plt.style.use('/home/ankitsingh/HR5_AGN/paper_style.mplstyle')
import tqdm 
from parallelbar import progress_map
import multiprocessing
import configparser
import os
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

# snaps.remove('296')

snaps.sort(reverse=True)

print(snaps)

# snaps=['99','34','35']



def Back(snapno):

    clusfile = f"{snapfiles}{snapno}.dat"
    clusters = pd.read_csv(clusfile,usecols=['HostHaloID'])
    clusters.sort_values('HostHaloID', inplace=True)

    print(outdir,snapno)
    with h5py.File(f"{outdir}clusters{snapno}.hdf5", "a") as f:

        if 'status' in f:
            del f['status']
            f.create_dataset('status',dtype=np.int32,data=1)
        else:
            f.create_dataset('status',dtype=np.int32,data=1)
                    
        # Open Galaxy find data files
        with open(f'{Fofd}/FoF.{snapno:0>5}/background_ptl.{snapno:0>5}', mode='rb') as file: # b -> binary

            
            hline=0
            sline=0
            kkk=0

            while True:
                
                fof = file.read(112)
                if not fof or (kkk>clusters.shape[0]-1):
                    #print(f"Done snap {snapno} with clusters: {kkk}")
                    break
                else:
                    #print(snapno,hline,kkk)
                    
                    # Check if the cluster is of interest
                    if (kkk<=clusters.shape[0]-1) & \
                        (hline == clusters['HostHaloID'].iloc[kkk]):
                        
                        # print(hline)
                        # Read fof halo
                        data2 = unpack('@6i11d',fof)
                                        
                        tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],\
                                        'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}
                    
                        #print(tsub)
                        posdm = np.empty((tsub['ndm'],3))
                        massdm = np.empty(tsub['ndm'])
                        veldm = np.empty((tsub['ndm'],3))

                        for j in range(tsub['ndm']):

                                                    
                                data3 = unpack('@13d1q1d1i1f',file.read(128))
                                tdm={'pos':data3[0:3],'vel':data3[3:6],'mass':data3[6],'dum0':data3[7],'tp':data3[8],'zp':data3[9],'mass0':data3[10],\
                                'tpp':data3[11],'indtab':data3[12],'id':data3[13],'potential':data3[14],'level':data3[15],'dum1':data3[16]}      

                                posdm[j,:] = np.array(tdm['pos']) 
                                veldm[j,:] = np.array(tdm['vel']) 
                                massdm[j] = tdm['mass']
                                                


                        posg = np.empty((tsub['ngas'],3))
                        massg = np.empty(tsub['ngas'])
                        velg = np.empty((tsub['ngas'],3))
                        tempg = np.empty(tsub['ngas'])
                        metalg = np.empty(tsub['ngas'])
                        feg = np.empty(tsub['ngas'])
                        hg = np.empty(tsub['ngas'])
                        og = np.empty(tsub['ngas'])
                            

                        if tsub['ngas']>0:
                            # Read gas particles
                            for j in range(tsub['ngas']):
                                                        
                                data4 = unpack('@4d4f1d5f1i2f1q4d',file.read(128))
                                tgas={'pos':data4[0:3],'dx':data4[3],'vel':data4[4:7],'dum0':data4[7],'density':data4[8],'temp':data4[9],'metal':data4[10],'fe':data4[11],\
                                    'h':data4[12],'o':data4[13],'level':data4[14],'mass':data4[15],'dum1':data4[16],'id':data4[17],'potential':data4[18],'f':data4[19:22]}

                                posg[j,:] = np.array(tgas['pos']) 
                                velg[j,:] = np.array(tgas['vel']) 
                                massg[j] = tgas['mass']
                                tempg[j] = tgas['temp']/tgas['density']
                                metalg[j] = tgas['metal']
                                feg[j] = tgas['fe']
                                og[j] = tgas['o']
                                hg[j] = tgas['h']

                                                

                        if tsub['nsink']>0:
                            # Read sink particles
                            for j in range(tsub['nsink']):
                                data5 = unpack('@20d2i',file.read(168))
                                tsink={'pos':data5[0:3],'vel':data5[3:6],'mass':data5[6],'tbirth':data5[7],'angm':data5[8:11],'ang':data5[11:14],'dmsmbh':data5[14:17],\
                                    'esave':data5[17],'smag':data5[18],'eps':data5[19],'id':data5[20],'dum0':data5[21]}


                        posstar = np.empty((tsub['nstar'],3))
                        massstar = np.empty(tsub['nstar'])
                        velstar = np.empty((tsub['nstar'],3))
                        zstar = np.empty(tsub['nstar'])

                        if tsub['nstar']>0:
                                                    
                            # Read star particles
                            for j in range(tsub['nstar']):
                                    
                                    data6 = unpack('@13d1q1d1i1f',file.read(128))
                                    
                                    tstar={'pos':data6[0:3],'vel':data6[3:6],'mass':data6[6],'dum0':data6[7],'tp':data6[8],'zp':data6[9],'mass0':data6[10],\
                                    
                                    'tpp':data6[11],'indtab':data6[12],'id':data6[13],'potential':data6[14],'level':data6[15],'dum1':data6[16]}

                                                        
                                    posstar[j,:] = np.array(tstar['pos']) 
                                    velstar[j,:] = np.array(tstar['vel']) 
                                    massstar[j] = tstar['mass']
                                    zstar[j] = tstar['zp']
                                    
                        # if f'{hline}/ICL/' in f:
                        #         del f[f'{hline}/ICL/'] 

                        # # Check if the dataset exists, delete it
                        # for dat in ['posstar','posdm','posgas','velstar',\
                        #     'veldm','velgas','massstar','massdm','massgas',\
                        #         'tgas','zgas','fegas','hgas','ogas','zstar']:

                            
                        #     try: 
                        #         if f'{hline}/ICL/{dat}' in f:
                        #             del f[f'{hline}/ICL/{dat}']                 
                        #     except:
                        #         pass

                        # write positions of stars in subhalo
                        f.create_dataset(f'{hline}/ICL/posstar',data=posstar/h0)
                        f.create_dataset(f'{hline}/ICL/posdm',data=posdm/h0)
                        f.create_dataset(f'{hline}/ICL/posgas',data=posg/h0)

                        # write positions of stars in subhalo
                        f.create_dataset(f'{hline}/ICL/velstar',data=velstar)
                        f.create_dataset(f'{hline}/ICL/veldm',data=veldm)
                        f.create_dataset(f'{hline}/ICL/velgas',data=velg)
                                                    
                        # Write mass

                        f.create_dataset(f'{hline}/ICL/massstar',data=massstar/h0)
                        f.create_dataset(f'{hline}/ICL/massdm',data=massdm/h0)
                        f.create_dataset(f'{hline}/ICL/massgas',data=massg/h0)


                        # write temperature
                        f.create_dataset(f'{hline}/ICL/tgas',data=tempg)

                        # write metallicities
                        f.create_dataset(f'{hline}/ICL/zgas',data=metalg)
                        f.create_dataset(f'{hline}/ICL/fegas',data=feg)
                        f.create_dataset(f'{hline}/ICL/ogas',data=og)
                        f.create_dataset(f'{hline}/ICL/hgas',data=hg)
                        f.create_dataset(f'{hline}/ICL/zstar',data=zstar)

                        # write attributes   
                        icl = f[f'/{hline}/ICL/']
                        for par in ['ndm','ngas','nsink','nstar','vel']:
                            
                            icl.attrs[par] = tsub[par]
                        for par in ['mtot','mdm','mgas','msink','mstar','pos']:
                            try:
                                icl.attrs[par] = tsub[par]/h0 
                            except:
                                icl.attrs[par] = np.array(tsub[par])/h0 

                        kkk=kkk+1

                        # print(f"Cluster found!! {kkk} with id {hline} at {snapno}")




                    else:
                        # Read fof halo
                        data2 = unpack('@6i11d',fof)
                                        
                        tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],\
                                        'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}
                    
                        #print(tsub)


                        # Read dm particles
                        for j in range(tsub['ndm']):               
                            file.read(128)
                                    
                        # Read gas particles
                        for j in range(tsub['ngas']):
                            file.read(128)

                        if tsub['nsink']>0:
                            # Read sink particles
                            for j in range(tsub['nsink']):
                                file.read(168)
                                        
                        if tsub['nstar']>0:
                            for j in range(tsub['nstar']):
                                file.read(128)
                            

                    hline=hline+1
                        

for snap in tqdm.tqdm(snaps):
  Back(snap)
