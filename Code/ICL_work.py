import numpy as np
from struct import unpack
import pandas as pd 
import h5py
from scipy.io import FortranFile
import matplotlib.pyplot as plt 
# plt.style.use('/home/ankitsingh/HR5_AGN/paper_style.mplstyle')
import tqdm
import sys
import configparser
  
import os 
import glob


files = glob.glob('*.dat')


# load up the parameter file
parser = configparser.ConfigParser()
parser.read('params.ini')


# clusfile = parser.get('Paths','clusfile')
Fofd = parser.get('Paths','Fofdir')
output = parser.get('Paths','outdir')


snapno=sys.argv[1]

print(f'Processing snap no. {snapno}')
clusfile = f"{snapno}.dat"
clusters = pd.read_csv(clusfile,usecols=['HostHaloID'])
clusters.sort_values('HostHaloID', inplace=True,ignore_index=True)

with h5py.File(f"{output}clusters{snapno}.hdf5", "w") as f:
            
            
            # Open Galaxy find data files
            with open(f'{Fofd}FoF.{snapno:0>5}/GALFIND.DATA.{snapno:0>5}', mode='rb') as file: # b -> binary

                hline=0
                sline=0
                kkk=0

                while True:
                    
                    fof = file.read(112)
                    if not fof or (kkk>clusters.shape[0]-1):
                        print("Done!! Smell the success.. :)",kkk,hline)
                        break
                    else:
                        #print(hline,kkk,clusters.loc[kkk,'HostHaloID'])
                        
                        # Check if the cluster is of interest
                        if clusters['HostHaloID'].iloc[kkk]==hline:
                                        
                        
                                        # Read fof halo
                                        data1 = unpack('@6i11d',fof)
                                        
                                        thalo={'nsub':data1[0],'ndm': data1[1],'nstar': data1[2],'nsink':data1[3],'ngas':data1[4],'npall':data1[5],'mtot':data1[6],\
                                        'mdm':data1[7],'mgas':data1[8],'msink':data1[9],'mstar':data1[10],'pos':data1[11:14],'vel':data1[14:17]}
                                        
                                        #print(thalo)

                                        # Read subhaloes
                                        for i in range(thalo['nsub']):
                                                data2 = unpack('@6i11d',file.read(112))
                                                tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],\
                                                'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}

                                                #print(hline,thalo['pos'],tsub['pos'])
                                                
                                                posdm = np.empty((tsub['ndm'],3))
                                                massdm = np.empty(tsub['ndm'])
                                                veldm = np.empty((tsub['ndm'],3))

                                                # Read dm particles
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
                                                        #print('gas',j,tsub)
                                                        

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
                                    
                                                
                                                

                                                # write positions of stars in subhalo
                                                f.create_dataset(f'{hline}/{sline}/posstar',data=posstar)
                                                f.create_dataset(f'{hline}/{sline}/posdm',data=posdm)
                                                f.create_dataset(f'{hline}/{sline}/posgas',data=posg)

                                                # write positions of stars in subhalo
                                                f.create_dataset(f'{hline}/{sline}/velstar',data=velstar)
                                                f.create_dataset(f'{hline}/{sline}/veldm',data=veldm)
                                                f.create_dataset(f'{hline}/{sline}/velgas',data=velg)
                                                
                                                # Write mass

                                                f.create_dataset(f'{hline}/{sline}/massstar',data=massstar)
                                                f.create_dataset(f'{hline}/{sline}/massdm',data=massdm)
                                                f.create_dataset(f'{hline}/{sline}/massgas',data=massg)


                                                # write temperature
                                                f.create_dataset(f'{hline}/{sline}/tgas',data=tempg)

                                                # write metallicities
                                                f.create_dataset(f'{hline}/{sline}/zgas',data=metalg)
                                                f.create_dataset(f'{hline}/{sline}/fegas',data=feg)
                                                f.create_dataset(f'{hline}/{sline}/ogas',data=og)
                                                f.create_dataset(f'{hline}/{sline}/hgas',data=hg)
                                                f.create_dataset(f'{hline}/{sline}/zstar',data=zstar)

                                                
                                                # get cluster group

                                                clus = f[f'/{hline}/']
                                                for par in ['nsub','nstar','nsink','ngas','mtot','mdm',\
                                                    'mgas','msink','mstar','pos','vel']:

                                                    clus.attrs[par] = thalo[par]

                                                gal = f[f'/{hline}/{sline}']
                                                for par in ['nstar','nsink','ngas','mtot','mdm',\
                                                    'mgas','msink','mstar','pos','vel']:

                                                    gal.attrs[par] = tsub[par]

                                        
                                                sline = sline + 1
                                        
                                        print(f"snap={snapno} ,Cluster found!! no.={kkk} ID={hline},{clusters['HostHaloID'].iloc[kkk]}")

                                                
                                        kkk=kkk+1

                                    
                                                
                        else:   
                                # Skip over uninteresting halos
                                
                                # Read fof halo
                                data1 = unpack('@6i11d',fof)
                                            
                                thalo={'nsub':data1[0],'ndm': data1[1],'nstar': data1[2],'nsink':data1[3],'ngas':data1[4],'npall':data1[5],'mtot':data1[6],\
                                            'mdm':data1[7],'mgas':data1[8],'msink':data1[9],'mstar':data1[10],'pos':data1[11:14],'vel':data1[14:17]}
                                            
                                #print(thalo)         
                                # Read subhaloes
                                for i in range(thalo['nsub']):
                                        
                                        data2 = unpack('@6i11d',file.read(112))
                                        tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],\
                                            'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}

                                        # Read dm particles
                                        for j in range(tsub['ndm']):               
                                            file.read(128)
                                            
                                        # Read gas particles
                                        if tsub['ngas']>0:
                                            for j in range(tsub['ngas']):
                                                file.read(128)

                                        if tsub['nsink']>0:
                                            # Read sink particles
                                            for j in range(tsub['nsink']):
                                                file.read(168)
                                                
                                        if tsub['nstar']>0:
                                            for j in range(tsub['nstar']):
                                                file.read(128)
                                    

                                        sline = sline + 1

                        hline = hline + 1
                                        


#Checking if all the clusters were found and stored in the hdf5 file
out = f"{output}clusters{snapno}.hdf5"
if os.path.exists(out):
            print(out)
            with h5py.File(out, 'r') as file:
                
                keys = list(file.keys())
                keys = [int(clus) for clus in keys]
                clusfile = f"{snapno}.dat"

                clusters = np.genfromtxt(clusfile,skip_header=1,dtype=int)

                clusters = clusters.tolist()

                if isinstance(clusters, int):
                    clusters = [clusters]
                # check the clusters that have been processed already
                diff2 = set(clusters) - set(keys)
                
                assert len(diff2) == 0, "The following clusters were not found in the hdf5 file: {}".format(diff2)
                    
