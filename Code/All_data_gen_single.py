import numpy as np
from struct import unpack
import pandas as pd 
import h5py
from scipy.io import FortranFile
# plt.style.use('/home/ankitsingh/HR5_AGN/paper_style.mplstyle')
import configparser
# from parallelbar import progress_map
import os 
import sys 
import tqdm

import concurrent.futures

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('/home/jwyoo/WOC_SIDM/hr5_metallicity/params.ini')

# clusfile = parser.get('Paths','clusfile')
Fofd = parser.get('Paths','Fofdir')
outdir = parser.get('Paths','outdir')


def read_icl(fout,hline,fof_icl,file_back):

    #print(hline)
    # Read fof halo
    data2 = unpack('@6i11d',fof_icl)
                                        
    tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}
                    
    #print(tsub)
    posdm = np.empty((tsub['ndm'],3))
    massdm = np.empty(tsub['ndm'])
    veldm = np.empty((tsub['ndm'],3))
    
    # Define the binary format for the structure with 3 doubles for pos, 3 doubles for vel, 1 double for mass, and 1 long long integer for id.
    data_format = '@6d1d1q'  # 3 doubles for pos, 3 doubles for vel, 1 double for mass, 1 long long integer for id
    for j in range(tsub['ndm']):

                                                    
        data3 = unpack(data_format, file_back.read(64))  # Read 64 bytes, since 6d + 1d + 1q = 8 elements, which equals 64 bytes

        # Map the unpacked data to the structure
        tdm = {'pos': data3[0:3],'vel': data3[3:6],'mass': data3[6],'id': data3[7]}
        posdm[j,:] = np.array(tdm['pos']) 
        veldm[j,:] = np.array(tdm['vel']) 
        massdm[j] = tdm['mass']
                                                
    

    posg = np.empty((tsub['ngas'],3))
    massg = np.empty(tsub['ngas'])
    velg = np.empty((tsub['ngas'],3))
                            
    data_format_gas = '@3d3f1f'
    if tsub['ngas']>0:
                            # Read gas particles
        for j in range(tsub['ngas']):
            #print(tsub['ngas'],j) 
            data4 = unpack(data_format_gas, file_back.read(40))  # Read 32 bytes, since 3d + 3f + 1d = 32 bytes
            tgas = {'pos': data4[0:3],'vel': data4[3:6],'mass': data4[6]}
            posg[j,:] = np.array(tgas['pos']) 
            velg[j,:] = np.array(tgas['vel']) 
            massg[j] = tgas['mass']
          
    
    
    posstar = np.empty((tsub['nstar'],3))
    massstar = np.empty(tsub['nstar'])
    velstar = np.empty((tsub['nstar'],3))

    if tsub['nstar']>0:
                                                    
        # Read star particles
        for j in range(tsub['nstar']):
                                    
            data5 = unpack(data_format, file_back.read(64))  # Read 64 bytes, since 6d + 1d + 1q = 8 elements, which equals 64 bytes

            # Map the unpacked data to the structure
            tstar = {'pos': data5[0:3],'vel': data5[3:6],'mass': data5[6],'id': data5[7]}

                                                        
            posstar[j,:] = np.array(tstar['pos']) 
            velstar[j,:] = np.array(tstar['vel']) 
            massstar[j] = tstar['mass']
                                    
    # check if subhalo is exists
    if f'{hline}/ICL/' in fout:

        # delete the old subhalo
        del fout[f'{hline}/ICL']

    # write positions of stars in subhalo
    fout.create_dataset(f'{hline}/ICL/posstar',data=posstar)
    fout.create_dataset(f'{hline}/ICL/posdm',data=posdm)
    fout.create_dataset(f'{hline}/ICL/posgas',data=posg)

    # write positions of stars in subhalo
    fout.create_dataset(f'{hline}/ICL/velstar',data=velstar)
    fout.create_dataset(f'{hline}/ICL/veldm',data=veldm)
    fout.create_dataset(f'{hline}/ICL/velgas',data=velg)
                                
    # Write mass

    fout.create_dataset(f'{hline}/ICL/massstar',data=massstar)
    fout.create_dataset(f'{hline}/ICL/massdm',data=massdm)
    fout.create_dataset(f'{hline}/ICL/massgas',data=massg)



    # write metallicities

    # write attributes   
    icl = fout[f'/{hline}/ICL/']
    for par in ['ndm','ngas','nsink','nstar','vel']:
                            
        icl.attrs[par] = tsub[par]
    for par in ['mtot','mdm','mgas','msink','mstar','pos']:
        try:
            icl.attrs[par] = tsub[par] 
        except:
            icl.attrs[par] = np.array(tsub[par])

   



def read_fof(fout,sline,hline,fof,file_fof):

    # Read fof halo
    data1 = unpack('@6i11d',fof)
                                            
    thalo={'nsub':data1[0],'ndm': data1[1],'nstar': data1[2],'nsink':data1[3],'ngas':data1[4],'npall':data1[5],'mtot':data1[6],'mdm':data1[7],'mgas':data1[8],'msink':data1[9],'mstar':data1[10],'pos':data1[11:14],'vel':data1[14:17]}
    
    #print('Halo: ',thalo)
                                            
    # Read subhaloes
    for i in tqdm.tqdm(range(thalo['nsub'])):
        data2 = unpack('@6i11d',file_fof.read(112))
        tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}

        #print(f"In halo {i} subhalo: {tsub}\n")
        tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}

        #print(tsub)
        posdm = np.empty((tsub['ndm'],3))
        massdm = np.empty(tsub['ndm'])
        veldm = np.empty((tsub['ndm'],3))

        # Define the binary format for the structure with 3 doubles for pos, 3 doubles for vel, 1 double for mass, and 1 long long integer for id.
        data_format = '@6d1d1q'  # 3 doubles for pos, 3 doubles for vel, 1 double for mass, 1 long long integer for id
        for j in range(tsub['ndm']):


            data3 = unpack(data_format, file_fof.read(64))  # Read 64 bytes, since 6d + 1d + 1q = 8 elements, which equals 64 bytes

            # Map the unpacked data to the structure
            tdm = {'pos': data3[0:3],'vel': data3[3:6],'mass': data3[6],'id': data3[7]}
            posdm[j,:] = np.array(tdm['pos'])
            veldm[j,:] = np.array(tdm['vel'])
            massdm[j] = tdm['mass']



        posg = np.empty((tsub['ngas'],3))
        massg = np.empty(tsub['ngas'])
        velg = np.empty((tsub['ngas'],3))

        data_format_gas = '@3d3f1f'
        if tsub['ngas']>0:
                            # Read gas particles
            for j in range(tsub['ngas']):
                #print(tsub['ngas'],j) 
                data4 = unpack(data_format_gas, file_fof.read(40))  # Read 32 bytes, since 3d + 3f + 1d = 32 bytes
                tgas = {'pos': data4[0:3],'vel': data4[3:6],'mass': data4[6]}
                posg[j,:] = np.array(tgas['pos'])
                velg[j,:] = np.array(tgas['vel'])
                massg[j] = tgas['mass']



        posstar = np.empty((tsub['nstar'],3))
        massstar = np.empty(tsub['nstar'])
        velstar = np.empty((tsub['nstar'],3))

        if tsub['nstar']>0:

            # Read star particles
            for j in range(tsub['nstar']):

                data5 = unpack(data_format, file_fof.read(64))  # Read 64 bytes, since 6d + 1d + 1q = 8 elements, which equals 64 bytes

                # Map the unpacked data to the structure
                tstar = {'pos': data5[0:3],'vel': data5[3:6],'mass': data5[6],'id': data5[7]}


                posstar[j,:] = np.array(tstar['pos'])
                velstar[j,:] = np.array(tstar['vel'])
                massstar[j] = tstar['mass']

        # check if subhalo is exists
        if f'{hline}/{sline}' in fout:

            # delete the old subhalo
            del fout[f'{hline}/{sline}']

                         

        # write positions of stars in subhalo
        fout.create_dataset(f'{hline}/{sline}/posstar',data=posstar)
        fout.create_dataset(f'{hline}/{sline}/posdm',data=posdm)
        fout.create_dataset(f'{hline}/{sline}/posgas',data=posg)

        # write positions of stars in subhalo
        fout.create_dataset(f'{hline}/{sline}/velstar',data=velstar)
        fout.create_dataset(f'{hline}/{sline}/veldm',data=veldm)
        fout.create_dataset(f'{hline}/{sline}/velgas',data=velg)
                                                    
        # Write mass

        fout.create_dataset(f'{hline}/{sline}/massstar',data=massstar)
        fout.create_dataset(f'{hline}/{sline}/massdm',data=massdm)
        fout.create_dataset(f'{hline}/{sline}/massgas',data=massg)



                                                    
        # get cluster group

        clus = fout[f'/{hline}/']
        for par in ['nsub','nstar','nsink','vel']:

            clus.attrs[par] = thalo[par]
                                                    
        for par in ['mtot','mdm','mgas','msink','mstar','pos']:

            try:
                clus.attrs[par] = thalo[par]
            except:
                clus.attrs[par] = np.array(thalo[par])

        gal = fout[f'/{hline}/{sline}']
        for par in ['nstar','nsink','ngas','vel']:
            gal.attrs[par] = tsub[par]

        for par in ['mtot','mdm','mgas','msink','mstar','pos']:

            try:
                gal.attrs[par] = tsub[par]
            except:
                gal.attrs[par] = np.array(tsub[par])
                                                        
        sline = sline + 1

    return sline 



def Make_hdf5(snapno,clusters):



    # clusfile = f"{snapfiles}{snapno}.dat"
    # clusters = pd.read_csv(clusfile,usecols=['HostHaloID'])
    # clusters.sort_values('HostHaloID', inplace=True,ignore_index=True)

    with h5py.File(f"{outdir}/clusters{snapno}.hdf5", "a") as fout:
                
                if 'status' in fout:
                    del fout['status']
                    fout.create_dataset('status',dtype=np.int32,data=0)
                else:
                    fout.create_dataset('status',dtype=np.int32,data=0)
                    
                    
                # Open Galaxy find data files
                with open(f'{Fofd}FoF.{snapno:0>5}/GALFIND.DATA.{snapno:0>5}', mode='rb') as file_fof: # b -> binary

                    with open(f'{Fofd}/FoF.{snapno:0>5}/background_ptl.{snapno:0>5}', mode='rb') as file_back: # b -> binary

                        hline=0
                        sline=0
                        kkk=0
                        while True:
                        
                            fof = file_fof.read(112)

                            fof_icl = file_back.read(112)

                            if not fof or (kkk>clusters.shape[0]-1):
                                print(f"Done!! snap={snapno},total clus={kkk},endded={hline}")
                                break
                            else:
                                # print(hline,sline,kkk,clusters.loc[kkk,'HostHaloID'])
                                
                                # Check if the cluster is of interest
                                if clusters[kkk]==hline:
                                
                                                print(f"found cluster: {hline} in {snapno}")
                                                print("Processing ICL") 
                                                read_icl(fout,hline,fof_icl,file_back)
                                                print("Processing FOF")
                                                sline = read_fof(fout,sline,hline,fof,file_fof)
                                            
                                                kkk=kkk+1

                                            
                                                        
                                else:   
                                        #print("me here1",sline,hline)
                                        skip_icl(fof_icl,file_back)
                                        # print("me here2",sline,hline)
                                        sline = skip_fof(sline,fof,file_fof)
                                        # print("me here3",sline,hline)
                                
                                hline = hline + 1 
                
                del fout['status']
                fout.create_dataset('status',dtype=np.int32,data=1)
                
                        


for snapno in range(1,31):
    print(f"Processing Snapshot No. {snapno}")
    Make_hdf5(snapno,np.array([0]))

