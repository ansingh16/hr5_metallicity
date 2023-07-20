import numpy as np
cimport numpy as np
from struct import unpack
import pandas as pd 
import h5py
from scipy.io import FortranFile
# plt.style.use('/home/ankitsingh/HR5_AGN/paper_style.mplstyle')
import configparser
# from parallelbar import progress_map
import os 
# import glob
# import random
import sys 
import glob
import tqdm 

from cython cimport boundscheck, wraparound
from libc.stdlib cimport malloc, free
cimport cython
from cython cimport boundscheck, wraparound
from libc.stdlib cimport malloc, free


# load up the parameter file
parser = configparser.ConfigParser()
parser.read('/home/ankitsingh/hr5_metalicity/params.ini')

# clusfile = parser.get('Paths','clusfile')
Fofd = parser.get('Paths','Fofdir')
outdir = parser.get('Paths','outdir')
snapfiles = parser.get('Paths','snapfiles')

h0=float(parser.get('Setting','h0'))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void read_icl( fout, int hline, bytes fof_icl, file_back):
    cdef bytes data2, data3, data4, data5, data6
    cdef int j
    cdef dict tsub, tdm, tgas, tsink, tstar
    cdef np.ndarray posdm, massdm, veldm, posg, massg, velg, tempg, metalg, feg, hg, og, posstar, massstar, velstar, zstar

    # print(hline)
    # Read fof halo
    data2 = unpack('@6i11d',fof_icl)
                                        
    tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}
                    
    #print(tsub)
    posdm = np.empty((tsub['ndm'],3))
    massdm = np.empty(tsub['ndm'])
    veldm = np.empty((tsub['ndm'],3))

    for j in range(tsub['ndm']):

                                                    
        data3 = unpack('@13d1q1d1i1f',file_back.read(128))
        tdm={'pos':data3[0:3],'vel':data3[3:6],'mass':data3[6],'dum0':data3[7],'tp':data3[8],'zp':data3[9],'mass0':data3[10],'tpp':data3[11],'indtab':data3[12],'id':data3[13],'potential':data3[14],'level':data3[15],'dum1':data3[16]}      

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
                                                        
            data4 = unpack('@4d4f1d5f1i2f1q4d',file_back.read(128))
            tgas={'pos':data4[0:3],'dx':data4[3],'vel':data4[4:7],'dum0':data4[7],'density':data4[8],'temp':data4[9],'metal':data4[10],'fe':data4[11],'h':data4[12],'o':data4[13],'level':data4[14],'mass':data4[15],'dum1':data4[16],'id':data4[17],'potential':data4[18],'f':data4[19:22]}

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
            data5 = unpack('@20d2i',file_back.read(168))
            tsink={'pos':data5[0:3],'vel':data5[3:6],'mass':data5[6],'tbirth':data5[7],'angm':data5[8:11],'ang':data5[11:14],'dmsmbh':data5[14:17],'esave':data5[17],'smag':data5[18],'eps':data5[19],'id':data5[20],'dum0':data5[21]}


    posstar = np.empty((tsub['nstar'],3))
    massstar = np.empty(tsub['nstar'])
    velstar = np.empty((tsub['nstar'],3))
    zstar = np.empty(tsub['nstar'])

    if tsub['nstar']>0:
                                                    
        # Read star particles
        for j in range(tsub['nstar']):
                                    
            data6 = unpack('@13d1q1d1i1f',file_back.read(128))
                                    
            tstar={'pos':data6[0:3],'vel':data6[3:6],'mass':data6[6],'dum0':data6[7],'tp':data6[8],'zp':data6[9],'mass0':data6[10],'tpp':data6[11],'indtab':data6[12],'id':data6[13],'potential':data6[14],'level':data6[15],'dum1':data6[16]}

                                                        
            posstar[j,:] = np.array(tstar['pos']) 
            velstar[j,:] = np.array(tstar['vel']) 
            massstar[j] = tstar['mass']
            zstar[j] = tstar['zp']
                                    


    # write positions of stars in subhalo
    fout.create_dataset(f'{hline}/ICL/posstar',data=posstar/h0)
    fout.create_dataset(f'{hline}/ICL/posdm',data=posdm/h0)
    fout.create_dataset(f'{hline}/ICL/posgas',data=posg/h0)

    # write positions of stars in subhalo
    fout.create_dataset(f'{hline}/ICL/velstar',data=velstar)
    fout.create_dataset(f'{hline}/ICL/veldm',data=veldm)
    fout.create_dataset(f'{hline}/ICL/velgas',data=velg)
                                
    # Write mass

    fout.create_dataset(f'{hline}/ICL/massstar',data=massstar/h0)
    fout.create_dataset(f'{hline}/ICL/massdm',data=massdm/h0)
    fout.create_dataset(f'{hline}/ICL/massgas',data=massg/h0)


    # write temperature
    fout.create_dataset(f'{hline}/ICL/tgas',data=tempg)

    # write metallicities
    fout.create_dataset(f'{hline}/ICL/zgas',data=metalg)
    fout.create_dataset(f'{hline}/ICL/fegas',data=feg)
    fout.create_dataset(f'{hline}/ICL/ogas',data=og)
    fout.create_dataset(f'{hline}/ICL/hgas',data=hg)
    fout.create_dataset(f'{hline}/ICL/zstar',data=zstar)

    # write attributes   
    icl = fout[f'/{hline}/ICL/']
    for par in ['ndm','ngas','nsink','nstar','vel']:
                            
        icl.attrs[par] = tsub[par]
    for par in ['mtot','mdm','mgas','msink','mstar','pos']:
        try:
            icl.attrs[par] = tsub[par]/h0 
        except:
            icl.attrs[par] = np.array(tsub[par])/h0 

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void skip_icl(bytes fof_icl, file_back):
    cdef bytes data2
    cdef dict tsub
    cdef int j
    # Read fof halo
    data2 = unpack('@6i11d',fof_icl)
                                        
    tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}
                    
                        
    # Read dm particles
    for j in range(tsub['ndm']):               
        
        file_back.read(128)

    if tsub['ngas']>0:
        # Read gas particles
        for j in range(tsub['ngas']):
            file_back.read(128)

    if tsub['nsink']>0:
        # Read sink particles
        for j in range(tsub['nsink']):
            file_back.read(168)
                                        
    if tsub['nstar']>0:
        for j in range(tsub['nstar']):
            file_back.read(128)
        


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int read_fof( fout, int sline, int hline, bytes fof,  file_fof):

    cdef int i, j
    cdef bytes data1, data2, data3, data4, data5, data6
    cdef dict thalo, tsub, tdm, tgas, tsink, tstar
    cdef np.ndarray posdm, massdm, veldm, posg, massg, velg, tempg, metalg, feg, hg, og, posstar, massstar, velstar, zstar

    # Read fof halo
    data1 = unpack('@6i11d', fof)

    thalo = {'nsub': data1[0], 'ndm': data1[1], 'nstar': data1[2], 'nsink': data1[3], 'ngas': data1[4], 'npall': data1[5], 'mtot': data1[6], 'mdm': data1[7], 'mgas': data1[8], 'msink': data1[9], 'mstar': data1[10], 'pos': data1[11:14], 'vel': data1[14:17]}

    # Read subhaloes
    for i in range(thalo['nsub']):
        data2 = unpack('@6i11d', file_fof.read(112))
        tsub = {'ndm': data2[0], 'ngas': data2[1], 'nsink': data2[2], 'nstar': data2[3], 'npall': data2[4], 'dum': data2[5], 'mtot': data2[6], 'mdm': data2[7], 'mgas': data2[8], 'msink': data2[9], 'mstar': data2[10], 'pos': data2[11:14], 'vel': data2[14:17]}

        posdm = np.empty((tsub['ndm'], 3))
        massdm = np.empty(tsub['ndm'])
        veldm = np.empty((tsub['ndm'], 3))

        # Read dm particles
        for j in range(tsub['ndm']):
            data3 = unpack('@13d1q1d1i1f', file_fof.read(128))
            tdm = {'pos': data3[0:3], 'vel': data3[3:6], 'mass': data3[6], 'dum0': data3[7], 'tp': data3[8], 'zp': data3[9], 'mass0': data3[10], 'tpp': data3[11], 'indtab': data3[12], 'id': data3[13], 'potential': data3[14], 'level': data3[15], 'dum1': data3[16]}
            posdm[j, :] = np.array(tdm['pos'])
            veldm[j, :] = np.array(tdm['vel'])
            massdm[j] = tdm['mass']

        posg = np.empty((tsub['ngas'], 3))
        massg = np.empty(tsub['ngas'])
        velg = np.empty((tsub['ngas'], 3))
        tempg = np.empty(tsub['ngas'])
        metalg = np.empty(tsub['ngas'])
        feg = np.empty(tsub['ngas'])
        hg = np.empty(tsub['ngas'])
        og = np.empty(tsub['ngas'])

        if tsub['ngas'] > 0:
            # Read gas particles
            for j in range(tsub['ngas']):
                data4 = unpack('@4d4f1d5f1i2f1q4d', file_fof.read(128))
                tgas = {'pos': data4[0:3], 'dx': data4[3], 'vel': data4[4:7], 'dum0': data4[7], 'density': data4[8], 'temp': data4[9], 'metal': data4[10], 'fe': data4[11], 'h': data4[12], 'o': data4[13], 'level': data4[14], 'mass': data4[15], 'dum1': data4[16], 'id': data4[17], 'potential': data4[18], 'f': data4[19:22]}
                posg[j, :] = np.array(tgas['pos'])
                velg[j, :] = np.array(tgas['vel'])
                massg[j] = tgas['mass']
                tempg[j] = tgas['temp'] / tgas['density']
                metalg[j] = tgas['metal']
                feg[j] = tgas['fe']
                og[j] = tgas['o']
                hg[j] = tgas['h']

        if tsub['nsink'] > 0:
            # Read sink particles
            for j in range(tsub['nsink']):
                data5 = unpack('@20d2i', file_fof.read(168))
                tsink = {'pos': data5[0:3], 'vel': data5[3:6], 'mass': data5[6], 'tbirth': data5[7], 'angm': data5[8:11], 'ang': data5[11:14], 'dmsmbh': data5[14:17], 'esave': data5[17], 'smag': data5[18], 'eps': data5[19], 'id': data5[20], 'dum0': data5[21]}

        posstar = np.empty((tsub['nstar'], 3))
        massstar = np.empty(tsub['nstar'])
        velstar = np.empty((tsub['nstar'], 3))
        zstar = np.empty(tsub['nstar'])

        if tsub['nstar'] > 0:

            # Read star particles
            for j in range(tsub['nstar']):
                data6 = unpack('@13d1q1d1i1f',file_fof.read(128))
                tstar={'pos':data6[0:3],'vel':data6[3:6],'mass':data6[6],'dum0':data6[7],'tp':data6[8],'zp':data6[9],'mass0':data6[10],'tpp':data6[11],'indtab':data6[12],'id':data6[13],'potential':data6[14],'level':data6[15],'dum1':data6[16]}

                                                            
                posstar[j,:] = np.array(tstar['pos']) 
                velstar[j,:] = np.array(tstar['vel']) 
                massstar[j] = tstar['mass']
                zstar[j] = tstar['zp']
                                        
                                                    
                                                    

        # write positions of stars in subhalo
        fout.create_dataset(f'{hline}/{sline}/posstar',data=posstar/h0)
        fout.create_dataset(f'{hline}/{sline}/posdm',data=posdm/h0)
        fout.create_dataset(f'{hline}/{sline}/posgas',data=posg/h0)

        # write positions of stars in subhalo
        fout.create_dataset(f'{hline}/{sline}/velstar',data=velstar)
        fout.create_dataset(f'{hline}/{sline}/veldm',data=veldm)
        fout.create_dataset(f'{hline}/{sline}/velgas',data=velg)
                                                    
        # Write mass

        fout.create_dataset(f'{hline}/{sline}/massstar',data=massstar/h0)
        fout.create_dataset(f'{hline}/{sline}/massdm',data=massdm/h0)
        fout.create_dataset(f'{hline}/{sline}/massgas',data=massg/h0)


        # write temperature
        fout.create_dataset(f'{hline}/{sline}/tgas',data=tempg)

        # write metallicities
        fout.create_dataset(f'{hline}/{sline}/zgas',data=metalg)
        fout.create_dataset(f'{hline}/{sline}/fegas',data=feg)
        fout.create_dataset(f'{hline}/{sline}/ogas',data=og)
        fout.create_dataset(f'{hline}/{sline}/hgas',data=hg)
        fout.create_dataset(f'{hline}/{sline}/zstar',data=zstar)

                                                    
        # get cluster group

        clus = fout[f'/{hline}/']
        for par in ['nsub','nstar','nsink','vel']:

            clus.attrs[par] = thalo[par]
                                                    
        for par in ['mtot','mdm','mgas','msink','mstar','pos']:

            try:
                clus.attrs[par] = thalo[par]/h0
            except:
                
                clus.attrs[par] = np.array(thalo[par])/h0

        gal = fout[f'/{hline}/{sline}']
        for par in ['nstar','nsink','ngas','vel']:
            gal.attrs[par] = tsub[par]

        for par in ['mtot','mdm','mgas','msink','mstar','pos']:

            try:
                gal.attrs[par] = tsub[par]/h0
            except:
                gal.attrs[par] = np.array(tsub[par])/h0
                                                        
        sline = sline + 1

    return sline 



cdef int skip_fof(int sline, bytes fof, file_fof):
    cdef int i, j
    cdef bytes data1, data2
    cdef dict thalo, tsub

    # Skip over uninteresting halos

    # Read fof halo
    data1 = unpack('@6i11d', fof)

    thalo = {'nsub': data1[0], 'ndm': data1[1], 'nstar': data1[2], 'nsink': data1[3], 'ngas': data1[4], 'npall': data1[5], 'mtot': data1[6], 'mdm': data1[7], 'mgas': data1[8], 'msink': data1[9], 'mstar': data1[10], 'pos': data1[11:14], 'vel': data1[14:17]}

    # Read subhaloes
    for i in range(thalo['nsub']):
        data2 = unpack('@6i11d', file_fof.read(112))
        tsub = {'ndm': data2[0], 'ngas': data2[1], 'nsink': data2[2], 'nstar': data2[3], 'npall': data2[4], 'dum': data2[5], 'mtot': data2[6], 'mdm': data2[7], 'mgas': data2[8], 'msink': data2[9], 'mstar': data2[10], 'pos': data2[11:14], 'vel': data2[14:17]}

        # Read dm particles
        for j in range(tsub['ndm']):
            file_fof.read(128)

        # Read gas particles
        if tsub['ngas'] > 0:
            for j in range(tsub['ngas']):
                file_fof.read(128)

        if tsub['nsink'] > 0:
            # Read sink particles
            for j in range(tsub['nsink']):
                file_fof.read(168)

        if tsub['nstar'] > 0:
            for j in range(tsub['nstar']):
                file_fof.read(128)

        sline = sline + 1

    return sline
