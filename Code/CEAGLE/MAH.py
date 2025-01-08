import numpy as np
import pandas as pd 
from scipy.io import FortranFile
import matplotlib.pyplot as plt 
from os.path import isfile, join
from os import listdir
import csv
from parallelbar import progress_map

# plt.style.use('/home/ankitsingh/HR5_AGN/paper_style.mplstyle')

import configparser



def read_snapfile(snapno):


    file = f'/scratch/ankitsingh/Galaxy_catalogs/galaxy_catalogue_{snapno}.parquet'


    snapi = pd.read_parquet(file)



    return snapi


def Traverse(clusID=0):

        cluster = snaplast.loc[snaplast['HostHaloID'].to_numpy()==clusID]

        # get total mass of first galaxy
                
        idin=cluster.loc[cluster['Mstar(Msun)'].to_numpy()==cluster['Mstar(Msun)'].max(),'ID'].values[0]
                    
        fm = cluster['HostMtot(Msun)'].to_numpy()[0]
                    
        # idin is the ID of the main galaxy in the cluster        

        with open(f'{output}/{clusID}_MAH.csv','w') as fout:

            fout.write('time,snap,host_flag,flag_prog,HostHaloID,MainGalID,ClusMass(Msun),Massfraction\n')
            
            # get the number of snapshot
            snapnum = [file.split('_')[3].split('.parquet')[0] for file in hr5files]
            
            
            iddec=idin
            i=0

            l=0
            thalf=0
            while i < (len(snapnum)-1): 
                
                time_simu = hr5outs['LBT'].to_numpy()[hr5outs['Snapshot'].to_numpy()==snapnum[i]][0]
                
                
                snapi = snapnum[i]
                
                #red = redshift.loc[redshift['snap']==int(snap),'z']
                # read the files
                
                sdati = read_snapfile(snapi)

                #print(snapnum[i],iddec,clusID)
                # flag for the progenitor

                try:
                    flagp = sdati.loc[sdati['ID'].to_numpy()==iddec,'flag_prog'].values[0]
                except IndexError:
                    break

                # Mass and ID of the Host FOF Halo of the main galaxy in the descendent snapshot
                FOFmassi=sdati['HostMtot(Msun)'].to_numpy()[sdati['ID'].to_numpy()==iddec][0]
                FOFIDi=sdati['HostHaloID'].to_numpy()[sdati['ID'].to_numpy()==iddec][0]
                central = sdati.loc[sdati['ID'].to_numpy()==iddec,'host_flag'].values[0]
                
                
                if FOFmassi/fm<0.5:
                            if l==0:
                                thalf = time_simu
                                l=l+1

                #print(sdati.loc[sdati['ID'].to_numpy()==iddec][['flag_prog', 'ID_prog_gen','step_descen_gen','host_flag']])

                # We are traversing the tree bottom up, checking the progenetor in the previous 
                # snapshot and calculating the Mass of the main galaxy
                fout.write(f'{time_simu},{snapi},{central},{flagp},{FOFIDi},{iddec},{FOFmassi:.3e},{FOFmassi/fm:.3f}\n')
                # print(f'{time_simu},{snapi},{central},{flagp},{FOFIDi},{iddec},{FOFmassi:.3e},{FOFmassi/fm:.3f}\n')
            
                if flagp!=3:

                    # snapi_p1 = snapnum[i+1]
                    # sdat_p1 = read_snapfile(snapi_p1)
                    
                    idprog = sdati['ID_prog'].to_numpy()[sdati['ID'].to_numpy()==iddec][0]
                    # All galaxies with descendent ID equal to the main galaxy of the cluster
                    #snapproj = idprog#sdat_p1['ID'].to_numpy()[sdat_p1['ID_descen'].to_numpy() == iddec]
                    
                    
                    iddec = idprog
                        
                    i=i+1
                    

                else:
                    # Case where geneuine progenetor is at different snapshot and maybe a satellite


                    # sdat_p1 = read_snapfile(snapnum[i])

                    idprog = sdati['ID_prog_gen'].to_numpy()[sdati['ID'].to_numpy()==iddec][0]
                    central = sdati.loc[sdati['ID'].to_numpy()==iddec,'host_flag'].values[0]

                    stepi_p1 = sdati['step_prog_gen'].to_numpy()[sdati['ID'].to_numpy()==iddec][0]

                    if stepi_p1==0:
                        break
                    else:
                        i = snapnum.index(str(stepi_p1).zfill(3))

                    iddec = idprog
                
                
                
        halfred = hr5outs.loc[hr5outs['LBT']==thalf,'Redshift'].values

        
        if halfred.size == 0:
            halfred=[0]

        return halfred  




# load up the parameter file
parser = configparser.ConfigParser()
parser.read('../params.ini')


Fofd = parser.get('Paths','Fofdir')
output = parser.get('Paths','outdir')
clusfile=parser.get('Paths','clusfile')
hr5cat=parser.get('Paths','galaxycats')
snapno = int(parser.get('Setting','snapno'))
hr5time=parser.get('Paths','hr5outs')

#get list of catalogs

hr5files = sorted([join(hr5cat,f) for f in listdir(hr5cat) if (isfile(join(hr5cat,f))& f.endswith('.parquet')) ],reverse=True)


clusters = pd.read_csv(clusfile)


snaplast = read_snapfile(snapno)


# Get snapshot corresponding to redshift
hr5outs = pd.read_csv(hr5time,dtype={'Snapshot':str,'Redshift':float,\
                    'LBT':float,'dx':float})



hr5files = sorted([join(hr5cat,f) for f in listdir(hr5cat) if (isfile(join(hr5cat,f))& f.endswith('.parquet')) ],reverse=True)




    
with open(f"{output}halfmass.txt", 'w') as hfile:
        
        hfile.write(f"HostHaloID,Redshift\n")
        cluslist = clusters['HostHaloID'].to_list()
        
        hred = progress_map(Traverse,cluslist, chunk_size=1,n_cpu=8)
        for cl, rr in zip(cluslist, hred):

            hfile.write(f"{cl},{rr}\n")
yo