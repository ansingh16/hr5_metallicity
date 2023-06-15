import pandas as pd 
from os.path import isfile, join
from os import listdir


# Function for reading the snapfiles
def read_snapfile(snapno):


    file = f'/scratch/ankitsingh/Galaxy_catalogs/galaxy_catalogue_{snapno}.parquet'


    snapi = pd.read_parquet(file)



    return snapi


# Function for traversing the merger tree for a cluster
def Traverse(output,idin=0,clusID=0):


        # Get snapshot corresponding to redshift
        hr5outs = pd.read_csv('./Time_data.dat',dtype={'Snapshot':str,'Redshift':float,\
                            'LBT':float,'dx':float})


        hr5cat = '/scratch/ankitsingh/Galaxy_catalogs/'
        
        hr5files = sorted([join(hr5cat,f) for f in listdir(hr5cat) if (isfile(join(hr5cat,f))& f.endswith('.csv')) ],reverse=True)

        # get the number of snapshot
        snapnum = [file.split('_')[3].split('.csv')[0] for file in hr5files]
        
        iddec=idin
        i=0
        l=0
        thalf=0
        while i < (len(snapnum)-1): 
                
                time = hr5outs['LBT'].to_numpy()[hr5outs['Snapshot'].to_numpy()==snapnum[i]][0]
                
                
                snapi = snapnum[i]
                snapi_p1 = snapnum[i+1]
                
                #red = redshift.loc[redshift['snap']==int(snap),'z']
                # read the files
                
                sdati = read_snapfile(snapi)

                sdat_p1 = read_snapfile(snapi_p1)

                snapdec = sdati['ID_prog'].to_numpy()[sdati['ID'].to_numpy()==iddec]

                
                snapproj = sdat_p1['ID'].to_numpy()[sdat_p1['ID_descen'].to_numpy() == iddec]

                
                if snapproj.shape[0]>=1:

                    # This means more than one progenetor of the galaxy
                    # We have to select minor and major mergers from these
                    # print(idin,iddec,snapdec['ID_prog'].values[0])
                    # choose main proginetor
                    iddec = snapdec[0]
                    
                    
                    
                    i=i+1

                else:

                    break 