import numpy as np
import pandas as pd 
import h5py
import yt 
from yt.units import Msun, Mpc
import tqdm 
import matplotlib.pyplot as plt 
#plt.style.use('/home/ankitsingh/hr5_agn/paper_style.mplstyle')

# low information printing for yt
# avoiding unnecessary information can save execution time 
#yt.set_log_level(40)


import configparser

# get all snapshots
import os 
import glob

import configparser





# load up the parameter file
parser = configparser.ConfigParser()
parser.read('params.ini')


# setting paths containg HR5 directories
#Fofd = parser.get('Paths','Fofdir')
output = parser.get('Paths','outdir')
r200 = parser.get('Paths','r200')

r200 = pd.read_csv(r200,skiprows=1,names=['snapshot' ,'groupnr' ,'subgroupnr' ,'CoP_x' ,'CoP_y', 'CoP_z' ,'R200', 'vx_CoP' ,'vy_CoP' ,'vz_CoP' \
                    ,'redshift' ,'lookbacktime'],sep=' ')[['snapshot','R200']]



# Function for saving png corresponding to each cluster
# at each snapshot
def cluster_projection(cluslast,clusID):
    
    # Function for making projection plots of the clusters
    # It reads the data from the hdf5 files and saves the plots
    # in the plotdir path defined in params.ini

    fig,ax = plt.subplots(1,1)
    ax.set_xlabel(r'Y (cMpc)')
    ax.set_ylabel(r'Z (cMpc)')
    
    # Get the position of unbound stars
    pos = f[f'/{clusID}/ICL/posstar']

    ax.scatter(pos[:,1],pos[:,2],s=0.2,color='r',alpha=0.4)

    for gal in f[f'/{clusID}/'].keys():

        # stars in galaxies are to be plotted in blue
        pos =  f[f'/{clusID}/{gal}/posstar']
        
        if gal!='ICL':
            ax.scatter(pos[:,1],pos[:,2],s=0.2,color='b',alpha=0.1)
        
    
    fig.savefig(f'{output}/{cluslast}/{snap}_Proj{clusID}.png')


# Function to make fits file containing the projected positions of the clusters
# called from savefits function
def make_fits(comp,proj,ds,width,res,kind,snap):


    # read the virial radius and make a sphere
    r200_snap = r200.loc[snap]['R200']
    sphere = ds.sphere([0, 0, 0], (r200_snap, "Mpc"))

    # make the fits file for a given projection and component
    prjpd_fits = yt.FITSParticleProjection(
            ds, proj, (comp, "particle_mass"),density=True, deposition="ngp",length_unit='Mpc',image_res=[res,res],data_source=sphere)

    # name the fits file
    prjpd_fits.change_image_name("particle_mass", f"{kind}_{comp}_density_{proj}")

    # update the header
    prjpd_fits.update_header(f"{kind}_{comp}_density_{proj}", "scale", f"{width}/{res} Mpc/pixel")

    # set the units
    prjpd_fits.set_unit(f"{kind}_{comp}_density_{proj}", "Msun/kpc**2")

    return prjpd_fits

# Function to save fits files
def savefits(snap,cluslast,clusID):

    # This function saves the FITS file from the hdf5 data.
    # It uses the yt for loading a generic particle dataset of stars and dm
    # but we right now use only stars for saving in the FITS file.
    # then saves it as fits. The units for projected mass density
    # is Msun/kpc**2

    clusatt = dict(f[f'/{clusID}/'].attrs.items())

    # Initialising global arrays for position and mass
    # of DM and stars
    
    # Stars
    posstarx = np.array([])
    posstary = np.array([])
    posstarz = np.array([])
    mstar= np.array([])
    posstarx_rest = np.array([])
    posstary_rest = np.array([])
    posstarz_rest = np.array([])
    mstar_rest= np.array([])
    posstarx_icl = np.array([])
    posstary_icl = np.array([])
    posstarz_icl = np.array([])
    mstar_icl = np.array([])

    # DM
    posdmy = np.array([])
    posdmz = np.array([])
    posdmx = np.array([])
    mdm= np.array([])
    posdmy_rest = np.array([])
    posdmz_rest = np.array([])
    posdmx_rest = np.array([])
    mdm_rest= np.array([])
    posdmx_icl = np.array([])
    posdmy_icl = np.array([])
    posdmz_icl = np.array([])
    mdm_icl = np.array([])

    # gas
    posgasy = np.array([])
    posgasz = np.array([])
    posgasx = np.array([])
    mgas= np.array([])
    posgasy_rest = np.array([])
    posgasz_rest = np.array([])
    posgasx_rest = np.array([])
    mgas_rest= np.array([])
    posgasx_icl = np.array([])
    posgasy_icl = np.array([])
    posgasz_icl = np.array([])
    mgas_icl = np.array([])



    # print('galaxy IDs?: ',f[f'/{clusID}/'].keys())
    # get total mass of first galaxy
    mid=next(iter(f[f'/{clusID}/'].keys()))
    # print('mid at start: ',mid)

    mostmass=f[f'/{clusID}/{mid}/'].attrs['mtot']

    # print('mass at start: ', mostmass)
    
    for gal in f[f'/{clusID}/'].keys():
        if(gal=='ICL'):
            break 
        print('mass: ',f[f'/{clusID}/{gal}'].attrs['mtot'])
        if f[f'/{clusID}/{gal}'].attrs['mtot']>mostmass:
                mostmass = f[f'/{clusID}/{gal}'].attrs['mtot']
                mid = gal 

    # print('mid at end: ',mid)
    # mid contains the ID of the BCG

    for gal in f[f'/{clusID}/'].keys():

        
        # Stars
        poss =  f[f'/{clusID}/{gal}/posstar']
        mstari =  f[f'/{clusID}/{gal}/massstar']
        posstarx = np.append(posstarx,poss[:,0]) 
        posstary= np.append(posstary,poss[:,1])
        posstarz = np.append(posstarz,poss[:,2])
        mstar = np.append(mstar,mstari)

        # DM
        posdm =  f[f'/{clusID}/{gal}/posdm']
        mdmi =  f[f'/{clusID}/{gal}/massdm']
        posdmx = np.append(posdmx,posdm[:,0])
        posdmy= np.append(posdmy,posdm[:,1])
        posdmz = np.append(posdmz,posdm[:,2])
        mdm = np.append(mdm,mdmi)

        # Gas
        posg =  f[f'/{clusID}/{gal}/posgas']
        mgasi =  f[f'/{clusID}/{gal}/massgas']
        posgasx = np.append(posgasx,posg[:,0])
        posgasy= np.append(posgasy,posg[:,1])
        posgasz = np.append(posgasz,posg[:,2])
        mgas = np.append(mgas,mgasi)


        

        # Check if this is BCG
        if gal==mid:
            data_bcg = {
                    ("star","particle_position_x"): poss[:,0]- clusatt['pos'][0],
                    ("star","particle_position_y"): poss[:,1]- clusatt['pos'][1],
                    ("star","particle_position_z"): poss[:,2]- clusatt['pos'][2],
                    ("star","particle_mass"): mstari,
                    ("dm","particle_position_x"): posdm[:,0]- clusatt['pos'][0],
                    ("dm","particle_position_y"): posdm[:,1]- clusatt['pos'][1],
                    ("dm","particle_position_z"): posdm[:,2]- clusatt['pos'][2],
                    ("dm","particle_mass"): mdmi,
                    ("gas","particle_position_x"): posg[:,0]- clusatt['pos'][0],
                    ("gas","particle_position_y"): posg[:,1]- clusatt['pos'][1],
                    ("gas","particle_position_z"): posg[:,2]- clusatt['pos'][2],
                    ("gas","particle_mass"): mgasi
                }
        
        # Rest of galaxies
        elif (gal!=mid) and (gal!='ICL'):

            
            # Stars
            poss =  f[f'/{clusID}/{gal}/posstar']
            mstari =  f[f'/{clusID}/{gal}/massstar']
            posstarx_rest = np.append(posstarx_rest,poss[:,0]) 
            posstary_rest= np.append(posstary_rest,poss[:,1])
            posstarz_rest = np.append(posstarz_rest,poss[:,2])
            mstar_rest = np.append(mstar_rest,mstari)

            # DM
            posdm =  f[f'/{clusID}/{gal}/posdm']
            mdmi =  f[f'/{clusID}/{gal}/massdm']
            posdmx_rest = np.append(posdmx_rest,posdm[:,0])
            posdmy_rest= np.append(posdmy_rest,posdm[:,1])
            posdmz_rest = np.append(posdmz_rest,posdm[:,2])
            mdm_rest = np.append(mdm_rest,mdmi)

            # Gas
            posg =  f[f'/{clusID}/{gal}/posgas']
            mgasi =  f[f'/{clusID}/{gal}/massgas']
            posgasx_rest = np.append(posgasx_rest,posg[:,0])
            posgasy_rest= np.append(posgasy_rest,posg[:,1])
            posgasz_rest = np.append(posgasz_rest,posg[:,2])
            mgas_rest = np.append(mgas_rest,mgasi)


        if gal=='ICL':

            
            # This is ICL case
            # Stars
            poss =  f[f'/{clusID}/{gal}/posstar']
            #print(poss.shape)
            mstari =  f[f'/{clusID}/{gal}/massstar']
            posstarx_icl = np.append(posstarx_icl,poss[:,0]) 
            posstary_icl= np.append(posstary_icl,poss[:,1])
            posstarz_icl = np.append(posstarz_icl,poss[:,2])
            mstar_icl = np.append(mstar_icl,mstari)

            # DM
            posdm =  f[f'/{clusID}/{gal}/posdm']
            mdmi =  f[f'/{clusID}/{gal}/massdm']
            posdmx_icl = np.append(posdmx_icl,posdm[:,0])
            posdmy_icl = np.append(posdmy_icl,posdm[:,1])
            posdmz_icl = np.append(posdmz_icl,posdm[:,2])
            mdm_icl = np.append(mdm_icl,mdmi)

            # Gas
            posg =  f[f'/{clusID}/{gal}/posgas']
            mgasi =  f[f'/{clusID}/{gal}/massgas']
            posgasx_icl = np.append(posgasx_icl,posg[:,0])
            posgasy_icl = np.append(posgasy_icl,posg[:,1])
            posgasz_icl = np.append(posgasz_icl,posg[:,2])
            mgas_icl = np.append(mgas_icl,mgasi)
        
        
        

        

    data_rest = {
        ("star","particle_position_x"): posstarx_rest- clusatt['pos'][0],
        ("star","particle_position_y"): posstary_rest- clusatt['pos'][1],
        ("star","particle_position_z"): posstarz_rest- clusatt['pos'][2],
        ("star","particle_mass"): mstar_rest,
        ("dm","particle_position_x"): posdmx_rest- clusatt['pos'][0],
        ("dm","particle_position_y"): posdmy_rest- clusatt['pos'][1],
        ("dm","particle_position_z"): posdmz_rest- clusatt['pos'][2],
        ("dm","particle_mass"): mdm_rest,
        ("gas","particle_position_x"): posgasx_rest- clusatt['pos'][0],
        ("gas","particle_position_y"): posgasy_rest- clusatt['pos'][1],
        ("gas","particle_position_z"): posgasz_rest- clusatt['pos'][2],
        ("gas","particle_mass"): mgas_rest
    }

    data_icl = {
        ("star","particle_position_x"): posstarx_icl- clusatt['pos'][0],
        ("star","particle_position_y"): posstary_icl- clusatt['pos'][1],
        ("star","particle_position_z"): posstarz_icl- clusatt['pos'][2],
        ("star","particle_mass"): mstar_icl,
        ("dm","particle_position_x"): posdmx_icl- clusatt['pos'][0],
        ("dm","particle_position_y"): posdmy_icl- clusatt['pos'][1],
        ("dm","particle_position_z"): posdmz_icl- clusatt['pos'][2],
        ("dm","particle_mass"): mdm_icl,
        ("gas","particle_position_x"): posgasx_icl- clusatt['pos'][0],
        ("gas","particle_position_y"): posgasy_icl- clusatt['pos'][1],
        ("gas","particle_position_z"): posgasz_icl- clusatt['pos'][2],
        ("gas","particle_mass"): mgas_icl
    }

    data_all = {
        ("star","particle_position_x"): posstarx- clusatt['pos'][0],
        ("star","particle_position_y"): posstary- clusatt['pos'][1],
        ("star","particle_position_z"): posstarz- clusatt['pos'][2],
        ("star","particle_mass"): mstar,
        ("dm","particle_position_x"): posdmx- clusatt['pos'][0],
        ("dm","particle_position_y"): posdmy- clusatt['pos'][1],
        ("dm","particle_position_z"): posdmz- clusatt['pos'][2],
        ("dm","particle_mass"): mdm,
        ("gas","particle_position_x"): posgasx- clusatt['pos'][0],
        ("gas","particle_position_y"): posgasy- clusatt['pos'][1],
        ("gas","particle_position_z"): posgasz- clusatt['pos'][2],
        ("gas","particle_mass"): mgas
    }

    # get width that can encompass all particles
    res=1024
    width=0.2
    result = None
    while result is None:
        try:
            # connect
            bbox = np.array([[-width,width], [-width, width], [-width, width]])

            ds_all = yt.load_particles(data_all, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
            result = yt.ParticleProjectionPlot(ds_all,'x',("star","particle_mass"))
            
        except:
            width=width+0.2

            
            pass
    
    ds_rest = yt.load_particles(data_rest, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
    ds_bcg = yt.load_particles(data_bcg, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
    ds_icl = yt.load_particles(data_icl, length_unit='Mpc', mass_unit='Msun', bbox=bbox)


    all_img=[]
    for cm in ['star','dm','gas']:

        for proj in ['x','y','z']:

            # print(cm,proj,"M here all")
            # FITS projection for all particles
            prjs = make_fits(cm,proj,ds_all,width,res,'all',snap)
            all_img.append(prjs)

            
            # print(cm,proj,"M here bcg")
            prjs = make_fits(cm,proj,ds_bcg,width,res,'bcg',snap)
            all_img.append(prjs)

            # print(cm,proj,"M here rest")
            prjs = make_fits(cm,proj,ds_rest,width,res,'rest',snap)
            all_img.append(prjs)

            # print(cm,proj,"M here icl")
            prjs = make_fits(cm,proj,ds_icl,width,res,'ICL',snap)
            all_img.append(prjs)
            



    # combined fits for all particles
    prj_fits_comb = yt.FITSImageData.from_images(all_img)
    prj_fits_comb.writeto(f'{output}/{cluslast}/{snap}_{clusID}_r200.fits', overwrite=True)




# Main code 

for snap in range(1,31):
    # save fits
    with h5py.File(f"{output}clusters{snap}.hdf5", "r") as f:
        savefits(snap,0,0)        

                     
