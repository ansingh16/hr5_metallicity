import numpy as np
import pandas as pd 
import h5py
import yt 
from yt.units import Msun, Mpc
import tqdm 
import matplotlib.pyplot as plt 
plt.style.use('/home/ankitsingh/hr5_agn/paper_style.mplstyle')

# low information printing for yt
# avoiding unnecessary information can save execution time 
yt.set_log_level(40)

import multiprocessing

import configparser

# get all snapshots
import os 
import glob

import configparser



# load up the parameter file
parser = configparser.ConfigParser()
parser.read('/home/ankitsingh/hr5_metalicity/params.ini')

output = parser.get('Paths','outdir')
snapfiles = parser.get('Paths','snapfiles')

# read all the snapshot dat files
files =  [filename for filename in os.listdir(snapfiles) if os.path.isfile(os.path.join(snapfiles, filename))]


snaps = [os.path.splitext(file)[0] for file in files]

# snaps.remove('296')

snaps.sort(reverse=True)
# Function for saving png corresponding to each cluster
# at each snapshot

# Function to make fits file containing the projected positions of the clusters
# called from savefits function
def make_fits(comp,proj,ds,width,res,kind):

    prjpd_fits = yt.FITSParticleProjection(
            ds, proj, (comp, "particle_mass"), deposition="ngp",length_unit='Mpc',image_res=[res,res])

    prjpd_fits.change_image_name("particle_mass", f"{kind}_{comp}_mass_{proj}")

    prjpd_fits.update_header(f"{kind}_{comp}_mass_{proj}", "scale", f"{width}/{res} Mpc/pixel")

    prjpd_fits.set_unit(f"{kind}_{comp}_mass_{proj}", "Msun")

    return prjpd_fits

# Function to save fits files
def savefits( snap,clusID,cluslast):

    # This function saves the FITS file from the hdf5 data.
    # It uses the yt for loading a generic particle dataset of stars and dm
    # but we right now use only stars for saving in the FITS file.
    # then saves it as fits. The units for projected mass density
    # is Msun/kpc**2


    hid = pd.read_csv(f'{snapfiles}{snap}.dat')
    
    if clusID in hid.to_numpy():
                                
            # read snap hdf5 file containing the data
            
            f = h5py.File(f"{output}clusters{snap}.hdf5", "r")
            print(f"For cluster = {cluslast}, Reading {snap} and Halo = {clusID}")
                                

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
            posstarx_sudo = np.array([])
            posstary_sudo = np.array([])
            posstarz_sudo = np.array([])
            mstar_sudo= np.array([])

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
            posdmx_sudo = np.array([])
            posdmy_sudo = np.array([])
            posdmz_sudo = np.array([])
            mdm_sudo = np.array([])

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
            posgasx_sudo = np.array([])
            posgasy_sudo = np.array([])
            posgasz_sudo = np.array([])
            mgas_sudo = np.array([])


            # get total mass of first galaxy
            mid=-1


            mostmass=0

            for gal in f[f'/{clusID}/'].keys():
                if f[f'/{clusID}/{gal}'].attrs['mstar']>mostmass:
                        mostmass = f[f'/{clusID}/{gal}'].attrs['mstar']
                        mid = gal 

            # mass of most massive galaxy
            mass_most_mass = f[f'/{clusID}/{mid}'].attrs['mstar']

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

                    mass_gal = f[f'/{clusID}/{gal}'].attrs['mstar']

                    # check if ratio of stellar mass to most massive gal mass

                    if (mass_gal/mass_most_mass)>0.2:
                        
                        # Stars
                        posstarx_sudo = np.append(posstarx_sudo,poss[:,0]) 
                        posstary_sudo= np.append(posstary_sudo,poss[:,1])
                        posstarz_sudo = np.append(posstarz_sudo,poss[:,2])
                        mstar_sudo = np.append(mstar_sudo,mstari)

                        # DM
                        posdmx_sudo = np.append(posdmx_sudo,posdm[:,0])
                        posdmy_sudo= np.append(posdmy_sudo,posdm[:,1])
                        posdmz_sudo = np.append(posdmz_sudo,posdm[:,2])
                        mdm_sudo = np.append(mdm_sudo,mdmi)

                        # gas

                        posgasx_sudo = np.append(posgasx_sudo,posg[:,0])
                        posgasy_sudo= np.append(posgasy_sudo,posg[:,1])
                        posgasz_sudo = np.append(posgasz_sudo,posg[:,2])
                        mgas_sudo = np.append(mgas_sudo,mgasi)


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


            data_sudo = {
                ("star","particle_position_x"): posstarx_sudo- clusatt['pos'][0],
                ("star","particle_position_y"): posstary_sudo- clusatt['pos'][1],
                ("star","particle_position_z"): posstarz_sudo- clusatt['pos'][2],
                ("star","particle_mass"): mstar_sudo,
                ("dm","particle_position_x"): posdmx_sudo- clusatt['pos'][0],
                ("dm","particle_position_y"): posdmy_sudo- clusatt['pos'][1],
                ("dm","particle_position_z"): posdmz_sudo- clusatt['pos'][2],
                ("dm","particle_mass"): mdm_sudo,
                ("gas","particle_position_x"): posgasx_sudo- clusatt['pos'][0],
                ("gas","particle_position_y"): posgasy_sudo- clusatt['pos'][1],
                ("gas","particle_position_z"): posgasz_sudo- clusatt['pos'][2],
                ("gas","particle_mass"): mgas_sudo
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
            ds_sudo = yt.load_particles(data_sudo, length_unit='Mpc', mass_unit='Msun', bbox=bbox)

            
            all_img=[]
            for cm in ['star','dm','gas']:

                for proj in ['x','y','z']:

                    # print(cm,proj,"M here all")
                    # FITS projection for all particles
                    # prjs = make_fits(cm,proj,ds_all,width,res,'all')
                    # all_img.append(prjs)

                    
                    # # print(cm,proj,"M here bcg")
                    # prjs = make_fits(cm,proj,ds_bcg,width,res,'bcg')
                    # all_img.append(prjs)

                    # # print(cm,proj,"M here rest")
                    # prjs = make_fits(cm,proj,ds_rest,width,res,'rest')
                    # all_img.append(prjs)

                    # # print(cm,proj,"M here icl")
                    # prjs = make_fits(cm,proj,ds_icl,width,res,'ICL')
                    # all_img.append(prjs)
                    
                    prjs = make_fits(cm,proj,ds_sudo,width,res,'sudo')
                    all_img.append(prjs)
                    



            # # combined fits for all particles
            # prj_fits_comb = yt.FITSImageData.from_images(all_img)
            # prj_fits_comb.writeto(f'/home/jwoo/{snap}_{clusID}.fits', overwrite=True)



#Main code
clusters = pd.read_csv(f'{snapfiles}296.dat')
# loop over clusters of interest
for clus in tqdm.tqdm([1664541]):#clusters['HostHaloID']):


            
            # read mass accretion files
            MAHfile = pd.read_csv(f"{output}/{clus}_MAH.csv")
            # # get snaps that come after half mass accumulation
            MAHfile = MAHfile[MAHfile['snap']>=int(min(snaps))]

            # # last snapshot is already done in previous analysis
            # beyond_half = beyond_half[beyond_half['snap'] != 296]
            
            if not os.path.exists(f"{output}/{clus}"):
                os.mkdir(f"{output}/{clus}")
            
            # Loop over the progenitors and make fits 
            args= [(snap,haloid,clus) for snap,haloid in zip(MAHfile['snap'].values,MAHfile['HostHaloID'].values)]
                        
            print(args)            
            
            # save fits
            with multiprocessing.Pool(1) as pool:
                pool.starmap(savefits,args)        

# savefits(296,1664541,1664541)
