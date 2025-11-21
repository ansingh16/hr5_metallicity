import numpy as np
import pandas as pd 
import h5py
import yt 
# from yt.units import Msun, Mpc
import tqdm 
import os
from astropy.io import fits
import tqdm
from projection_utils import make_ds_for_view,make_fits, fibonacci_sphere
#plt.style.use('/home/ankitsingh/hr5_agn/paper_style.mplstyle')

# low information printing for yt
# avoiding unnecessary information can save execution time 
yt.set_log_level(40)


# Function to save fits files
def savefits(cluslast,clusID):

    

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



        
        galids = list(f[f'/{clusID}/'].keys())
            
            
        galids.remove('ICL')

        mid=-1


        mostmass=0
        totmass=0
        #print("starting",mostmass)
        for gal in galids:
            totmass = f[f'/{clusID}/{gal}'].attrs['mstar'] + totmass
            if f[f'/{clusID}/{gal}'].attrs['mstar']>mostmass:
                mostmass = f[f'/{clusID}/{gal}'].attrs['mstar']
                mid = gal 
            #print(mid,mostmass)
        totmass = totmass + f[f'/{clusID}/ICL'].attrs['mstar']
        print(f"Total mass of cluster: {np.log10(totmass)} and cluster mstar: {np.log10(f[f'/{clusID}/'].attrs['mstar'])}")
        print('mid at end: ',mid)

        for gal in f[f'/{clusID}/'].keys():

            # print(gal)
            
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

        print(data_all[("star","particle_mass")].sum())
        # get width that can encompass all particles
        res=1024
        
    
        # Example integration into your block
        view_dirs = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]

        if not TEST:
            num_dirs = 10
        else:
            num_dirs = 10

        # include the three Cartesian axes, then fill the rest
        num_remaining = num_dirs - len(view_dirs)
        if num_remaining > 0:
            fib_dirs = fibonacci_sphere(num_remaining)
            for v in fib_dirs:
                view_dirs.append(tuple(v.tolist()))

        # produce projections for each viewing direction
        view_index = 0
        for view_dir in tqdm.tqdm(view_dirs):
            tag = f"view{view_index:03d}"

                        
            # Labels for subsets and particle components
            subsets = ["all", "bcg", "rest", "icl"]
            components = ["star", "dm", "gas"]

            
            # Create yt datasets for each subset
            print("processing all")
            ds_all_view,width, bbox  = make_ds_for_view(data_dict=data_all, view_dir=view_dir,bbox=None)
            print("processing bcg")
            ds_bcg_view  = make_ds_for_view(data_bcg,  view_dir,bbox)
            print("processing rest")
            ds_rest_view = make_ds_for_view(data_rest, view_dir,bbox)
            print("processing icl")
            ds_icl_view  = make_ds_for_view(data_icl,  view_dir,bbox)

            

            
            # Store all FRB arrays
            arrays = []

            for subset, ds_view in zip(subsets, [ds_all_view, ds_bcg_view, ds_rest_view, ds_icl_view]):
                for comp in components:
                    arr = make_fits(comp, "x", ds_view, width, res, f"{subset}_{tag}")
                    arrays.append((subset, comp, arr))

            # Create a multi-extension FITS file
            hdulist = fits.HDUList([fits.PrimaryHDU()])  # Primary HDU is empty

            for subset, comp, arr in arrays:
                hdu_name = f"{subset}_{comp}"
                hdu = fits.ImageHDU(data=arr, name=hdu_name)
                hdulist.append(hdu)

            for hdu in hdulist[1:]:  # skip the primary
                hdu.header['SNAP'] = snap
                hdu.header['CLUSID'] = clusID
                hdu.header['TAG'] = tag
                hdu.header['WIDTH'] = 2*width
                hdu.header['RES'] = res
                hdu.header['VIEW'] = str(view_dir)
            # hdulist.writeto("combined.fits", overwrite=True)

            # combine these images into a single FITS file for this view and write
            try:
                # prj_fits_comb = yt.FITSImageData.from_images(imgs)
                if TEST:
                    out_dir = f"{output}/test/{cluslast}"
                    
                else:
                    out_dir = f"{output}/{cluslast}"
                
                # check if output directory exists
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir, exist_ok=True)
                
                hdulist.writeto(f'{out_dir}/{snap}_{clusID}_{tag}.fits', overwrite=True)
            except Exception as e:
                print(f"Error writing FITS for snap {snap} cluster {clusID} view {view_index} with error {e}")

            # increment view counter and continue — do not save directions anywhere
            view_index += 1
                



    


cos='CDM'
halo='12'
TEST=True
print(f'Processing halo {halo} for cosmo {cos}')
output = f'/scratch/ankitsingh/Galaxy_catalogs/ICL_data/CEAGLE/Data/Output_new/{cos}/halo{halo}/'

# check if output directory exists
if not os.path.exists(output):
    os.makedirs(output)

# Main code 
if TEST:
    snapshots = [1]
else:
    snapshots = list(range(32))

for snap in snapshots:
    clusters_file = f"{output}clusters{snap}.hdf5"
    if not os.path.exists(clusters_file):
        print(f"Clusters file not found for snap {snap}: {clusters_file}")
        continue

    with h5py.File(clusters_file, "r") as f:
            # print keys
            clusters= list(f.keys())
            # remove status key
            clusters.remove('status')
            print(f"working on snap {snap} for halo {halo} for cosmo {cos}")
            for clusID in clusters:
                cluslast = snap
                savefits(cluslast, clusID)
                
        