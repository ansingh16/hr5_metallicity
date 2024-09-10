import h5py
import numpy as np
import configparser
import re 
import yt 
import astropy.units as u
import pandas as pd 
from scipy import stats
from scipy.optimize import curve_fit
import tqdm 
import matplotlib.pyplot as plt
plt.style.use('../paper_style.mplstyle')
import logging
import math


yt.mylog.setLevel(logging.ERROR)


# load up the parameter file
parser = configparser.ConfigParser()
parser.read('../params.ini')


outdir = parser.get('Paths','outdir')
h0 = float(parser.get('Setting','h0'))

# Function to pad lists to the same length
def pad_list(series):
    max_len = max(len(x) for x in series)
    return np.array([np.pad(x, (0, max_len - len(x)), 'constant', constant_values=np.nan) for x in series])

# Compute the element-wise median with padding
def median_list_column(series):
    padded_array = pad_list(series)
    return np.nanmedian(padded_array, axis=0)



class Cluster:
    """
    A class for the cluster data
    
    Attributes
    ----------

    clusID: int
        cluster id at the snapshot
    clus_mdm: float 
        total dark matter mass of the cluster
    clus_mgas: float
        total gas mass of the cluster
    clus_msink: float 
        total sink particle mass of cluster 
    clus_mstar: float 
        total stellar mass of cluster
    clus_mtot: float 
        total mass of the cluster
    clus_ngas: int
        total number of gas particles in cluster 
    clus_nsink: int
        total number of sink particles in cluster
    clus_nstar: int 
        total number of star particles in cluster
    clus_nsub: int 
        total number of subhalos (galaxies) cluster
    clus_pos: array of float
        position of cluster in cMpc in the simulation
    clus_vel: array of float
        velocity of the cluster in km/s 
    f: hdf5 file type
        excess to the raw hdf5 file
    snap: int 
        snapshot number 

    Methods
    -------

    get_all_parts: 
        function to get all particles of a type in cluster
    get_alldat_gal: 
        get all data for a galaxy
    get_galids: 
        get all galaxy ids in a cluster
    save_yt_dataset: 
        get yt dataset for a cluster

    """

    def __init__(self,snapno,clusno):
        """
        Construct all necessary attributes for a cluster.

        Parameters
        ----------
        snapno : int
            snapshot number
        clusno : int
            cluster number
        """

        self.snap=snapno
        self.clusID=clusno

        self.f = h5py.File(f"{outdir}clusters{self.snap}.hdf5", "r")

        self.bcgid = self._BCG_ID()

        attrs = self.f[f'/{self.clusID}/'].attrs

        for att in attrs.keys():
            setattr(self, f"clus_{att}",attrs[att]) 
        
    def _BCG_ID(self):
        """
        Get the ID of the central galaxy and return it.

        :return: The ID of the central galaxy.
        """
        # get the ID of central galaxy
        # get total mass of first galaxy
        gallis = self.get_galids()
        
        mid=-1


        mostmass=0
        #print("starting",mostmass)
        for gal in gallis:
            if self.f[f'/{self.clusID}/{gal}'].attrs['mstar']>mostmass:
                    mostmass = self.f[f'/{self.clusID}/{gal}'].attrs['mstar']
                    mid = gal 
            #print(mid,mostmass)

        return mid 
    
    def get_galids(self): 
        """
        Returns a list of galaxy IDs from the `self.f` dictionary that corresponds to the cluster with ID `self.clusID`.

        :return: A list of galaxy IDs.
        """

        galids = list(self.f[f'/{self.clusID}/'].keys())
        
        
        galids.remove('ICL')

        # get low resolution galaxies
        low_res = []
        for gal in galids:
            if self.f[f'/{self.clusID}/{gal}/'].attrs['mstar']<2.0e9:
               low_res.append(gal)
        # remove low resolution galaxies
        galids = [gal for gal in galids if gal not in low_res]

        return galids
    
    

    def get_alldat_gal(self,galist):
        """
        This function takes in a list of galaxy IDs or a single galaxy ID and returns a Galaxy object or a list of Galaxy objects respectively. 
        
        Args:
        - galist (list or int): A list of galaxy IDs or a single galaxy ID.
        
        Returns:
        - outgal (list): A list of Galaxy objects if galist is a list.
        - gal (Galaxy object): A Galaxy object if galist is an int.
        """

        if isinstance(galist, list):
            
            outgal=[]
            for galid in galist:
                # print(f"\rProcessing galaxy {galid}", end='')
                gal = Galaxy(self.snap,self.clusID,galid)
                
                for part in ['gas','star','dm']:
                            
                    vars=list(self.f[f'/{self.clusID}/{galid}/'].keys())
                    partvar = [var for var in vars if re.search(part, var)]
                    for var in partvar:
                        stm=var.replace(part,'')
                        setattr(gal,f'{part}_{stm}',self.f[f'/{self.clusID}/{galid}/{var}'])
                            
                outgal.append(gal)
            return outgal
        else:
            # here galist is an int
            # print(f"\rProcessing galaxy {galist}", end='')
            gal = Galaxy(self.snap,self.clusID,galist)
            
            for part in ['gas','star','dm']:
                            
                vars=list(self.f[f'/{self.clusID}/{galist}/'].keys())
                partvar = [var for var in vars if re.search(part, var)]
                for var in partvar:
                    stm=var.replace(part,'')
                    setattr(gal,f'{part}_{stm}',self.f[f'/{self.clusID}/{galist}/{var}'])

            return gal

    def get_all_parts(self,partype):
        """
        This function returns an instance of the Galaxy class containing all parts of the specified type.
        
        Args:
            partype (str): The type of part to retrieve.
        
        Returns:
            Galaxy: An instance of the Galaxy class containing all parts of the specified type.
        """

        gal = Galaxy(self.snap, self.clusID)
        galids = self.get_galids()

        for galid in galids:
            galtmp = self.get_alldat_gal(galid)
            vars = [var for var in self.f[f'/{self.clusID}/{galid}/'].keys() if re.search(partype, var)]
            
            for var in vars:
                stm = var.replace(partype, '')
                dat_comb = np.concatenate((getattr(gal, f'{partype}_{stm}'), getattr(galtmp, f'{partype}_{stm}')), axis=0)
                setattr(gal, f'{partype}_{stm}', dat_comb)

        return gal
                
            

    def save_yt_dataset(self):
        """
        Saves the yt dataset of the BCG, ICM, and rest of the galaxies.
        
        
        Returns:
        --------
        ds_all: yt dataset
            Dataset containing all the particles.
        ds_rest: yt dataset
            Dataset containing the rest of the galaxies particles.
        ds_bcg: yt dataset
            Dataset containing the BCG particles.
        ds_icm: yt dataset
            Dataset containing the ICM particles.
        """

        # get galaxies as BCG, ICM and rest of the galaxies
        #BCG
        mid = self._BCG_ID()
        bcg = self.get_alldat_gal(mid)
        # Rest
        galids= list(self.f[f'/{self.clusID}/'].keys())
        
        allgal = self.get_alldat_gal(galids)
        restgal = [galid for galid in galids if (galid != mid) & (galid !='ICL')]  
        restgal = self.get_alldat_gal(restgal)
        #ICM
        icm = self.get_alldat_gal('ICL')
        
        
        #get variable list
        varbls=list(vars(icm))
        varbls = [x for x in varbls if re.search('_',x)]
        starvar = [var.replace('star_','') for var in varbls if re.search('star_',var)and (var[0]!=r'_')]
        gasvar = [var.replace('gas_','') for var in varbls if re.search('gas_',var)and (var[0]!=r'_')]
        dmvar = [var.replace('dm_','') for var in varbls if re.search('dm_',var)and (var[0]!=r'_')]
        
        #These variables are to be added to the dataset apart from pos,vel,mass
        var_s =set(starvar)-set(['pos_com','mass','vel'])
        var_g =set(gasvar)-set(['pos_com','mass','vel'])
        var_d = set(dmvar)-set(['pos_com','mass','vel'])
        
        
        #fill data in icm and bcg
        data_bcg = {}
        data_icm={}
        data_rest={}
        data_all={}
        data_dict={'BCG':data_bcg,'ICM':data_icm,'Rest':data_rest,'All':data_all}

        component_dict = {'BCG':bcg,'ICM':icm,'Rest':restgal,'All':allgal}

        for component in component_dict.keys():
            print(f"filling {component} data")

            glx = component_dict[component]
            
            if not isinstance(glx,list):
                # it is bcg or icm
                
                for part in ['gas','star','dm']:
                        data_dict[component][(f"{part}","particle_mass")] = getattr(glx,f'{part}_mass')[:]
                        for i,dir in enumerate(['x','y','z']):
                            data_dict[component][(f"{part}",f"particle_position_{dir}")] = getattr(glx,f'{part}_pos_com')[:,i]
                        
                        if part=='star':
                            for var in var_s:
                                data_dict[component][(f"{part}",f"{var}")] = getattr(glx,f'{part}_{var}')[:]
                        if part=='gas':
                            for var in var_g:
                                data_dict[component][(f"{part}",f"{var}")] = getattr(glx,f'{part}_{var}')[:]
                        if part=='dm':
                            for var in var_d:
                                data_dict[component][(f"{part}",f"{var}")] = getattr(glx,f'{part}_{var}')[:]
            else:
                # it is rest of galaxies and all galaxies list
                for part in ['gas','star','dm']:
                    
                    
                    
                    data_dict[component][(f"{part}","particle_mass")] = np.concatenate([getattr(gal,f'{part}_mass') for gal in glx],axis=0)
                    
                    

                    for i,dir in enumerate(['x','y','z']):
                            if component == 'All':
                                data_dict[component][(f"{part}",f"particle_position_{dir}")] = np.concatenate([getattr(gal,f'{part}_pos_com') for gal in glx],axis=0)[:,i]

                            elif component == 'Rest':
                                data_dict[component][(f"{part}",f"particle_position_{dir}")] = np.concatenate([getattr(gal,f'{part}_pos') - self.clus_pos for gal in glx],axis=0)[:,i]
                            
                    if part=='star':
                            for var in var_s:
                                data_dict[component][(f"{part}",f"{var}")] = np.concatenate([getattr(gal,f'{part}_{var}') for gal in glx],axis=0)[:]
                    if part=='gas':
                            for var in var_g:
                                data_dict[component][(f"{part}",f"{var}")]= np.concatenate([getattr(gal,f'{part}_{var}') for gal in glx],axis=0)[:]
                    if part=='dm':
                            for var in var_d:
                                data_dict[component][(f"{part}",f"{var}")]= np.concatenate([getattr(gal,f'{part}_{var}') for gal in glx],axis=0)[:]
            
        
        jjj=0
        for dat in data_dict.keys():

                print(f"{dat} has {data_dict[dat]['star','particle_mass'].shape} particles")
                jjj = jjj+data_dict[dat]['star','particle_mass'].shape[0]


        # get width that can encompass all particles
        res=1024
        width=2
        result = None

        while result is None:
            try:
                    
                    bbox = np.array([[-width,width], [-width, width], [-width, width]])
                    ds_all = yt.load_particles(data_all, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
                    result = yt.ParticleProjectionPlot(ds_all,'x',("star","particle_mass"))
            except:        
                    width=width+0.2
        
        # print(result, width)
        
        bbox = np.array([[-width,width], [-width, width], [-width, width]])

        ds_all = yt.load_particles(data_dict['All'], length_unit='Mpc', mass_unit='Msun', bbox=bbox)
        result = yt.ParticleProjectionPlot(ds_all,'x',("star","particle_mass"))
                
        

        ds_rest = yt.load_particles(data_dict['Rest'], length_unit='Mpc', mass_unit='Msun', bbox=bbox)
        ds_bcg = yt.load_particles(data_dict['BCG'], length_unit='Mpc', mass_unit='Msun', bbox=bbox)
        ds_icm = yt.load_particles(data_dict['ICM'], length_unit='Mpc', mass_unit='Msun', bbox=bbox)

        return ds_all,ds_rest,ds_bcg,ds_icm




# Class for handling galaxy data
class Galaxy(Cluster):
    """
    A class for galaxy data

    dm_mass: float arr(gal_ndm)
        array of dark matter mass 
    dm_pos: float arr(gal_ndm,3)
        array of dark matter position 
    dm_pos_com: float arr(gal_ndm,3)
        array of dark matter position in center of mass 
    dm_vel: float arr(gal_ndm,3)
        array of dark matter velocity 
    galID: int
        galaxy ID
    gal_mdm: float
        total dark matter mass of galaxy
    gal_mgas: float
        total gas mass of galaxy
    gal_msink: float
        total sink mass of galaxy
    gal_mstar: float
        total stellar mass of galaxy
    gal_mtot: float
        total mass of galaxy
    gal_ngas: int
        total number of gas particles in galaxy
    gal_nsink: int
        total number of sink particles in galaxy
    gal_nstar: int
        total number of stellar particles in galaxy
    gal_pos: float arr(3)
        galaxy position 
    gal_vel: float arr(3)
        array of gas particle velocity 
    gas_fe: float arr(gal_ngas)
        array of fe 
    gas_h: float arr(gal_ngas)
        array of h 
    gas_mass: float arr(gal_gas)
        array of gas particle mass 
    gas_o: float arr(gal_ngas)
        array of o 
    gas_pos: float arr(gal_ngas,3)
        array of gas particle position 
    gas_pos_com: float arr(gal_ngas,3)
        array of gas particle position in center of mass 
    gas_t: float arr(gal_ngas)
        array of gas temperature 
    gas_vel: float arr(gal_ngas,3)
        array of gas particle velocity 
    gas_z: float arr(gal_ngas)
        array of gas metallicity 
    rcom_dm:float:arr(gal_ndm)
        array of dark matter particle distances from origin in COM frame
    rcom_gas: float arr(gal_ngas)
        array of gas particle distances from origin in COM frame
    rcom_star: float arr(gal_nstar)
        array of stellar particle distances from origin in COM frame
    star_mass: float arr(gal_nstar)
        array of stellar particle mass
    star_pos: float arr(gal_nstar,3)
        array of stellar particle position size
    star_pos_com: float arr(gal_nstar,3)
        array of stellar particle position in COM frame 
    star_vel: float arr(gal_nstar,3)
        array of stellar particle velocity 
    star_z: float arr(gal_nstar)
        array of stellar particle metallicty 
    
    """
    def __init__(self,snap,clusID,galid=None):
            super().__init__(snap,clusID)

            """
            Parameters
            ----------
            snap: int
                snapshot
            clusID: int
                cluster ID
            galid: int
                galaxy ID
            """
            
            
            if galid==None:
                self.galID = None
                pass
            else:
                attrs = self.f[f'/{self.clusID}/{galid}'].attrs

                for att in attrs.keys():
                    setattr(self, f"gal_{att}",attrs[att]) 

            self.galID     = galid
            self.star_pos  = np.empty((0, 3))
            self.star_vel  = np.empty((0, 3))
            self.star_mass = np.empty(0)
            self.star_z    = np.empty(0)
            self.dm_pos    = np.empty((0, 3))
            self.dm_vel    = np.empty((0, 3))
            self.dm_mass   = np.empty(0)
            self.gas_pos   = np.empty((0, 3))
            self.gas_mass  = np.empty(0)
            self.gas_vel   = np.empty((0, 3))
            self.gas_fe    = np.empty(0)
            self.gas_h     = np.empty(0)
            self.gas_o     = np.empty(0)
            self.gas_t     = np.empty(0)
            self.gas_z     = np.empty(0)
            
            
    @property
    def star_pos(self):
        return self._star_pos

    @star_pos.setter
    def star_pos(self, value):
        self._star_pos = value
        self._update('star')
        
    @property
    def gas_pos(self):
        return self._gas_pos

    @gas_pos.setter
    def gas_pos(self, value):
        self._gas_pos = value
        self._update('gas')
    
    @property
    def dm_pos(self):
        return self._dm_pos

    @dm_pos.setter
    def dm_pos(self, value):
        self._dm_pos = value
        self._update('dm')

    
        

    def _update(self,part):
        """
        Update the positions of the particles and calculate the center of mass of the selected particle type.
        
        :param part: A string representing the particle type to update. It can be 'star', 'gas', or 'dm'.
        
        :return: None
        """
        
        if self.galID != None:
            if part=='star':
                self.star_pos_com = self.star_pos-self.gal_pos
                self.rcom_star = np.linalg.norm(self.star_pos_com,axis=1)

                
            elif part=='gas':
                self.gas_pos_com = self.gas_pos-self.gal_pos
                self.rcom_gas = np.linalg.norm(self.gas_pos_com,axis=1)
            
            elif part=='dm':
                self.dm_pos_com = self.dm_pos-self.gal_pos
                self.rcom_dm = np.linalg.norm(self.dm_pos_com,axis=1)

    # function to get yt dataset for a single galaxy
    def get_yt_dataset(self):
        """
        Returns a yt dataset for a single galaxy.
        """

        gal = self.get_alldat_gal(self.galID)
        #get variable list
        varbls=list(vars(gal))
        varbls = [x for x in varbls if re.search('_',x)]
        starvar = [var.replace('star_','') for var in varbls if re.search('star_',var)and (var[0]!=r'_')]
        gasvar = [var.replace('gas_','') for var in varbls if re.search('gas_',var)and (var[0]!=r'_')]
        dmvar = [var.replace('dm_','') for var in varbls if re.search('dm_',var)and (var[0]!=r'_')]
        
        #These variables are to be added to the dataset apart from pos,vel,mass
        var_s =set(starvar)-set(['pos_com','mass','vel'])
        var_g =set(gasvar)-set(['pos_com','mass','vel'])
        var_d = set(dmvar)-set(['pos_com','mass','vel'])
        
        # the id of the galaxy in a list
        glx=self.get_alldat_gal(self.galID)
        
        # dictionary to fill the data
        data_dict = {} 
        
        # Loop over the particle types and fill the data dictionary
        # if self.gas_mass>0:
        #     par_types = ['gas','star','dm']
        # else:
        #     par_types = ['star','dm']

        for part in ['gas','star','dm']:
                        data_dict[(f"{part}","particle_mass")] = getattr(glx,f'{part}_mass')[:]
                        for i,dir in enumerate(['x','y','z']):
                            data_dict[(f"{part}",f"particle_position_{dir}")] = getattr(glx,f'{part}_pos_com')[:,i]
                        
                        if part=='star':
                            for var in var_s:
                                data_dict[(f"{part}",f"{var}")] = getattr(glx,f'{part}_{var}')[:]
                        if part=='gas':
                            for var in var_g:
                                data_dict[(f"{part}",f"{var}")] = getattr(glx,f'{part}_{var}')[:]
                        if part=='dm':
                            for var in var_d:
                                data_dict[(f"{part}",f"{var}")] = getattr(glx,f'{part}_{var}')[:]
        

        data_all = data_dict
        
        # if the galaxy is the BCG, increase the width of the plot
        if self.galID == self.bcgid:
            width = 1
        else: 
            width = 0.2
        result = None

        # Loop to encompass all the particles in the domain
        while result is None:
            try:
                    
                    bbox = np.array([[-width,width], [-width, width], [-width, width]])
                    ds_all = yt.load_particles(data_all, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
                    result = yt.ParticleProjectionPlot(ds_all,'x',("star","particle_mass"))
            except:        
                    width=width+0.2
        
        # bounding box final
        bbox = np.array([[-width,width], [-width, width], [-width, width]])

        # load the dataset and return
        ds_all = yt.load_particles(data_all, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
        
        return ds_all
    

    def _half_mass_radius(self, particle_type):
        # Determine which attributes to use based on the particle type
        if particle_type == 'star':
            pos = self.star_pos[:]
            mass = self.star_mass[:]
            total_mass = self.gal_mstar
        elif particle_type == 'gas':
            pos = self.gas_pos[:]
            mass = self.gas_mass[:]
            total_mass = self.gal_mgas
        elif particle_type == 'dm':
            pos = self.dm_pos[:]
            mass = self.dm_mass[:]
            total_mass = self.gal_mdm
        else:
            raise ValueError("Invalid particle type. Choose from 'star', 'gas', or 'dm'.")

        # Calculate the center of mass (COM) for the chosen particle type
        com = np.average(pos, weights=mass, axis=0)

        # Calculate the distance from COM for each particle
        pos_gal = pos - com
        r_sc = np.linalg.norm(pos_gal, axis=1)


        # Store the r_sc in the class as an attribute based on particle type
        if particle_type == 'star':
            self.r_star_sc = r_sc
        elif particle_type == 'gas':
            self.r_gas_sc = r_sc
        elif particle_type == 'dm':
            self.r_dm_sc = r_sc

        # Create a DataFrame with distances and masses
        gal_data = pd.DataFrame({'rcom': r_sc, 'mass': mass})

        # Sort the DataFrame by the distance from the center
        df_sorted = gal_data.sort_values(by='rcom').reset_index(drop=True)

        # Calculate the cumulative mass
        df_sorted['cumulative_mass'] = df_sorted['mass'].cumsum()

        # Find the half-mass value
        half_mass = total_mass / 2

        # Find the index where the cumulative mass first reaches or exceeds the half mass
        half_mass_index = df_sorted[df_sorted['cumulative_mass'] >= half_mass].index[0]

        # Set the half-mass radius
        half_mass_radius = df_sorted.loc[half_mass_index, 'rcom']

        return half_mass_radius



    
    
    # function for metallicity gradient
    def get_metal_slope(self, r_rhalf_max=2, r_bin_width=0.1, var=None):
        """
        Parameters
        ----------
        r_rhalf_max: float
            Maximum radius of the bin
        r_bin_width: float
            Width of the bin
        var: str
            The variable to be used can be either 'star' or 'gas'
        """

        # Get half-mass radius and select metallicity and distance arrays based on particle type
        half_mass_radius = self._half_mass_radius(var)
        if var == 'star':
            met, dist = self.star_z[:], self.r_star_sc
        elif var == 'gas':
            met, dist = self.gas_z[:], self.r_gas_sc
        else:
            raise ValueError("Invalid particle type. Choose from 'star', 'gas', or 'dm'.")

        # Normalize metallicity and calculate r/r_half
        gal_data = pd.DataFrame({
            'r_rhalf': dist / half_mass_radius,
            var: met / 0.02
        })

        # Bin the data
        bins = np.arange(0, r_rhalf_max, r_bin_width)
        gal_data['binned'] = pd.cut(gal_data['r_rhalf'], bins, labels=False)

        # Calculate median r/r_half and metallicity for each bin
        binned_data = gal_data.groupby('binned').agg({'r_rhalf': 'median', var: 'median'}).dropna()

        # Perform linear regression
        slope, _, _, _, std_err = stats.linregress(binned_data['r_rhalf'], binned_data[var])

        return binned_data['r_rhalf'].values, binned_data[var].values, slope, std_err


# function for metallicity gradient
    def get_feh_slope(self, r_rhalf_max=2, r_bin_width=0.1):
        """
        Parameters
        ----------
        r_rhalf_max: float
            Maximum radius of the bin
        r_bin_width: float
            Width of the bin
        """

        # Get half-mass radius and select metallicity and distance arrays based on particle type
        half_mass_radius = self._half_mass_radius('gas')
        feh, dist = self.gas_fe[:]/self.gas_h[:], self.r_gas_sc
        
        if half_mass_radius>0:
            # Normalize metallicity and calculate r/r_half
            gal_data = pd.DataFrame({
                'r_rhalf': dist / half_mass_radius,
                'feh': np.log10(feh) + 4.5
            })

            
            # Bin the data
            bins = np.arange(0, r_rhalf_max, r_bin_width)

            gal_data['binned'] = pd.cut(gal_data['r_rhalf'], bins, labels=False)

            # Calculate median r/r_half and metallicity for each bin
            binned_data = gal_data.groupby('binned').agg({'r_rhalf': 'median', 'feh': 'median'}).dropna()

            
            
            slope, _, _, _, std_err = stats.linregress(binned_data['r_rhalf'], binned_data['feh'])
            
        
            return binned_data['r_rhalf'].values, binned_data['feh'].values, slope, std_err

        else:
             return np.nan, np.nan, np.nan, np.nan


class Analysis:
    def __init__(self, snap):
        self.snap = snap
        self.slope_df = None 
        self.median_gradient_met = None
        self.median_gradient_feh = None
        self.median_slope_Zg = None
        self.median_slope_Zs = None
    
    def get_slope_data(self,galids=None,clusids=None,rmax=4,rbin_width=0.1,var='met',dump='False'):
        """
        Parameters
        ----------
        galids: list
            List of galaxy IDs
        clusids: list
            List of cluster IDs corresponding to the galaxy IDs
        var: str
            The variable to be used can be either 'met' or 'feh'
        dump: str
            Whether to dump the data or not
        plot: str
            Whether to plot the data or not

        """
        
        
        # if data is not to be dumped load from {outdir}
        if dump=='False':
            main_df = pd.read_json(f'{outdir}/Gradient_{var}.json', orient='lines')

            # if feh is to be plotted
            if var == 'feh':         
                # Save DataFrame to JSON file
                slope_df_with_gas = main_df[main_df['R_g'].notna()]
                median_Rg = median_list_column(slope_df_with_gas['R_g'])
                median_feh = median_list_column(slope_df_with_gas['feh'])

                
                self.median_gradient_feh = pd.DataFrame({'median_Rg':median_Rg,'median_feh':median_feh})

                

            # if metallicity is to be plotted
            if var == 'met':
                                
                                
                # Compute the median list for columns R and Z
                median_Rs = median_list_column(main_df['R_s'])
                median_Zs = median_list_column(main_df['Z_s'])

                slope_df_with_gas = main_df[main_df['R_g'].notna()]
                median_Rg = median_list_column(slope_df_with_gas['R_g'])
                median_Zg = median_list_column(slope_df_with_gas['Z_g'])

                self.median_gradient_met = pd.DataFrame({'median_Rg':median_Rg,'median_Zg':median_Zg,'median_Rs':median_Rs,'median_Zs':median_Zs})


            

        elif dump=='True':

            # setup main dataframe that will contain data
            if var == 'met':
                main_df = pd.DataFrame(columns=['ID','R_s','Z_s','slope_s','std_e_s','R_g','Z_g','slope_g','std_e_g'])
            elif var=='feh':
                main_df = pd.DataFrame(columns=['ID','R_s','feh','slope_feh','std_e_feh'])


            # Loop over the galaxy IDs and cluster IDs
            for galid,clusid in tqdm.tqdm(zip(galids,clusids),total=len(galids)):
            
                    # initiate cluster object and galaxy object
                    clus = Cluster(self.snap,clusid)
                    gal = clus.get_alldat_gal(galid)
                    
                    # check to see required variable is 'feh' or 'met'
                    if var == 'feh':

                        # get metallicity gradient only for galaxy with gas
                        if gal.gal_mgas>0:
                            R_g,feh,slope_feh,std_e_feh = gal.get_feh_slope(r_rhalf_max=rmax,r_bin_width=rbin_width)
                        else:
                            R_g,feh,slope_feh,std_e_feh = np.nan,np.nan,np.nan,np.nan

                        # append the main_df
                        main_df = main_df.append({'ID':galid,'clusID':clusid,'R_g':R_g,'feh':feh,'slope_feh':slope_feh,'std_e_feh':std_e_feh},ignore_index=True)

                        
                    elif var == 'met':
                        
                        # get metallicity gradient only for galaxies
                        R_s,Z_s,slope_s,std_e_s = gal.get_metal_slope(r_rhalf_max=rmax,r_bin_width=rbin_width,var='star')

                        # get metallicity gradient only for galaxy with gas
                        if gal.gal_mgas>0:
                            R_g,Z_g,slope_g,std_e_g = gal.get_metal_slope(r_rhalf_max=4,r_bin_width=0.3,var='gas')
                        else:
                            R_g,Z_g,slope_g,std_e_g = np.nan,np.nan,np.nan,np.nan
                        
                       
                        # append the main_df
                        main_df = main_df.append({'ID':galid,'clusID':clusid,'R_s':R_s,'Z_s':Z_s,'slope_s':slope_s,'std_e_s':std_e_s,'R_g':R_g,'Z_g':Z_g,'slope_g':slope_g,'std_e_g':std_e_g},ignore_index=True)
            


            # if feh is to be plotted
            if var == 'feh':         
                # Save DataFrame to JSON file
                slope_df_with_gas = main_df[main_df['R_g'].notna()]
                median_Rg = median_list_column(slope_df_with_gas['R_g'])
                median_feh = median_list_column(slope_df_with_gas['feh'])

                
                self.median_gradient_feh = pd.DataFrame({'median_Rg':median_Rg,'median_feh':median_feh})

                

            # if metallicity is to be plotted
            if var == 'met':
                                
                                
                # Compute the median list for columns R and Z
                median_Rs = median_list_column(main_df['R_s'])
                median_Zs = median_list_column(main_df['Z_s'])

                slope_df_with_gas = main_df[main_df['R_g'].notna()]
                median_Rg = median_list_column(slope_df_with_gas['R_g'])
                median_Zg = median_list_column(slope_df_with_gas['Z_g'])

                self.median_gradient_met = pd.DataFrame({'median_Rg':median_Rg,'median_Zg':median_Zg,'median_Rs':median_Rs,'median_Zs':median_Zs})

            
            main_df.to_json(f'{outdir}/Gradient_{var}.json',orient='records')
            
        self.slope_df = main_df

        if var == 'feh':
            self.median_slope_feh = main_df['slope_feh'].median()
        elif var == 'met':
            self.median_slope_Zs = main_df['slope_s'].median()
            self.median_slope_Zg = main_df['slope_g'].median()
            

        
        
    def plot_gradient(self,var=None,save=True):
        
        if var == 'met':
                fig,ax = plt.subplots(2,1,figsize=(8,10),sharex=True)

        elif var == 'feh':
                fig,ax = plt.subplots(1,1)

        # plot for metallicity
        if var == 'met':
                    for i in tqdm.tqdm(range(len(self.slope_df))):
    
                        R_g = self.slope_df['R_g'].iloc[i]
                        R_s = self.slope_df['R_s'].iloc[i]
                        Z_s = self.slope_df['Z_s'].iloc[i]
                        Z_g = self.slope_df['Z_g'].iloc[i]
                        
                        
                        ax[0].plot(R_s,np.log10(Z_s),color='grey',alpha=0.1)
                        
                        # check if R_g list is not nan

                        if isinstance(R_g, np.ndarray) or isinstance(R_g, list):
                            ax[1].plot(R_g, np.log10(Z_g), color='grey', alpha=0.1)


                    # set labels and plot
                    ax[1].set_xlabel(r'$\mathrm{R/R_{\mathrm{half}}}$')
                    ax[0].set_ylabel(r'$\mathrm{log \ Z_{star}(Zsun)}$')
                    ax[1].set_ylabel(r'$\mathrm{log \ Z_{gas}(Zsun)}$')

                    ax[0].plot(self.median_gradient_met.median_Rs,np.log10(self.median_gradient_met.median_Zs),marker='o',markersize=8,label=r'$\mathrm{Median \ Z_{star}}$')

                    
                    ax[1].plot(self.median_gradient_met.median_Rg,np.log10(self.median_gradient_met.median_Zg),marker='o',markersize=8,label=r'$\mathrm{Median \ Z_{gas}}$')
                    
                    ax[0].legend()
                    ax[1].legend()

                    # ax[1].set_ylim(-0.5,1)
                    ax[1].set_xlim(0,4)

                    fig.subplots_adjust(hspace=0)

                    
        elif var == 'feh':
                    # plot for feh
                    for i in tqdm.tqdm(range(len(self.slope_df))):

                        R_g = self.slope_df['R_g'].iloc[i]
                        feh = self.slope_df['feh'].iloc[i]

                        # check if R_g list or numpy array
                        if isinstance(R_g, np.ndarray) or isinstance(R_g, list):
                            ax.plot(R_g, feh, color='grey', alpha=0.1)



                    ax.plot(self.median_gradient_feh.median_Rg,self.median_gradient_feh.median_feh,marker='o',markersize=8,label='MedianFe/H')

                    ax.set_xlabel(r'$\mathrm{R/R_{\mathrm{half}}}$')
                    ax.set_ylabel(r'$\mathrm{[F/H]}$')


        if save:    

            fig.savefig(f'../Plots/Gradient_{var}.png')
            fig.savefig(f'../Plots/Gradient_{var}.pdf')



    