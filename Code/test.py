import pandas as pd
import HR5_cluster as hr5
import os 
import matplotlib.pyplot as plt
plt.style.use('../paper_style.mplstyle')
import numpy as np
import configparser
# load up the parameter file
parser = configparser.ConfigParser()
parser.read('../params_hr5.ini')

outdir = parser.get('Paths','outdir')
galcats = parser.get('Paths','galaxycats')
morphs = parser.get('Paths','morphs')
snapdir = parser.get('Paths','snapfiles')
# get snaps
snaps = [ int(os.path.splitext(file)[0]) for file in os.listdir(snapdir)if  int(os.path.splitext(file)[0])>50]

fig,ax = plt.subplots(1,1)

vars = {'feh':r'$\rm [Fe/H]$','Zs':r'$\rm Z_{\rm star}$','Zg':r'$\rm Z_{\rm gas}$'}

# define three distinct colors for four curves
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

 # get snapshot data
time_df = pd.read_csv('../Data/Time_data.csv')
time_df.columns = ['snapshot','redshift','LBT','dx']
time_df.set_index('snapshot',inplace=True)

# redshift you are interested in
redshifts = [0.625,1,2,3,5]

# get values of snapshots near these redshift
snaps= [(time_df['redshift']-red).abs().idxmin()for red in redshifts]

redshifts = time_df.loc[snaps]['redshift'].tolist()


for i,var in enumerate(['sfr']):#enumerate(vars.keys()):
    
  

    for j,snap in enumerate(snaps):#sorted(snaps):
        out = f'{outdir}/Slope_{snap}.json'

        snapm = str(snap).zfill(3)
        # get only elliptical galaxies
        morpho = pd.read_csv(f'{morphs}/morphology_{snapm}.csv')

        e_galaxies = morpho[morpho['sersicn']>2.5]

        red = redshifts[j]

        Ana = hr5.Analysis(snap)

        # chech if the json file exists
        if not os.path.exists(f'{outdir}/Gradient_{var}_{snap}.json'):
                Ana.get_slope_data(galids=e_galaxies.ID,clusids=e_galaxies.clusID,rmax=4,rbin_width=0.3,var=f'{var}',dump_data=True,use_cache=False)
        else:
               Ana.get_slope_data(galids=e_galaxies.ID,clusids=e_galaxies.clusID,rmax=4,rbin_width=0.3,var=f'{var}',dump_data=False,use_cache=True)
        
        grad = getattr(Ana, f'median_gradient_{var}')
            
        ax.plot(grad['median_Rs'],grad[f'median_{var}'],marker='o',markersize=8,label=f'{red:0.3f}',color=colors[j])

        ax.set_ylabel(f'{vars[var]}')


        ax.set_xlabel(r'$\mathrm{R/R_{\mathrm{half}}}$')
            
        
        


ax.legend(loc=3,title='Redshift')      
# fig.subplots_adjust(hspace=0)
# fig.savefig('../Plots/feh_evolution_redshift.png',bbox_inches='tight')      


    
