import pandas as pd
import numpy as np 


morpho = pd.read_csv('../Data/morphology.csv')


dtype = {'ID':np.int32, 'mstar': np.float32,'galstarmass': np.float32, 'rmssize':np.float32, 'asym':np.float32,'einasn':np.float16, 'coni':np.float32, 'rhalf':np.float32, 'sersicn':np.float16, 'vrot':np.float32, 'vsig':np.float32, 'sfr':np.float32}


keys = list(dtype.keys())
values = list(dtype.values())

morph = np.genfromtxt('../Data/galmorf.out_296.txt',dtype=values)

morph = pd.DataFrame(morph)
morph.columns = keys

# value of sersicn from morph for galaxies with ID in morpho

morpho = morpho.merge(morph[['ID', 'sersicn']], on='ID', how='inner')

morpho.reset_index(drop=True,inplace=True)

morpho.drop(columns='Unnamed: 0',inplace=True)

morpho = pd.read_csv('../Data/morphology.csv')

