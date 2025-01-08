
import pandas as pd 

dat296 = pd.read_parquet('/scratch/ankitsingh/Galaxy_catalogs/galaxy_catalogue_296.parquet')

dat296 = dat296.loc[dat296['pure']==1]
datfil = dat296.loc[dat296['HostMtot(Msun)']>5.0e13]

datfil[['HostHaloID','HostMtot(Msun)']].drop_duplicates(keep='first').to_csv('./groups5e13.csv',index=False)
