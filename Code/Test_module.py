import HR5_cluster as CL
import re 
import yt 

clus = CL.Cluster(296,1561636)

#gal = clus.get_alldat_gal([1613370])
# print(gal[0].star_pos[:,0])
#gal = clus.get_all_parts('gas')
# gal = clus.get_all_parts('dm')

# print(gal._rcom_star)
#print(gal.rcom_gas)
# print(gal._rcom_dm)

ds_all,ds_rest,ds_bcg,ds_icm = clus.save_yt_dataset(1561636)
proj = yt.ParticleProjectionPlot(ds_all,'x',("star","particle_mass"))
