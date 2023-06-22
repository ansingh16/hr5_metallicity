# import HR5_cluster as CL
# import re 
# import yt 

# clus = CL.Cluster(296,1561636)

# gal = clus.get_alldat_gal(1613370)
# print(dir(gal))
# print(gal.gal_pos)
#gal = clus.get_all_parts('gas')
# gal = clus.get_all_parts('dm')

# print(gal._rcom_star)
#print(gal.rcom_gas)
# print(gal._rcom_dm)

# ds_all,ds_rest,ds_bcg,ds_icm = clus.save_yt_dataset(1561636)
# proj = yt.ParticleProjectionPlot(ds_all,'x',("star","particle_mass"))

from math import radians, sin, cos

def cpu_bench(number):
    product = 1.0
    for elem in range(number):
        angle = radians(elem)
        product *= sin(angle)**2 + cos(angle)**2
    return product

tasks = [1000000 + i for i in range(100)]



from parallelbar import progress_map

if __name__=='__main__':
    result = progress_map(cpu_bench, tasks,n_cpu=2,chunk_size=1,core_progress=True)
