import HR5_cluster as hr5
import pandas as pd
from multiprocessing import Pool, cpu_count
import tqdm

# Get all the IDs of clusters present at the given snapshot
snapshot = 296

# Define the instance of the class Cluster with the snapshot
clus296 = pd.read_csv('../Data/groups5e13.csv')

snap296 = pd.read_parquet('/scratch/ankitsingh/Galaxy_catalogs/galaxy_catalogue_296.parquet')

def process_cluster(clusid):
    print(f"Processing cluster {clusid}")
    morpho = []
    clus = hr5.Cluster(snapshot, clusid)
    for galid in clus.get_galids():
            try:
                gal = clus.get_alldat_gal(galid)
                morpho.append({'ID': gal.galID,'clusID': clusid, 'morpho': gal.get_morphology()})
            except Exception as e:
                gale = snap296.loc[snap296['ID'] == galid]
                print(f"Galaxy {galid} in {clusid} with Mstar{gale['Mstar(Msun)'].values[0]} not found")

    return morpho

if __name__ == '__main__':
    # Use all available CPU cores
    with Pool(4) as pool:
        results = list(tqdm.tqdm(pool.imap(process_cluster, clus296.HostHaloID.values), total=len(clus296)))

    # Flatten the list of lists
    morpho = [item for sublist in results for item in sublist]

    # Convert to DataFrame
    morpho = pd.DataFrame(morpho)

    # Save to CSV
    morpho.to_csv('../Data/morphology.csv', index=False)
