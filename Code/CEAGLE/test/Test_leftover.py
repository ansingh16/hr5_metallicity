import h5py 
import glob 
import os 
import configparser



# load up the parameter file
parser = configparser.ConfigParser()
parser.read('/home/ankitsingh/hr5_metalicity/params.ini')

outdir = parser.get('Paths','outdir')
snapfiles = parser.get('Paths','snapfiles')
files = glob.glob(f'{snapfiles}*.dat')
snap_files = sorted([int(os.path.basename(snap).split('.')[0]) for snap in files],reverse=True)

print(snap_files)

tobepreocessed = []

#with open('../to_be_processed.txt', 'w') as out:
for snapno in snap_files:
        try:   
            with h5py.File(f"{outdir}/clusters{snapno}.hdf5", "r") as f:
                            
                if 'status' in f:
                    print(f'present in {snapno} with satus {f["status"][()]}')
                            
        except:
            #out.write(f'{snapno}\n')
            print(f'cant open {snapno}')
