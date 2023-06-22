#!/bin/sh
#PBS -N icl
#PBS -V
#PBS -q normal
#PBS -l nodes=1:ppn=10
#PBS -m abe
#PBS -r n
#PBS -A inhouse

module load gnu7

cd /home/ankitsingh/hr5_metallicity/Code/
ulimit -s unlimited

source ~/.bashrc
source /home/ankitsingh/miniconda3/bin/activate ICL

echo `date`
mpirun python /home/ankitsingh/hr5_metallicity/Code/ICL_work.py
echo `date`

