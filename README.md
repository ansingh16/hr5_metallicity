This repository contains the code for analysing the data extracted from [Horizon-Run 5 simulation](https://arxiv.org/abs/2006.01039).

The main code for analysis is present in **Code** directory. The jupyter file **HR5_metallicity.ipynb** contains a tutorial outlining the structure of the data files, methods and parameters available from the HR5_module available in the same folder.

Before starting the analysis following steps are recommended:

- Install [miniconda package management system](https://docs.conda.io/en/latest/miniconda.html) and create a python environment using the file icl_env.txt

        
        conda create --name hr5_metal --file icl_env.txt
        conda activate hr5_metal
- Change the paths given in the **params.ini** to user specific locations

For any querries plase feel free to contact *ankitsingh@kias.re.kr*
