## One-time setup instructions for rocket.hpc
- Create a directory in which to store the code, typically in your home directory (`cd ~`)
- In that directory, download this copy of the code: 
```
git clone https://github.com/HartharnSam/SPINS_main.git
chmod u=rwx ./*
```
From the SPINS\_main directory, run (line by line to check for errors):
```
module load intel/2021a Autoconf/2.71-GCCcore-10.3.0 Automake/1.16.3-GCCcore-10.3.0 libtool/2.4.6-GCCcore-10.3.0
patch --verbose -p1 < ./spins-autoconf.patch
cd systems
./makemake.sh rocket
cd ..
./make_deps.sh rocket
```
### Building your casefile
Then build your casefile as per the SPINS instructions:
- enter the src directory (`cd src`)
- `make cases/case_directory/yourcase.x
- This requires a file called your_case.cpp in the cases/case_directory directory. There are several cases included with the code so you may want to start with one of those. dipole.x is a good first test, and then something more complex with mode1_shoal.x . 
- To do a quick test, there should be an executable yourcase.x, move to the directory folder (`cd cases/case\_directory`) and execute the code to the first timestep (`./yourcase.x`), and then cancel with ctrl+c. 

### Running your case
- Be careful not to run in your home directory on rocket, cases typically fill up the 5GB quota and then cancel. Instead, copy your .x and spins.conf file to a new case directory in /nobackup/$USERNAME
- Set up your batch of simulations, each in a directory with an .x file and specific spins.conf file, and list them in a text file case_list
- Run the file you can find in SLURM_script/run_spins.sh:
` sbatch ./run\_spins.sh`:wq


 


