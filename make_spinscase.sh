# .bash

cd systems
# Load modules
module load slurm/15.08.6
module load intel/mkl/64/11.3.4/2016.4.258 openmpi/intel-opa/intel-hfi/64/1.10.4 intel/compiler/64/16.0.4/2016.4.258

./makemake.sh oswald

cd ../
./make_deps.sh -j oswald

cd src
make cases/dipole/dipole.x
