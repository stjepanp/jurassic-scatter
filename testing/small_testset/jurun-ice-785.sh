#!/bin/bash -x
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --time=1:00:00
#SBATCH --partition=gpus
#SBATCH --gres=gpu:1
#SBATCH --account=slmet

#module load Python/3.8.5
#module load OpenMPI/4.1.0rc1
#module load CUDA/11.3
#module load Nsight-Systems/2021.1.1
#module load Nsight-Compute/2020.3.0
#module load ParaStationMPI/5.4.7-1
#module load mpi-settings/CUDA
#module load GCC/9.3.0
#module load GCCcore/.10.3.0
#module load GCCcore/.9.3.0
#module load Intel/2020.2.254-GCC-9.3.0
#module load GSL/2.6 


module load Python/3.8.5
module load GCC/9.3.0
module load OpenMPI/4.1.0rc1
module load CUDA/11.3
module load Nsight-Systems/2021.1.1
module load Nsight-Compute/2020.3.0
module load ParaStationMPI/5.4.7-1
module load mpi-settings/CUDA
module load GCCcore/.10.3.0
module load GCCcore/.9.3.0
module load GSL/2.6 

export OMP_NUM_THREADS=16

n=$(< aux/submission_index)
m=$(( n + 1 ))
echo $m > aux/submission_index
echo last: $n, new: $m

### set up
#src=/p/home/jusers/pozgaj1/juwels/GSP-stjepanp/Project/Repos/jurassic-scatter/src
src=../../src
w1=785
w2=798

### run forward model
time srun $src/formod cloud-${w1}-${w2}.ctl obs33.tab atm.tab submissions/rad-${m}.tab AEROFILE aero.tab DIRLIST aux/first_002
