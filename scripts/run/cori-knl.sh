#!/bin/bash

#SBATCH --qos=debug
#SBATCH --time=0:05:00
#SBATCH -N1 -n64
#SBATCH --constraint=knl
#SBATCH -J test_job
#SBATCH --error job.err
#SBATCH --output job.out

export CRAY_CPU_TARGET=mic-knl
module load PrgEnv-intel
module load boost

srun -n 64 /global/homes/t/terzani/Pic/forked_ALaDyn/bin/ALaDyn >> opic.txt 2>> epic.txt

