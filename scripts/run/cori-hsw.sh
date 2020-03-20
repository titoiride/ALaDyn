#!/bin/bash

#SBATCH --qos=debug
#SBATCH --time=0:05:00
#SBATCH -N20 -n320
#SBATCH --constraint=haswell
#SBATCH -J test_job
#SBATCH --error job.err
#SBATCH --output job.out

export CRAY_CPU_TARGET=haswell
module load PrgEnv-intel
module load boost

srun -n 320 ./ALaDyn >> opic.txt 2>> epic.txt

