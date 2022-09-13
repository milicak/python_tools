#!/bin/bash
#SBATCH -A gooogr 
#SBATCH -J "SOSE" 
#SBATCH --time=24:00:00
##SBATCH -p longq
##SBATCH -p shortq
##SBATCH -p gpuq
##SBATCH -p bigmemq
##SBATCH -p mixq
##SBATCH -p defq
##SBATCH -p b224q
#SBATCH -p core40q
#SBATCH -n 40
#SBATCH --output=%j.out
#SBATCH --error=%j.err


source ~/loadmodule_python 
python subset_calc_era5.py
