#!/bin/bash
#
#SBATCH --nodes=1
#
#SBATCH --ntasks=1
#SBATCH --partition=astro_long 
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-11:59:59
#SBATCH --cpus-per-task=20
#SBATCH --mem=32768
srun ./mu1e5-100by100by300.out > mu1e5-100by100by300.raw 

