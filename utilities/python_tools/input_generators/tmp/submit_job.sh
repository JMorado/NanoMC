#!/bin/bash

#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH -J mc
#SBATCH -p batch


./nanomc_uvt.exe < input_file > output
