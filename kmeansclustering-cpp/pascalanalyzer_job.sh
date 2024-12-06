#!/bin/bash
#SBATCH --job-name=PaScal_job
#SBATCH --output=PaScal_job%j.out
#SBATCH --error=PaScal_job%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-0:10
#SBATCH --partition=amd-512

pascalanalyzer -t aut -c 1,2,4,8,16,32,64,128 -i 7,13,25,50,100,200,400,800 -o main_read_txt_my_vector_no_parallelized.json ./main_read_txt


