###########################################
# Run FastPG on docker 
# log in username@login.bioinf.unc.edu
# get a interactive node, and run bash script to get into the docker environment
# srun -N 1 --cpus-per-task 8 --pty bash
#singularity shell -B $PWD --pwd $PWD -C fastpg-latest.simg
# then run R

input_file <- "Zhengmix4eq.2k.30PC.list"
input_file <- "Zhengmix8eq.2k.30PC.list"
input_file <- "pbmc68k.3k.30PC.list"

# import data
dat <- readRDS(file=file.path("data", paste(input_file, ".rds", sep="")))

# set parameters for FastPG
k <- 20
num_threads <- 8

# run FastPG