###########################################
# Run FastPG on docker 
# log in username@login.bioinf.unc.edu
# get a interactive node, and run bash script to get into the docker environment
# srun -N 1 --cpus-per-task 8 --pty bash
#singularity shell -B $PWD --pwd $PWD -C fastpg-latest.simg
# then run R

input_file <- "Zhengmix4eq.30PC"
input_file <- "Zhengmix8eq.30PC"
input_file <- "pbmc68k.30PC"

# import data
dat <- readRDS(file=paste(input_file, ".rds", sep=""))

# set parameters for FastPG
k <- 20
num_threads <- 8

# run FastPG
fg.clusters <- FastPG::fastCluster( dat, k, num_threads )

# extract cell assignment
fg_cl <- fg.clusters$communities

# save cell assignment
saveRDS(fg_cl, file=paste(input_file, ".fg.cl", sep=""))




