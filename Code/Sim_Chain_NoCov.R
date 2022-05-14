library(LTNLDA)
library(phyloseq)
library(topicmodels)

#load the simulated dataset
load("../Data/ps_sim.rda")

#set the number of subcommunities according to slurm task id
K = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#set C to truth
C = 9

#run the LTNLDA model - default values correspond to true values except for K and C
results = LTNLDA(ps_sim,K = K, C = C)
save(results,file = paste("Sim_C9_TrueK4_K",toString(K),".rda"))