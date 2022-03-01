library(LTNLDA)
library(phyloseq)
library(topicmodels)

#load the simulated dataset
load("../Data/ps_sim.rda")

#set the number of subcommunities according to slurm task id
K = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#set C to truth
C = 9

#run LDA on the dataset with the same K, for the same number of iterations
dtm = otu_table(ps_sim)
results_LDA = LDA(t(dtm),k = K, method = "Gibbs", 
                  control = list(burnin = 10000, thin = 10, iter = 1000))
lda_post = posterior(results_LDA)
save(lda_post,file = paste("Sim_cov_LDA_K_",toString(K),".rda"))

#run the LTNLDA model - default values correspond to true values except for K and C
results = LTNLDA_cov(ps_sim,K = K, C = C)
save(results,file = paste("Sim_Cov_C9_TrueK4_K",toString(K),".rda"))