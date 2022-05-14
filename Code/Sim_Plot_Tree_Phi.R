library(ggplot2)
library(phyloseq)


set.seed(1)

####################################
# Set up long form data structures #
####################################

props = NULL
samps = NULL
tops = NULL
mod = NULL
Subcommunities = NULL

#########
# K = 4 #
#########

#assume this data structure isn't in the package yet
#load regular LTN-LDA
load("../../Simulations/Response/Sim_C10_TrueK4_K 4 .rda")

post_phi_dk = model$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_kv = model$Mean_Post_Beta_k

#load independent covariance LTN-LDA
load("../../Simulations/Response/Sim_Unif_Tree_C10_TrueK4_K 4 .rda")
post_phi_unif_dk = model$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_unif_kv = model$Mean_Post_Beta_k

#correct label switching 
dat_perm = c(4,2,1,3)

#permutation for K = 4
top_perm = 1:K
for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    mod = c(mod,"LTN-LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,post_phi_unif_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,dat_perm[k])
    mod = c(mod,"LTN-LDA: Uniform Tree")
    Subcommunities = c(Subcommunities,toString(K))

  }
}

#########
# K = 5 #
#########

#assume this data structure isn't in the package yet
#load regular LTN-LDA
load("../../Simulations/Response/Sim_C10_TrueK4_K 5 .rda")

post_phi_dk = model$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_kv = model$Mean_Post_Beta_k

#load independent covariance LTN-LDA
load("../../Simulations/Response/Sim_Unif_Tree_C10_TrueK4_K 5 .rda")
post_phi_unif_dk = model$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_unif_kv = model$Mean_Post_Beta_k

#permutation for K = 4
top_perm = c(4,5,1,2,3)

#correct label switching
dat_perm = c(1,2,5,3,4)


for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    mod = c(mod,"LTN-LDA")
    Subcommunities = c(Subcommunities,toString(K))

    props = c(props,post_phi_unif_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,dat_perm[k])
    mod = c(mod,"LTN-LDA: Uniform Tree")
    Subcommunities = c(Subcommunities,toString(K))

  }
}

#########
# K = 7 #
#########

#assume this data structure isn't in the package yet
#load regular LTN-LDA
load("../../Simulations/Response/Sim_C10_TrueK4_K 7 .rda")

post_phi_dk = model$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_kv = model$Mean_Post_Beta_k

#load independent covariance LTN-LDA
load("../../Simulations/Response/Sim_Unif_Tree_C10_TrueK4_K 7 .rda")
post_phi_unif_dk = model$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_unif_kv = model$Mean_Post_Beta_k

#permutation for K = 4
top_perm = c(5,7,4,6,2,3,1)

#correct label switching
dat_perm = c(3,5,1,6,4,2,7)


for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    mod = c(mod,"LTN-LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,post_phi_unif_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,dat_perm[k])
    mod = c(mod,"LTN-LDA: Uniform Tree")
    Subcommunities = c(Subcommunities,toString(K))
    
  }
}

#########
# K = 10 #
#########

#assume this data structure isn't in the package yet
#load regular LTN-LDA
load("../../Simulations/Response/Sim_C10_TrueK4_K 10 .rda")

post_phi_dk = model$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_kv = model$Mean_Post_Beta_k

#load independent covariance LTN-LDA
load("../../Simulations/Response/Sim_Unif_Tree_C10_TrueK4_K 10 .rda")
post_phi_unif_dk = model$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_unif_kv = model$Mean_Post_Beta_k

#permutation for K = 4
top_perm = c(7,8,3,9,5,6,1,2,4,10)

#correct label switching
dat_perm = c(9,6,2,1,5,3,7,8,4,10)


for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    mod = c(mod,"LTN-LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,post_phi_unif_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,dat_perm[k])
    mod = c(mod,"LTN-LDA: Uniform Tree")
    Subcommunities = c(Subcommunities,toString(K))
    
  }
}

###########################
# Make dataframe and plot #
###########################

df = data.frame(props,samps,tops,mod,Subcommunities)

top_perm = 1:K
for(i in 1:nrow(df)){
  df$tops[i] = top_perm[df$tops[i]]
}

p1 = ggplot(df, aes(x=samps, y=props)) +
  geom_line(aes(col = Subcommunities)) +
  scale_color_manual( values = c("firebrick","darkgreen","blue","black")) +
  facet_grid(tops ~ mod, scale = "free_x") +
  xlab("") +
  ylab("") +
  ggtitle("Posterior Mean Estimates of ASV-Subcommunity Distributions") + 
  theme_classic()
p1
