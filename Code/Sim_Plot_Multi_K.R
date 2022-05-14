library(ggplot2)
library(phyloseq)

#TRy 4,7,10?

set.seed(1)

####################################
# Set up long form data structures #
####################################

props = NULL
samps = NULL
tops = NULL
model = NULL
Subcommunities = NULL

#########
# K = 4 #
#########

#assume this data structure isn't in the package yet
#load blocked covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_Cov_C9_TrueK4_K 4 .rda")

post_phi_dk = results$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_kv = results$Mean_Post_Beta_k

#load independent covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_C9_TrueK4_K 4 .rda")
post_phi_ind_dk = results$Mean_Post_Phi_d

post_beta_ind_kv = results$Mean_Post_Beta_k

#load LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_cov_LDA_K_ 4 .rda")

#correct label switching for lda
dat_perm = c(2,3,1,4)

#correct label switching for ind ltn
ind_perm = c(4,2,3,1)

#permutation for K = 4
top_perm = 1:K
for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LTN LDA Cov")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,post_phi_ind_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,ind_perm[k])
    model = c(model,"LTN LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,lda_post$topics[d,k])
    samps = c(samps,d)
    tops = c(tops,dat_perm[k])
    model = c(model,"LDA")
    Subcommunities = c(Subcommunities,toString(K))
  }
}

#########
# K = 5 #
#########

#assume this data structure isn't in the package yet
#load blocked covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_Cov_C9_TrueK4_K 5 .rda")

post_phi_dk = results$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_kv = results$Mean_Post_Beta_k

#load independent covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_C9_TrueK4_K 5 .rda")
post_phi_ind_dk = results$Mean_Post_Phi_d

post_beta_ind_kv = results$Mean_Post_Beta_k

#load LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_cov_LDA_K_ 5 .rda")

#correct label switching for lda
dat_perm = c(3,2,1,4,5)

#correct label switching for ind ltn
ind_perm = c(3,1,5,4,2)

#permutation for K = 4
top_perm = c(1,4,3,2,5)

for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LTN LDA Cov")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,post_phi_ind_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,ind_perm[k])
    model = c(model,"LTN LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,lda_post$topics[d,k])
    samps = c(samps,d)
    tops = c(tops,dat_perm[k])
    model = c(model,"LDA")
    Subcommunities = c(Subcommunities,toString(K))
  }
}

#########
# K = 7 #
#########

#assume this data structure isn't in the package yet
#load blocked covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_Cov_C9_TrueK4_K 7 .rda")

post_phi_dk = results$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_kv = results$Mean_Post_Beta_k

#load independent covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_C9_TrueK4_K 7 .rda")
post_phi_ind_dk = results$Mean_Post_Phi_d

post_beta_ind_kv = results$Mean_Post_Beta_k

#load LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_cov_LDA_K_ 7 .rda")

#correct label switching for lda
dat_perm = c(1,4,2,6,7,5,3)

#correct label switching for ind ltn
ind_perm = c(5,6,3,4,1,2,7)

#permutation for K = 4
top_perm = c(4,2,3,1,7,6,5)

for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LTN LDA Cov")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,post_phi_ind_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,ind_perm[k])
    model = c(model,"LTN LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,lda_post$topics[d,k])
    samps = c(samps,d)
    tops = c(tops,dat_perm[k])
    model = c(model,"LDA")
    Subcommunities = c(Subcommunities,toString(K))
  }
}

##########
# K = 10 #
##########

#assume this data structure isn't in the package yet
#load blocked covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_Cov_C9_TrueK4_K 10 .rda")

post_phi_dk = results$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]

post_beta_kv = results$Mean_Post_Beta_k

#load independent covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_C9_TrueK4_K 10 .rda")
post_phi_ind_dk = results$Mean_Post_Phi_d

post_beta_ind_kv = results$Mean_Post_Beta_k

#load LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_cov_LDA_K_ 10 .rda")

#correct label switching for lda
dat_perm = c(6,7,10,9,8,5,4,1,2,3)

#correct label switching for ind ltn
ind_perm = c(2,3,1,8,10,9,7,4,6,5)

#permutation for K = 4
top_perm = c(10,7,4,9,3,6,2,8,5,1)

for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LTN LDA Cov")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,post_phi_ind_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,ind_perm[k])
    model = c(model,"LTN LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,lda_post$topics[d,k])
    samps = c(samps,d)
    tops = c(tops,dat_perm[k])
    model = c(model,"LDA")
    Subcommunities = c(Subcommunities,toString(K))
  }
}

###########################
# Make dataframe and plot #
###########################

df = data.frame(props,samps,tops,model,Subcommunities)

top_perm = 1:K
for(i in 1:nrow(df)){
  df$tops[i] = top_perm[df$tops[i]]
}

p1 = ggplot(df, aes(x=samps, y=props)) +
  geom_line(aes(col = Subcommunities)) +
  scale_color_manual( values = c("firebrick","darkgreen","blue","black")) +
  facet_grid(tops ~ model, scale = "free_x") +
  xlab("") +
  ylab("") +
  ggtitle("Posterior Mean Estimates of ASV-Subcommunity Distributions") + 
  theme_classic()
p1