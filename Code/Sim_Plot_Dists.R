library(gridExtra)
library(ggplot2)
library(phyloseq)


set.seed(1)

####################################
# Set up long form data structures #
####################################

#long-form data structures
vals = NULL #contains values of distributions 
input = NULL #contains the leaf number
samp = NULL #contains the sample number
Distribution = NULL #contains a distributoin id
tops = NULL #which subcommunity, corrects for label switching

########
K = 4 #
########

#assume this data structure isn't in the package yet
#load blocked covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_Cov_C9_TrueK4_K 5 .rda")

post_phi_dk = results$Mean_Post_Phi_d
K = dim(post_phi_dk)[2]
D = dim(post_phi_dk)[1]
V = dim(post_beta_kv)[2]

post_beta_kv = results$Mean_Post_Beta_k

post_beta_kdv = results$Mean_Post_Beta_kd


#load independent covariance LTN-LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_C9_TrueK4_K 4 .rda")
post_phi_ind_dk = results$Mean_Post_Phi_d
K = dim(post_phi_ind_dk)[2]
D = dim(post_phi_ind_dk)[1]
V = dim(post_beta_ind_kv)[2]

post_beta_ind_kv = results$Mean_Post_Beta_k

post_beta_ind_kdv = results$Mean_Post_Beta_kd

#correct label switching for ind ltn
ind_perm = 1:K

#load LDA
load("../../Simulations/Blocked_Covariance/Case_Studies/Sim_cov_LDA_K_ 4 .rda")

#correct label switching for lda
# #k = 4
dat_perm = 1:K

top_perm = 1:K
for(k in 1:K){
  # ord = order(post_beta_kv[k,])
  
  for(d in 1:D){
    vals = c(vals,post_beta_kdv[k,d,])
    input = c(input,1:V)
    samp = c(samp,rep(d,V))
    Distribution = c( Distribution,rep("LTN LDA Samples",V))
    tops = c(tops,rep(top_perm[k],V))
  }
  
  # vals = c(vals,lda_post$terms[dat_perm[k],])
  vals = c(vals,rep(0,V))
  input = c(input,1:V)
  samp = c(samp,rep(D + 1,V))
  Distribution = c( Distribution,rep("LDA Average",V))
  tops = c(tops,rep(k,V))
  
  vals = c(vals,post_beta_kv[k,])
  input = c(input,1:V)
  samp = c(samp,rep(D + 2,V))
  Distribution = c( Distribution,rep("LTN LDA Average",V))
  tops = c(tops,rep(top_perm[k],V))
  
}

#make dataframe
Sim_dist_df = data.frame(input,vals,samp, Distribution,tops)

#make plots of sample-to-sample heterogeneity in distributions
p1 = ggplot(Sim_dist_df, aes(x=input, y=vals, group=samp, col = Distribution)) +
  geom_line(aes(size = Distribution,col = Distribution)) +
  scale_color_manual( values = c("black","firebrick","skyblue1")) +
  scale_size_manual( values = c(1,1,0.1)) +
  facet_wrap(Sim_dist_df$tops,ncol = 1) +
  xlab("") +
  ylab("") +
  ggtitle("Posterior Mean Estimates of ASV-Subcommunity Distributions for Cov LTNLDA") + 
  theme_classic()


####################################
# Set up long form data structures #
####################################

#long-form data structures
vals = NULL #contains values of distributions 
input = NULL #contains the leaf number
samp = NULL #contains the sample number
Distribution = NULL #contains a distributoin id
tops = NULL #which subcommunity, corrects for label switching

for(k in 1:K){
  # ord = order(post_beta_kv[k,])
  
  for(d in 1:D){
    vals = c(vals,post_beta_ind_kdv[k,d,])
    input = c(input,1:V)
    samp = c(samp,rep(d,V))
    Distribution = c( Distribution,rep("LTN LDA Samples",V))
    tops = c(tops,rep(ind_perm[k],V))
  }
  
  # vals = c(vals,lda_post$terms[dat_perm[k],])
  vals = c(vals,rep(0,V))
  input = c(input,1:V)
  samp = c(samp,rep(D + 1,V))
  Distribution = c( Distribution,rep("LDA Average",V))
  tops = c(tops,rep(k,V))
  
  vals = c(vals,post_beta_ind_kv[k,])
  input = c(input,1:V)
  samp = c(samp,rep(D + 2,V))
  Distribution = c( Distribution,rep("LTN LDA Average",V))
  tops = c(tops,rep(ind_perm[k],V))
  
}

#make dataframe
Sim_dist_df = data.frame(input,vals,samp, Distribution,tops)

#make plots of sample-to-sample heterogeneity in distributions
p2 = ggplot(Sim_dist_df, aes(x=input, y=vals, group=samp, col = Distribution)) +
  geom_line(aes(size = Distribution,col = Distribution)) +
  scale_color_manual( values = c("black","firebrick","skyblue1")) +
  scale_size_manual( values = c(1,1,0.1)) +
  facet_wrap(Sim_dist_df$tops,ncol = 1) +
  xlab("") +
  ylab("") +
  ggtitle("Posterior Mean Estimates of ASV-Subcommunity Distributions for Ind LTNLDA") + 
  theme_classic()

#plot both plots side-by-side
grid.arrange(p1,p2,widths = c(1,1))