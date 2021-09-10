#some of these are unnecessary
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(statmod)
library(matrixStats)
library(rmutil)
library(zoo)
library(pgdraw)
library(RcppDist)
library(BayesianGLasso)
library(topicmodels)
library(gridExtra)
library(randtests)


set.seed(1)

#Define Parameters
true_K = 4
D = 50
N = 10000
alpha = 1

#load an already existing edge matrix
load("../Data/tree_mat.rda")

#rename the data structure because of legacy code
toy_tree = tree.edge

#find characteristics of tree
#find the root node
root = setdiff(toy_tree[,1],toy_tree[,2]) 
#find the internal nodes and reorder
internal_nodes = unique(toy_tree[,1])
internal_nodes = sort(internal_nodes)
internal_nodes_C = internal_nodes - 1
#find A, the total number of nodes
A = max(internal_nodes)
#find the set of leaves
leaves = setdiff(toy_tree[,2],toy_tree[,1])
#find the number of leaves
V = length(leaves)

#Generate some tree data structures

#descendants[[i]] contains the two immediate descendants of node i
descendants = NULL
descendants_mat = matrix(0,ncol=2,nrow=max(tree.edge))
for(i in 1:max(tree.edge)){
  if (sum(which(tree.edge[,1]==i))>0){
    descendants[[i]] = tree.edge[which(tree.edge[,1] == i),2]
    descendants_mat[i,] = descendants[[i]]
  }
}
descendants_mat_C = descendants_mat - 1


#parents[i] contains the parent of node i
parents = NULL
for (i in 1:max(tree.edge)){
  if (sum(which(tree.edge[,2]==i))>0){
    parents[i] = tree.edge[which(tree.edge[,2]==i),1]
  }
}

#ancestors[[i]] contains all of the ancestors of node i
ancestors = NULL
ancestors_C = NULL
for (i in 1:max(tree.edge)){
  up=NULL
  parent = parents[i]
  while (is.na(parents[parent])==FALSE){
    up = c(up,parent)
    parent = parents[parent]
  }
  ancestors[[i]] = c(up,parent) #Adds the root of the tree as well
  ancestors_C[[i]] = ancestors[[i]]-1
}

#layers[[i]] containts the nodes in the i^th layer of the tree
#layers[[1]] is the root node, layers[[2]] is the root nodes children, etc
layers = NULL
#initialize layer 2 and 1
#have to do an n-tuple (n>1) first, o/w this explodes
layers[[2]] = descendants[[root]]
layers[[1]] = root
for (i in 3:max(tree.edge)){
  descend = NULL
  for (j in 1:length(layers[[i-1]])){
    descend = c(descend,descendants[[layers[[i-1]][j]]])
  }
  if ((sum(descend)>0)==TRUE){
    layers[[i]] = descend
  } else{
    break
  }
}

#left_leaves[[node]] contains the leaves left-descended from node
#right_leaves[[node]] contains the leaves right-descended from node
left_leaves = NULL
right_leaves = NULL
left_leaves[[max(internal_nodes)+1]] = rep(0,5)
right_leaves[[max(internal_nodes)+1]] = rep(0,5)
for (node in internal_nodes){
  left_descend = NULL
  right_descend = NULL
  
  descend = descendants[[node]]
  left = descend[1]
  right = descend[2]
  
  #if the descendant is a leaf we can termiante
  if((left %in% leaves)==TRUE){
    left_descend = left
  } else {
    #cycle through all of the leaves and see which are left descendants
    for (nodes in leaves){
      if((left %in% ancestors[[nodes]])==TRUE){
        left_descend = c(left_descend,nodes)
      }
    }
  }
  left_leaves[[node]] = left_descend
  
  #if the descendant is a leaf we can termiante
  if((right %in% leaves)==TRUE){
    right_descend = right
  } else {
    #cycle through all of the leaves and see which are right descendants
    for (nodes in leaves){
      if((right %in% ancestors[[nodes]])==TRUE){
        right_descend = c(right_descend,nodes)
      }
    }
  }
  right_leaves[[node]] = right_descend
}
left_leaves[[max(internal_nodes)+1]] = NULL
right_leaves[[max(internal_nodes)+1]] = NULL


#need to find, for each leaf
# the nodes which have to succeed
# the nodes which have to fail
#in order for the leaf to be selected

#leaf_success[[leaf]] contains the nodes from which leaf is left-descended
leaf_success = NULL
leaf_success_C = NULL
leaf_success[[max(leaves)+1]] = rep(0,5)
leaf_success_C[[max(leaves)+1]] = rep(0,5)
for (leaf in leaves){
  node_list=NULL
  node_list = c(leaf,ancestors[[leaf]])
  successes = NULL
  for (node in ancestors[[leaf]]){
    if ((descendants[[node]][1] %in% node_list)==TRUE){
      successes = c(successes,node)
    }
  }
  if (is.null(successes)==FALSE){
    leaf_success[[leaf]] = successes
    leaf_success_C[[leaf]] = successes - 1
  } else {
    leaf_success_C[[leaf]] = -1
  }
}
leaf_success[[max(leaves)+1]]  = NULL
leaf_success_C[[max(leaves)+1]]  = NULL

#leaf_failures[[leaf]] contains the nodes from which leaf is right-descended
leaf_failures = NULL
leaf_failures[[max(leaves)+1]] = rep(0,5)
leaf_failures_C = NULL
leaf_failures_C[[max(leaves)+1]] = rep(0,5)
for (leaf in leaves){
  node_list=NULL
  node_list = c(leaf,ancestors[[leaf]])
  failures = NULL
  for (node in ancestors[[leaf]]){
    if ((descendants[[node]][2] %in% node_list)==TRUE){
      failures = c(failures,node)
    }
  }
  if (is.null(failures)==FALSE){
    leaf_failures[[leaf]] = failures
    leaf_failures_C[[leaf]] = failures - 1
  } else {
    leaf_failures_C[[leaf]] = -1
  }
}
leaf_failures[[max(leaves)+1]] = NULL
leaf_failures_C[[max(leaves)+1]] = NULL



#come up with a mapping from internal nodes to 1:p
#node_map[internal_nodes] in {1,2,dots,p}
p = length(internal_nodes)
node_map = rep(0,A)
for(x in 1:p){
  node_map[internal_nodes[x]] = x
}

###################
#                 #
###################
# Data Generation #
###################
#                 #
###################

#useless data structure, part of legacy code.  Would have to change C code to get rid of
Phi = 0.01*diag(p) 
#covariance prior for mu_k
Lambda = 1*diag(p)

#generate tree shrinkage according to inverse-gamma distributions
#a1 = 10, a2 =10^4, b = 10, threshold is C = 10
threshold = 10
true_gam_shape_p = rep(0,p)
gam_shape_p = rep(0,p)
gam_rate_p = rep(0,p)
for(a in 1:p){
  num_leaves = length(c(left_leaves[[internal_nodes[a]]], right_leaves[[internal_nodes[a]]]))
  #modelled shrinkage
  if(num_leaves<threshold){ #threshold is C = 10
    gam_shape_p[a] = 10
  } else {
    gam_shape_p[a] = 10^4
  }
  #true shrinkage
  if(num_leaves<threshold){
    true_gam_shape_p[a] = 10
  } else {
    true_gam_shape_p[a] = 10^4
  }
  gam_rate_p[a] = 10
}

#generate true Sigma_k from the above hyperparameters
true_Sigma_ppk = array(0,dim=c(p,p,true_K))
for(k in 1:true_K){
  for(a in 1:p){
    true_Sigma_ppk[a,a,k] = 1/rgamma(1,shape = true_gam_shape_p[a],rate = gam_rate_p[a])
  }
}

#invert Sigma
true_W_ppk =  array(0,dim=c(p,p,true_K))
for(k in 1:true_K){
  true_W_ppk[,,k] = solve(true_Sigma_ppk[,,k])
}

#generate true mu_k from a N(0, Lambda) distribution
true_mu_pk = matrix(0,nrow=p,ncol=true_K)
for(k in 1:true_K){
  true_mu_pk[,k] = Lambda %*% matrix(rnorm(p,0,1),nrow=p)
}
#Set the top level nodes at each topic so that the four subcommunities
#aggregate at different parts of the tree
true_mu_pk[1,1] = 2
true_mu_pk[2,1] = 2
true_mu_pk[3,1] = 2

true_mu_pk[1,2] = 2
true_mu_pk[2,2] = 2
true_mu_pk[3,2] = -2

true_mu_pk[1,3] = -2
true_mu_pk[22,3] = 2
true_mu_pk[23,3] = 2

true_mu_pk[1,4] = -2
true_mu_pk[22,4] = 2
true_mu_pk[23,4] = -2



#generate true phi_d  from a Dir(alpha) dist
true_phi_dk = matrix(0,nrow=D,ncol=true_K)
for (d in 1:D){
  for (k in 1:true_K){
    true_phi_dk[d,] = rgamma(true_K,1,1)
  }
  true_phi_dk[d,] = true_phi_dk[d,]/sum(true_phi_dk[d,])
}

#genreate true psi_dk from a N(mu_k,Sigma_k) dist
true_psi_pdk = array(0,dim=c(p,D,true_K))
for(k in 1:true_K){
  for(d in 1:D){
    true_psi_pdk[,d,k] = chol(true_Sigma_ppk[,,k]) %*% matrix(rnorm(p,0,1),nrow=p) + true_mu_pk[,k]
  }
}

#transform psi values into theta values
true_theta_kda = array(0,dim=c(true_K,D,A))
for (k in 1:true_K){
  for (d in 1:D){
    for (a in 1:p){
      true_theta_kda[k,d,internal_nodes[a]] = exp(true_psi_pdk[a,d,k])/(1+exp(true_psi_pdk[a,d,k]))
    }
  }
}

#find the implicit beta-distributions from the theta
true_beta_kdv = array(0,dim=c(true_K,D,V))
for (d in 1:D){
  for (k in 1:true_K){
    for (leaf in leaves){ #for each leaf
      true_beta_kdv[k,d,leaf] = prod(true_theta_kda[k,d,leaf_success[[leaf]]])*prod((1-true_theta_kda[k,d,leaf_failures[[leaf]]]))
    }
  }
}

#generate true subcommunity assignments
true_ta = NULL
for (d in 1:D){
  true_ta[[d]] = sample(1:true_K,N,replace=TRUE,prob=true_phi_dk[d,])
}

##################
# Generate tokens #
##################
words = matrix(0,nrow=D,ncol=N)
for (d in 1:D){ #for each sample
  for (n in 1:N){ #for each token
    ta = true_ta[[d]][n] #find the true subcommunity assignment
    
    #generate a vector of sample-specific subcommunity probability distributions
    #this is equivalent to using the beta distributions, but finds the
    #implied probability distributions directly from the theta because of legacy code
    prob_vec = rep(0,length(leaves))
    for (leaf in leaves){ 
      prob_vec[leaf] = prod(true_theta_kda[ta,d,leaf_success[[leaf]]])*prod((1-true_theta_kda[ta,d,leaf_failures[[leaf]]]))
    }
    
    #draws the token from the beta distribution
    words[d,n] = sample(leaves,1,prob = prob_vec)
  }
}

#convert the vectors of tokens into a counts matrix
dtm = matrix(0,nrow=D,ncol=V)
for(d in 1:D){
  for(v in 1:V){
    dtm[d,v] = length(which(words[d,]==v))
  }
}


#######################################
# Make long data structures for plots #
#######################################
props = NULL #contains proportions of subcommunities in samples
samps = NULL #contains sample ID
tops = NULL #contains permutation to correct for label switching
model = NULL #contains model ID
Subcommunities = NULL #contains subcommunity ID

#########
# K = 4 #
#########


load("../Results/Simulations/Markov Chains/Sim_C10_TrueK4_K 4 .rda")#load dataset

#unpack list
nc_dnt = results[[5]]
K = dim(nc_dnt)[3]
chain_phi_dki = results[[4]]
iterations = dim(chain_phi_dki)[3]
mu_chain_k_ip = results[[2]]

#calculate posterior mean phi_d
post_phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

#calculate posterior mean mu_k
post_mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  for(a in 1:p){
    post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
  }
}

#find implied theta_ka values
post_theta_ka = matrix(0,nrow=K,ncol=A)
for (k in 1:K){
  for (a in 1:p){
    post_theta_ka[k,internal_nodes[a]] = exp(post_mu_pk[a,k])/(1+exp(post_mu_pk[a,k]))
  }
}

#find implied beta_k values
post_beta_kv = matrix(0,nrow=K,ncol=V)
for (k in 1:K){
  for (leaf in leaves){ #for each leaf
    post_beta_kv[k,leaf] = prod(post_theta_ka[k,leaf_success[[leaf]]])*prod((1-post_theta_ka[k,leaf_failures[[leaf]]]))
  }
}

#do LDA analysis using VEM algorithm - results can change slightly on each run
lda_mod = LDA((dtm),k = K, method = "VEM")
lda_post = posterior(lda_mod)

#correct for label switching between lda and ltn-lda
#minimize SSE between lda beta_k and ltn-lda beta_k
perm = permut(1:K)

#calcu SSE between each beta_k for lda and beta_t for ltn-lda
SSE = matrix(0,nrow=K,ncol=K)
for (k in 1:K){
  for (t in 1:K){
    SSE[k,t] = sum((lda_post$terms[t,]-post_beta_kv[k,])^2)
  }
}

#use the previous results to find total SEE for each permutation
SSE_total = NULL
for (i in 1:nrow(perm)){
  temp = 0
  for (k in 1:K){
    temp = temp + SSE[k,perm[i,k]]
  }
  SSE_total[i] = temp
}
dat_perm = perm[which(SSE_total == min(SSE_total)),]
#this reorders LDA post to match up with LTN post

#permutations to correct for label switching as K changes
top_perm = c(7,8,6,9)

#put data into long-form data structures
for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LTN LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,lda_post$topics[d,dat_perm[k]])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LDA")
    Subcommunities = c(Subcommunities,toString(K))
  }
}

#########
# K = 5 #
#########

load("../Results/Simulations/Markov Chains/Sim_C10_TrueK4_K 5 .rda") #load datat

#unpack list
nc_dnt = results[[5]]
K = dim(nc_dnt)[3]
chain_phi_dki = results[[4]]
mu_chain_k_ip = results[[2]]

#calculate posterior mean phi_d
post_phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

#calculate posterior mean mu_k
post_mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  for(a in 1:p){
    post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
  }
}

#find implied theta_ka values
post_theta_ka = matrix(0,nrow=K,ncol=A)
for (k in 1:K){
  for (a in 1:p){
    post_theta_ka[k,internal_nodes[a]] = exp(post_mu_pk[a,k])/(1+exp(post_mu_pk[a,k]))
  }
}

#find implied beta_k values
post_beta_kv = matrix(0,nrow=K,ncol=V)
for (k in 1:K){
  for (leaf in leaves){ #for each leaf
    post_beta_kv[k,leaf] = prod(post_theta_ka[k,leaf_success[[leaf]]])*prod((1-post_theta_ka[k,leaf_failures[[leaf]]]))
  }
}

#do LDA analysis using VEM algorithm - results can change slightly on each run
lda_mod = LDA((dtm),k = K, method = "VEM")
lda_post = posterior(lda_mod)

#correct for label switching between lda and ltn-lda
#minimize SSE between lda beta_k and ltn-lda beta_k
perm = permut(1:K)

#calcu SSE between each beta_k for lda and beta_t for ltn-lda
SSE = matrix(0,nrow=K,ncol=K)
for (k in 1:K){
  for (t in 1:K){
    SSE[k,t] = sum((lda_post$terms[t,]-post_beta_kv[k,])^2)
  }
}

#use the previous results to find total SEE for each permutation
SSE_total = NULL
for (i in 1:nrow(perm)){
  temp = 0
  for (k in 1:K){
    temp = temp + SSE[k,perm[i,k]]
  }
  SSE_total[i] = temp
}
dat_perm = perm[which(SSE_total == min(SSE_total)),]
#this reorders LDA post to match up with LTN post

#permutations to correct for label switching as K changes
top_perm = c(9,1,7,8,6)

#load data into long-form datastructure
for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LTN LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,lda_post$topics[d,dat_perm[k]])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LDA")
    Subcommunities = c(Subcommunities,toString(K))
  }
}

#########
# K = 7 #
#########

load("../Results/Simulations/Markov Chains/Sim_C10_TrueK4_K 7 .rda") #load dataset

nc_dnt = results[[5]]
K = dim(nc_dnt)[3]
chain_phi_dki = results[[4]]
mu_chain_k_ip = results[[2]]

#calculate posterior mean phi_d
post_phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

#calculate posterior mean mu_k
post_mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  for(a in 1:p){
    post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
  }
}

#find implied theta_ka values
post_theta_ka = matrix(0,nrow=K,ncol=A)
for (k in 1:K){
  for (a in 1:p){
    post_theta_ka[k,internal_nodes[a]] = exp(post_mu_pk[a,k])/(1+exp(post_mu_pk[a,k]))
  }
}

#find implied beta_k values
post_beta_kv = matrix(0,nrow=K,ncol=V)
for (k in 1:K){
  for (leaf in leaves){ #for each leaf
    post_beta_kv[k,leaf] = prod(post_theta_ka[k,leaf_success[[leaf]]])*prod((1-post_theta_ka[k,leaf_failures[[leaf]]]))
  }
}

#do LDA analysis using VEM algorithm - results can change slightly on each run
lda_mod = LDA((dtm),k = K, method = "VEM")
lda_post = posterior(lda_mod)

#correct for label switching between lda and ltn-lda
#minimize SSE between lda beta_k and ltn-lda beta_k
perm = permut(1:K)

#calcu SSE between each beta_k for lda and beta_t for ltn-lda
SSE = matrix(0,nrow=K,ncol=K)
for (k in 1:K){
  for (t in 1:K){
    SSE[k,t] = sum((lda_post$terms[t,]-post_beta_kv[k,])^2)
  }
}

#use the previous results to find total SEE for each permutation
SSE_total = NULL
for (i in 1:nrow(perm)){
  temp = 0
  for (k in 1:K){
    temp = temp + SSE[k,perm[i,k]]
  }
  SSE_total[i] = temp
}
dat_perm = perm[which(SSE_total == min(SSE_total)),]
#this reorders LDA post to match up with LTN post

#permutations to correct for label switching as K changes
top_perm = c(1,8,9,2,4,6,7)

#put data into long-form data structures
for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LTN LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,lda_post$topics[d,dat_perm[k]])
    samps = c(samps,d)
    tops = c(tops,k)
    model = c(model,"LDA")
    Subcommunities = c(Subcommunities,toString(K))
  }
}

##########
# K = 10 #
##########

load("../Results/Simulations/Markov Chains/Sim_C10_TrueK4_K 10 .rda") #load data

#unpack results
nc_dnt = results[[5]]
K = dim(nc_dnt)[3]
chain_phi_dki = results[[4]]
mu_chain_k_ip = results[[2]]

#calculate posterior mean phi_d
post_phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

#calculate posterior mean mu_k
post_mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  for(a in 1:p){
    post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
  }
}

#find implied theta_ka values
post_theta_ka = matrix(0,nrow=K,ncol=A)
for (k in 1:K){
  for (a in 1:p){
    post_theta_ka[k,internal_nodes[a]] = exp(post_mu_pk[a,k])/(1+exp(post_mu_pk[a,k]))
  }
}

#find implied beta_k values
post_beta_kv = matrix(0,nrow=K,ncol=V)
for (k in 1:K){
  for (leaf in leaves){ #for each leaf
    post_beta_kv[k,leaf] = prod(post_theta_ka[k,leaf_success[[leaf]]])*prod((1-post_theta_ka[k,leaf_failures[[leaf]]]))
  }
}

#do LDA analysis using VEM algorithm - results can change slightly on each run
lda_mod = LDA((dtm),k = K, method = "VEM")
lda_post = posterior(lda_mod)

#correct for label switching between lda and ltn-lda
#minimize SSE between lda beta_k and ltn-lda beta_k
perm = permut(1:K)

#calcu SSE between each beta_k for lda and beta_t for ltn-lda
SSE = matrix(0,nrow=K,ncol=K)
for (k in 1:K){
  for (t in 1:K){
    SSE[k,t] = sum((lda_post$terms[t,]-post_beta_kv[k,])^2)
  }
}

#use the previous results to find total SEE for each permutation
SSE_total = NULL
for (i in 1:nrow(perm)){
  temp = 0
  for (k in 1:K){
    temp = temp + SSE[k,perm[i,k]]
  }
  SSE_total[i] = temp
}
dat_perm = perm[which(SSE_total == min(SSE_total)),]
#this reorders LDA post to match up with LTN post

#permutations to correct for label switching as K changes
top_perm = 1:K

#load data into long-form structures
for(k in 1:K){
  for(d in 1:D){
    props = c(props,post_phi_dk[d,k])
    samps = c(samps,d)
    tops = c(tops,top_perm[k])
    model = c(model,"LTN LDA")
    Subcommunities = c(Subcommunities,toString(K))
    
    props = c(props,lda_post$topics[d,dat_perm[k]])
    samps = c(samps,d)
    tops = c(tops,k)
    model = c(model,"LDA")
    Subcommunities = c(Subcommunities,toString(K))
  }
}

###########################
# Make dataframe and plot #
###########################


Sim_abun_df = data.frame(props,samps,tops,model,Subcommunities)

#make new permutation of subcommunities to put most significant ones first
top_perm = c(5,6,7,8,9,1,2,3,4,10)
for(i in 1:nrow(Sim_abun_df)){
  Sim_abun_df$tops[i] = top_perm[Sim_abun_df$tops[i]]
}

save(Sim_abun_df,file="../Results/Simulations/Markov Chains/Sim_abun_df.rda")

##########################
# Add distribution plots #
##########################

load("../Results/Simulations/Markov Chains/Sim_C10_TrueK4_K 10 .rda") #load data

#unpack results
nc_dnt = results[[5]]
K = dim(nc_dnt)[3]
chain_phi_dki = results[[4]]
psi_chain_k_ipd = results[[3]]
mu_chain_k_ip = results[[2]]
Sigma_chain_k_ipp = results[[1]]

#find mean posterior psi_pdk
post_psi_pdk = array(0,dim=c(p,D,K))
for(k in 1:K){
  post_psi_pd = matrix(0,nrow=p,ncol=D)
  for(d in 1:D){
    post_psi_pd[,d] = apply(psi_chain_k_ipd[[k]][1:iterations,,d],2,mean)
  }
  post_psi_pdk[,,k] = post_psi_pd
}

#find implited theta_kda values
post_theta_kda = array(0,dim=c(K,D,A))
for (k in 1:K){
  for (d in 1:D){
    for (a in 1:p){
      post_theta_kda[k,d,internal_nodes[a]] = exp(post_psi_pdk[a,d,k])/(1+exp(post_psi_pdk[a,d,k]))
    }
  }
}

#find the implicit beta-distributions
post_beta_kdv = array(0,dim=c(K,D,V))
for (d in 1:D){
  for (k in 1:K){
    for (leaf in leaves){ #for each leaf
      post_beta_kdv[k,d,leaf] = prod(post_theta_kda[k,d,leaf_success[[leaf]]])*prod((1-post_theta_kda[k,d,leaf_failures[[leaf]]]))
    }
  }
}

#find psoterior mean mu_k
post_mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  for(a in 1:p){
    post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
  }
}

#find theta_ka values
post_theta_ka = matrix(0,nrow=K,ncol=A)
for (k in 1:K){
  for (a in 1:p){
    post_theta_ka[k,internal_nodes[a]] = exp(post_mu_pk[a,k])/(1+exp(post_mu_pk[a,k]))
  }
}

#find the implicit beta-distributions
post_beta_kv = matrix(0,nrow=K,ncol=V)
for (k in 1:K){
  for (leaf in leaves){ #for each leaf
    post_beta_kv[k,leaf] = prod(post_theta_ka[k,leaf_success[[leaf]]])*prod((1-post_theta_ka[k,leaf_failures[[leaf]]]))
  }
}

#find posterior mean sigma_ppk
post_Sigma_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  temp = matrix(0,nrow=p,ncol=p)
  for(i in 1:iterations){
    temp = temp + Sigma_chain_k_ipp[[k]][i,,]
  }
  post_Sigma_ppk[,,k] = temp/iterations
}

#find posterior mean phi_d
post_phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

#do LDA analysis using VEM algorithm - results can change slightly on each run
lda_mod = LDA((dtm),k = K, method = "VEM")
lda_post = posterior(lda_mod)

#correct for label switching between lda and ltn-lda
#minimize SSE between lda beta_k and ltn-lda beta_k
perm = permut(1:K)

#calcu SSE between each beta_k for lda and beta_t for ltn-lda
SSE = matrix(0,nrow=K,ncol=K)
for (k in 1:K){
  for (t in 1:K){
    SSE[k,t] = sum((lda_post$terms[t,]-post_beta_kv[k,])^2)
  }
}

#use the previous results to find total SEE for each permutation
SSE_total = NULL
for (i in 1:nrow(perm)){
  temp = 0
  for (k in 1:K){
    temp = temp + SSE[k,perm[i,k]]
  }
  SSE_total[i] = temp
}
dat_perm = perm[which(SSE_total == min(SSE_total)),]
#this reorders LDA post to match up with LTN post


#######################################
# Plot sample-to-sample heterogeneity #
#######################################

#long-form data structures
vals = NULL #contains values of distributions 
input = NULL #contains the leaf number
samp = NULL #contains the sample number
Distribution = NULL #contains a distributoin id
tops = NULL #which subcommunity, corrects for label switching

for(k in 1:K){
  # ord = order(post_beta_kv[k,])
  
  for(d in 1:D){
    vals = c(vals,post_beta_kdv[k,d,])
    input = c(input,1:V)
    samp = c(samp,rep(d,V))
    Distribution = c( Distribution,rep("LTN LDA Samples",V))
    tops = c(tops,rep(k,V))
  }
  
  vals = c(vals,lda_post$terms[dat_perm[k],])
  input = c(input,1:V)
  samp = c(samp,rep(D + 1,V))
  Distribution = c( Distribution,rep("LDA Average",V))
  tops = c(tops,rep(k,V))
  
  vals = c(vals,post_beta_kv[k,])
  input = c(input,1:V)
  samp = c(samp,rep(D + 2,V))
  Distribution = c( Distribution,rep("LTN LDA Average",V))
  tops = c(tops,rep(k,V))
  
}

#make dataframe
Sim_dist_df = data.frame(input,vals,samp, Distribution,tops)

#reorder for K = 10, same as above
Sim_dist_df$tops_f = factor(Sim_dist_df$tops, levels=c('7','8','6','9','10','4','5','1','2','3'))
top_perm = c(5,6,7,8,9,1,2,3,4,10)
for(i in 1:nrow(df)){
  Sim_dist_df$tops[i] = top_perm[Sim_dist_df$tops[i]]
}

save(Sim_dist_df,file="../Results/Simulations/Markov Chains/Sim_dist_df.rda")
