library(Rcpp)
library(RcppArmadillo)
library(pgdraw)
library(RcppDist) #probably unnecessary, but part of legacy code
library(topicmodels)
library(combinat)

taskID = 1 #used for setting random seed and saving results

set.seed(taskID)

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


##################
# MAke test set #
#################
#####################
# Generate test set #
#####################

#this set is used to estimate the sample-specific parameters for the test set

#set some test-set specific parameters - number of sampels and reads/sampel
D_t = D
N_t = N/2 #the test sample only corresponds to one half of of the held out sample
#this lets us do document-completion

#generate phi_d on test set from Dir(alpha) dist
test_phi_d = matrix(0,nrow=D_t,ncol=true_K)
for (d in 1:D_t){
  test_phi_d[d,] = rgamma(true_K,shape=alpha,rate=1)
  test_phi_d[d,] = test_phi_d[d,]/sum(test_phi_d[d,])
}

#draw true subcommunity assignments for test set
test_true_z_dn = matrix(0,nrow=D_t,ncol=N_t)
for (d in 1:D_t){
  test_true_z_dn[d,] = sample(1:true_K,N_t, prob = test_phi_d[d,],replace=TRUE)
}

#genreate true psi_dk for test set from N(mu_k,Sigma_k) using true mu_k and Sigma_k
test_psi_pdk = array(0,dim=c(p,D_t,true_K))
for(k in 1:true_K){
  for(d in 1:D_t){
    test_psi_pdk[,d,k] = chol(true_Sigma_ppk[,,k]) %*% matrix(rnorm(p,0,1),nrow=p) + true_mu_pk[,k]
  }
}

#Find implied theta values
test_theta_kda = array(0,dim=c(true_K,D_t,A))
for (k in 1:true_K){
  for (d in 1:D_t){
    for (a in 1:p){
      test_theta_kda[k,d,internal_nodes[a]] = exp(test_psi_pdk[a,d,k])/(1+exp(test_psi_pdk[a,d,k]))
    }
  }
}

#find the implicit beta-distributions
test_beta_kdv = array(0,dim=c(true_K,D_t,V))
for (d in 1:D_t){
  for (k in 1:true_K){
    for (leaf in leaves){ #for each leaf
      test_beta_kdv[k,d,leaf] = prod(test_theta_kda[k,d,leaf_success[[leaf]]])*prod((1-test_theta_kda[k,d,leaf_failures[[leaf]]]))
    }
  }
}

#Using sample-specific beta and subcommunity assignments, generate vectors of tokens
test_word_dn = matrix(0,nrow=D_t,ncol=N_t)
for(d in 1:D_t){
  for (n in 1:ncol(test_true_z_dn)){
    test_word_dn[d,n] = sample(1:V,1,prob=test_beta_kdv[test_true_z_dn[d,n],d,])
  }
}

#transform vector of tokens into various data structures
test_docs = NULL
test_docs_C = NULL
test_docs_mat_C = matrix(0,nrow=D_t,ncol=N_t)
for (d in 1:D_t){
  test_docs[[d]] = test_word_dn[d,]
  test_docs_C[[d]] = test_docs[[d]]-1
  test_docs_mat_C[d,] = test_docs_C[[d]]
}


#####################
# Generate New Data #
#####################

#this extra set is used to evalute the perplexity of the held out set

#find true subcommunity assignments from the test-set specific phi_d
extra_true_z_dn = matrix(0,nrow=D_t,ncol=N_t)
for (d in 1:D_t){
  extra_true_z_dn[d,] = sample(1:true_K,N_t, prob = test_phi_d[d,],replace=TRUE)
}

#generate vectors of tokens using the beta and subcommunity assignments
extra_word_dn = matrix(0,nrow=D_t,ncol=N_t)
for(d in 1:D_t){
  for (n in 1:ncol(extra_word_dn)){
    extra_word_dn[d,n] = sample(1:V,1,prob=test_beta_kdv[extra_true_z_dn[d,n],d,])
  }
}

#convert form of token vector
extra_docs = NULL
for (d in 1:D_t){
  extra_docs[[d]] = extra_word_dn[d,]
}

#convert token vectors into counts matrix
extra_dtm = matrix(0,nrow=D_t,ncol=V)
for(d in 1:D_t){
  for(v in 1:V){
    extra_dtm[d,v] = length(which(extra_docs[[d]]==v))
  }
}

###################
# Initialization! #
###################

#fake K
K= true_K #set number of topics the model expects to find - here set to truth

#take list of token vectors for the training set
docs = NULL
docs_C = NULL
docs_mat_C = matrix(0,nrow=D,ncol=N)
for (d in 1:nrow(words)){
  docs[[d]] = words[d,]
  docs_C[[d]] = docs[[d]]-1
  docs_mat_C[d,] = docs_C[[d]]
}

# Initialize subcommunity assignment vectors according to discrete-uniform
# find the token-subcommunity matrix wt
#wt[k,v] is the number of sequencing reads v assigned to subcommunity k
wt = matrix(0, K, V) # initialize word-topic count matrix
ta = lapply(docs, function(x) rep(0, length(x))) # initialize topic assignment list
ta_C = lapply(docs, function(x) rep(0, length(x))) 
ta_mat_C_dn = matrix(0,nrow=D,ncol=N)
for(d in 1:length(docs)){ # for each sample
  for(w in 1:length(docs[[d]])){ # for each token in sample d
    ta[[d]][w] = sample(1:K, 1) # randomly assign subcommunity to token w.
    ta_C[[d]][w] = ta[[d]][w]-1
    ta_mat_C_dn[d,w] = ta_C[[d]][w]
    
    ti = ta[[d]][w] # find the subcommunity assignment of token w in sample d
    wi = docs[[d]][w] # find the token w in sample d
    wt[ti,wi] = wt[ti,wi]+1 # update token-subcommunity count matrix by 1 in appropriate place
  }
}

#find the sample-subcommunity count matrix dt
#dt[d,k] is the number of tokens assigned to subcommunity k in sample d
dt = matrix(0, length(docs), K)
for(d in 1:length(docs)){ # for each sample d
  for(t in 1:K){ # for each subcommunity t
    dt[d,t] = sum(ta[[d]]==t) # count tokens in sample d assigned to subcommunity t   
  }
}

#find counts of reads by sample and subcommunity
#wc_dwt[d,w,k] contains the number of times read w has subcommunity k in sample d
wc_dwt = array(0,dim=c(D,V,K)) 
for(d in 1:length(docs)){ # for each sample
  for (k in 1:K){
    ta.temp = docs[[d]][which(ta[[d]]==k)]
    for(w in 1:length(ta.temp)){ # for each token
      ti = k #find the assigned topic
      wi = ta.temp[w] #find the relavent read
      wc_dwt[d,wi,ti] = wc_dwt[d,wi,ti]+1 # updatecount matrix     
    }
  }
}

#find counts of read assigned to subcommunity descended from specific nodes
#nc_dnt[d,a,k] is the number of read descened from node a in sample d assigned to subcommunity k
nc_dnt =  array(0,dim=c(D,max(tree.edge),K)) 
nc_dnt_prop = array(0,dim=c(D,max(tree.edge),K)) #for memory purposes in legacy code
for(d in 1:length(docs)){ # for each sample
  for (k in 1:K){
    nc_dnt[d,1:max(leaves),k]=wc_dwt[d,,k] #count on the leaves is equal to wc_dwt
    nc_dnt_prop[d,1:max(leaves),k]=wc_dwt[d,,k]#count on the leaves is equal to wc_dwt 
    #work way up tree iteratively, summing counts descended from each node
    for(i in length(layers):1){
      for (j in 1:length(layers[[i]])){
        if (length(descendants[[layers[[i]][j]]])>0){
          nc_dnt[d,layers[[i]][j],k] = sum(nc_dnt[d,descendants[[layers[[i]][j]]],k])
          nc_dnt_prop[d,layers[[i]][j],k] = sum(nc_dnt_prop[d,descendants[[layers[[i]][j]]],k])
        }
      }
    }
  }
}

#Now we initialize parameter values

#Initialize kappa according to nc_dnt
kappa_pdk = array(0,dim=c(p,D,K))
for(k in 1:K){
  for(a in 1:p){
    for(d in 1:D){
      kappa_pdk[a,d,k] = nc_dnt[d,descendants[[internal_nodes[a]]][1],k] - nc_dnt[d,internal_nodes[a],k]/2
    }
  }
}

#initialize phi from a Dir(1) distribution
phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    phi_dk[d,] = rgamma(K,1,1)
  }
  phi_dk[d,] = phi_dk[d,]/sum(phi_dk[d,])
}

#initialize psi from a N(0,I) dist
psi_pdk = array(0,dim=c(p,D,K))
for(k in 1:K){
  for(d in 1:D){
    psi_pdk[,d,k] = chol(diag(p)) %*% matrix(rnorm(p,0,1),nrow=p)
  }
}


#find implicit theta values
theta_kda = array(0,dim=c(K,D,A))
for (k in 1:K){
  for (d in 1:D){
    for (a in 1:p){
      theta_kda[k,d,internal_nodes[a]] = exp(psi_pdk[a,d,k])/(1+exp(psi_pdk[a,d,k]))
    }
  }
}

#find the implicit beta-distributions
beta_kdv = array(0,dim=c(K,D,V))
for (d in 1:D){
  for (k in 1:K){
    for (leaf in leaves){ #for each leaf
      beta_kdv[k,d,leaf] = prod(theta_kda[k,d,leaf_success[[leaf]]])*prod((1-theta_kda[k,d,leaf_failures[[leaf]]]))
    }
  }
}

#initialize v from PG distributions according to nc_dnt
v_pdk = array(0,dim=c(p,D,K))
for(k in 1:K){
  for(d in 1:D){
    for(a in 1:p){
      if(nc_dnt[d,internal_nodes[a],k]<1){
        v_pdk[a,d,k] = 0
      } else {
        v_pdk[a,d,k] = pgdraw(nc_dnt[d,internal_nodes[a],k],psi_pdk[a,d,k])
      }
    }
  }
}

#initialize Sigma according to modelled hyperparameter a1, a2, b
Sigma_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  for(a in 1:p){
    Sigma_ppk[a,a,k] = 1/rgamma(1,shape = gam_shape_p[a],rate = gam_rate_p[a])
  }
}
#find the inverses
W_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  W_ppk[,,k] = solve(Sigma_ppk[,,k])
}

#initialize mu according to N(0,I) dist
mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  mu_pk[,k] = diag(p) %*% matrix(rnorm(p,0,1),nrow=p)
}

#compute Lamba Inverse - this was harder before Lambda was diagonal
Lambda_inv = solve(Lambda)


#Initiealize chains and set chain lengths
warmup = 10000 #number of warmup iterations
iterations = 1000 #number of iterations recorded after warmup
thin = 10 #amount by which we thin recorded iteratoins.  Run thin*iterations after warmup

#Pre-allocate chains
chain_phi_dki = array(0,dim=c(D,K,iterations))
psi_chain_k_ipd = NULL
mu_chain_k_ip = NULL
Sigma_chain_k_ipp = NULL
for(k in 1:K){
  psi_chain_k_ipd[[k]] = array(0,dim=c(iterations,p,D))
  mu_chain_k_ip[[k]] = matrix(0,nrow=iterations,ncol=p) 
  Sigma_chain_k_ipp[[k]] = array(0,dim = c(iterations,p,p))
}
p_z = rep(0,K)

#pre-allocate output
results = NULL
results[[5]] = nc_dnt
results[[4]] = chain_phi_dki
results[[3]] =   psi_chain_k_ipd
results[[2]] = mu_chain_k_ip
results[[1]] = Sigma_chain_k_ipp

#source C function
sourceCpp("PG_C_Sampler_Strict_Sparse_Tree_Shrink.cpp")

#this is here to run the pgdraw function in the C code --- probably a better way of doing this
f_pg = function(b,c){
  return(pgdraw(b,c))
}

#run the Gibbs sampler!
begin = Sys.time()
results = LTN(results, f_pg, Sigma_ppk, W_ppk, mu_pk, v_pdk, psi_pdk, kappa_pdk, theta_kda, beta_kdv, Lambda_inv, Phi, gam_shape_p, gam_rate_p, chain_phi_dki, psi_chain_k_ipd, mu_chain_k_ip, Sigma_chain_k_ipp, nc_dnt, dt, descendants_mat_C, ta_C, docs_C, ancestors_C, internal_nodes_C, leaf_success_C, leaf_failures_C, K, p, D, V, alpha, iterations, warmup, thin)
end = Sys.time()

#record output as different data strucures
nc_dnt = results[[5]]
chain_phi_dki = results[[4]]
psi_chain_k_ipd = results[[3]]
mu_chain_k_ip = results[[2]]
Sigma_chain_k_ipp = results[[1]]

#find posterior mean estimates for the psi
post_psi_pdk = array(0,dim=c(p,D,K))
for(k in 1:K){
  post_psi_pd = matrix(0,nrow=p,ncol=D)
  for(d in 1:D){
    post_psi_pd[,d] = apply(psi_chain_k_ipd[[k]][1:iterations,,d],2,mean)
  }
  post_psi_pdk[,,k] = post_psi_pd
}

#find posteiror mean estimates for the mu
post_mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  for(a in 1:p){
    post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
  }
}

#find posterior mean estimates for the Sigma
post_Sigma_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  temp = matrix(0,nrow=p,ncol=p)
  for(i in 1:iterations){
    temp = temp + Sigma_chain_k_ipp[[k]][i,,]
  }
  post_Sigma_ppk[,,k] = temp/iterations
}

#find posterior mean estimates for the phi
post_phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

####################################################
# Generate and do Document Completion on test set  #
####################################################

#transform vector of tokens into various data structures
test_docs = NULL
test_docs_C = NULL
test_docs_mat_C = matrix(0,nrow=D_t,ncol=N_t)
for (d in 1:D_t){
  test_docs[[d]] = test_word_dn[d,]
  test_docs_C[[d]] = test_docs[[d]]-1
  test_docs_mat_C[d,] = test_docs_C[[d]]
}


# Initialize subcommunity assignment vectors according to discrete-uniform
# find the token-subcommunity matrix wt
#wt[k,v] is the number of sequencing reads v assigned to subcommunity k
wt = matrix(0, K, V) # initialize count matrix
ta = lapply(test_docs, function(x) rep(0, length(x))) # initialize sub assignment list
ta_C = lapply(test_docs, function(x) rep(0, length(x))) 
ta_mat_C_dn = matrix(0,nrow=D_t,ncol=N_t)
for(d in 1:length(test_docs)){ # for each sample
  for(w in 1:length(test_docs[[d]])){ # for each token in sample d
    ta[[d]][w] = sample(1:K, 1) # randomly assign subcommunity to token w.
    ta_C[[d]][w] = ta[[d]][w]-1
    ta_mat_C_dn[d,w] = ta_C[[d]][w]
    
    ti = ta[[d]][w] # find the subcommunity assignment of token w in sample d
    wi = test_docs[[d]][w] # find the token w in sample d
    wt[ti,wi] = wt[ti,wi]+1 # update token-subcommunity count matrix by 1 in appropriate place
  }
}

#find the sample-subcommunity count matrix dt
#dt[d,k] is the number of tokens assigned to subcommunity k in sample d
dt = matrix(0, length(test_docs), K)
for(d in 1:length(test_docs)){ # for each sample d
  for(t in 1:K){ # for each subcommunity t
    dt[d,t] <- sum(ta[[d]]==t) # count tokens in sample d assigned to subcommunity t   
  }
}

#find counts of reads by sample and subcommunity
#wc_dwt[d,w,k] contains the number of times read w has subcommunity k in sample d
wc_dwt = array(0,dim=c(D_t,V,K)) 
for(d in 1:length(test_docs)){ # for each sample
  for (k in 1:K){
    ta.temp = test_docs[[d]][which(ta[[d]]==k)]
    for(w in 1:length(ta.temp)){ # for each token
      ti = k # find the assigned topic
      wi = ta.temp[w] #find the relavent read
      wc_dwt[d,wi,ti] = wc_dwt[d,wi,ti]+1 # update count matrix      
    }
  }
}

#find counts of read assigned to subcommunity descended from specific nodes
#nc_dnt[d,a,k] is the number of read descened from node a in sample d assigned to subcommunity k
nc_dnt =  array(0,dim=c(D_t,max(tree.edge),K)) 
nc_dnt_prop = array(0,dim=c(D_t,max(tree.edge),K)) #for memory purposes in legacy code
for(d in 1:length(test_docs)){ # for each sample
  for (k in 1:K){
    nc_dnt[d,1:max(leaves),k]=wc_dwt[d,,k] #count on the leaves is equal to wc_dwt
    nc_dnt_prop[d,1:max(leaves),k]=wc_dwt[d,,k] #count on the leaves is equal to wc_dwt
    #work way up tree iteratively, summing counts descended from each node
    for(i in length(layers):1){
      for (j in 1:length(layers[[i]])){
        if (length(descendants[[layers[[i]][j]]])>0){
          nc_dnt[d,layers[[i]][j],k] = sum(nc_dnt[d,descendants[[layers[[i]][j]]],k])
          nc_dnt_prop[d,layers[[i]][j],k] = sum(nc_dnt_prop[d,descendants[[layers[[i]][j]]],k])
        }
      }
    }
  }
}

#Now we initialize parameter values

#Initialize kappa according to nc_dnt
kappa_pdk = array(0,dim=c(p,D_t,K))
for(k in 1:K){
  for(a in 1:p){
    for(d in 1:D_t){
      kappa_pdk[a,d,k] = nc_dnt[d,descendants[[internal_nodes[a]]][1],k] - nc_dnt[d,internal_nodes[a],k]/2
    }
  }
}

#initialize phi from a Dir(1) distribution
phi_dk = matrix(0,nrow=D_t,ncol=K)
for(d in 1:D_t){
  for(k in 1:K){
    phi_dk[d,] = rgamma(K,1,1)
  }
  phi_dk[d,] = phi_dk[d,]/sum(phi_dk[d,])
}

#initialize psi from a N(0,I) dist
psi_pdk = array(0,dim=c(p,D_t,K))
for(k in 1:K){
  for(d in 1:D_t){
    psi_pdk[,d,k] = chol(diag(p)) %*% matrix(rnorm(p,0,1),nrow=p)
  }
}


#find implied theta values
theta_kda = array(0,dim=c(K,D_t,A))
for (k in 1:K){
  for (d in 1:D_t){
    for (a in 1:p){
      theta_kda[k,d,internal_nodes[a]] = exp(psi_pdk[a,d,k])/(1+exp(psi_pdk[a,d,k]))
    }
  }
}

#find implied beta values
beta_kdv = array(0,dim=c(K,D_t,V))
for (d in 1:D_t){
  for (k in 1:K){
    for (leaf in leaves){ #for each leaf
      beta_kdv[k,d,leaf] = prod(theta_kda[k,d,leaf_success[[leaf]]])*prod((1-theta_kda[k,d,leaf_failures[[leaf]]]))
    }
  }
}

#initialize v from PG distributions according to nc_dnt
v_pdk = array(0,dim=c(p,D_t,K))
for(k in 1:K){
  for(d in 1:D_t){
    for(a in 1:p){
      if(nc_dnt[d,internal_nodes[a],k]<1){
        v_pdk[a,d,k] = 0
      } else {
        v_pdk[a,d,k] = pgdraw(nc_dnt[d,internal_nodes[a],k],psi_pdk[a,d,k])
      }
    }
  }
}

#initialize Sigma equal to posterior means from training set.  
# These do not update in the perplexity gibbs sampler
Sigma_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  Sigma_ppk[,,k] = post_Sigma_ppk[,,k]
}
#find the inverses
W_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  W_ppk[,,k] = solve(Sigma_ppk[,,k])
}

#initialize mu according to posterior means from training set
#these do not update in the gibbs sampler
mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  mu_pk[,k] = post_mu_pk[,k]
}

#compute Lamba Inverse
Lambda_inv = solve(Lambda)


#Initialize chains and set chain lengths
warmup = 10000 #number of warmup iterations
iterations = 1000 #number of iterations recorded after warmup
thin = 10 #amount by which we thin recorded iteratoins.  Run thin*iterations after warmup

#pre-allocate chains
chain_phi_dki = array(0,dim=c(D_t,K,iterations))
psi_chain_k_ipd = NULL
mu_chain_k_ip = NULL
Sigma_chain_k_ipp = NULL
for(k in 1:K){
  psi_chain_k_ipd[[k]] = array(0,dim=c(iterations,p,D_t))
  mu_chain_k_ip[[k]] = matrix(0,nrow=iterations,ncol=p) 
  Sigma_chain_k_ipp[[k]] = array(0,dim = c(iterations,p,p))
}
p_z = rep(0,K)

#pre-allocate results
results = NULL
results[[5]] = nc_dnt
results[[4]] = chain_phi_dki
results[[3]] =   psi_chain_k_ipd
results[[2]] = mu_chain_k_ip
results[[1]] = Sigma_chain_k_ipp

#source C function
sourceCpp("PG_C_Sampler_Strict_Sparse_Tree_Shrink_Perp.cpp")

#this is here to run the pgdraw function in the C code --- probably a better way of doing this
f_pg = function(b,c){
  return(pgdraw(b,c))
}

#run the modified Gibbs sampler to estimate sample-speficic parameters on the test set
begin = Sys.time()
results = LTN(results, f_pg, Sigma_ppk, W_ppk, mu_pk, v_pdk, psi_pdk, kappa_pdk, theta_kda, beta_kdv, Lambda_inv, Phi, gam_shape_p, gam_rate_p, chain_phi_dki, psi_chain_k_ipd, mu_chain_k_ip, Sigma_chain_k_ipp, nc_dnt, dt, descendants_mat_C, ta_C, test_docs_C, ancestors_C, internal_nodes_C, leaf_success_C, leaf_failures_C, K, p, D_t, V, alpha, iterations, warmup, thin)
end = Sys.time()

#pull result out of list
nc_dnt = results[[5]]
chain_phi_dki = results[[4]]
psi_chain_k_ipd = results[[3]]
mu_chain_k_ip = results[[2]]
Sigma_chain_k_ipp = results[[1]]

#find posterior mean estimates for the phi on the test set
test_post_phi_dk = matrix(0,nrow=D_t,ncol=K)
for(d in 1:D_t){
  for(k in 1:K){
    test_post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

#find the implied theta  and beta values at each iteration in the test set
post_theta_kdai = array(0,dim=c(K,D_t,A,iterations))
post_beta_kdvi = array(0,dim=c(K,D_t,V,iterations))
for(i in 1:iterations){
  #convert psi to theta
  for (k in 1:K){
    for (d in 1:D_t){
      for (a in 1:p){
        post_theta_kdai[k,d,internal_nodes[a],i] = exp(psi_chain_k_ipd[[k]][i,a,d])/(1+exp(psi_chain_k_ipd[[k]][i,a,d]))
      }
    }
  }
  #convert theta to beta
  for (d in 1:D_t){
    for (k in 1:K){
      for (leaf in leaves){ #for each leaf
        post_beta_kdvi[k,d,leaf,i] = prod(post_theta_kdai[k,d,leaf_success[[leaf]],i])*prod((1-post_theta_kdai[k,d,leaf_failures[[leaf]],i]))
      }
    }
  }
}


#######################
# Evaluate Perplexity #
#######################

#estimate perplexity on the extra set
doc_comp_prob = rep(0,D_t)
for (d in 1:D_t){
  total = 0
  for (i in 1:iterations){
    subtotal =0 
    for (w in 1:length(extra_docs[[d]])){
      wid = extra_docs[[d]][w]
      
      subtotal = subtotal + log(sum(chain_phi_dki[d,,i]*post_beta_kdvi[,d,wid,i]))
      
    }
    total = subtotal+total
  }
  doc_comp_prob[d] = total/(iterations)
}
den = 0 
for (d in 1:length(extra_docs)){
  den = den + length(extra_docs[[d]])
}
ltn_perp = exp(-sum(doc_comp_prob)/den)

##########
# topicmodels
# use the topic models package to estimate perplexity for LDA
# same number of iterations as LTN-LDA
# use a Gibbs sampler and a variational EM algorithm

LDA_out = LDA(dtm,k = K, method = "Gibbs", control = list(iter = 20000))
gibbs_perp = perplexity(LDA_out, newdata = extra_dtm)

LDA_out = LDA(dtm,k = K, method = "VEM")
vem_perp = perplexity(LDA_out, newdata = extra_dtm)

print(ltn_perp)
print(gibbs_perp)
print(vem_perp)

perp = c(ltn_perp,gibbs_perp,vem_perp)

#save the perplexity results, name according to taskID
save(perp,file = paste("../Results/Simulations/Perplexity/perp",toString(taskID),".rda",sep=""))

# ##############
# # optional code for L_p distances
# #set p
# p_norm = 2
# 
# #find L_[ distance for phi for various permutations of 1:K to correct for label switching
# perms = permn(1:K)
# L_p_vec = rep(0,length(perms))
# for(i in 1:length(perms)){
#   
#   L_p = 0
#   for(k in 1:K){
#     L_p = L_p + sum((true_phi_dk[,k] - post_phi_dk[,perms[[i]][k]])^p_norm)
#   }
#   L_p_vec[i] = L_p^(1/p_norm)
# }
# #find the permutaiton which minimizes label switching
# perm = perms[[which(L_p_vec == min(L_p_vec))]]
# 
# #find L_p distances for the three paraemters
# L_phi = sum((true_phi_dk - post_phi_dk[,perm])^p_norm)^(1/p_norm)
# L_beta_kdv = sum((true_beta_kdv - post_beta_kdv[perm,,])^p_norm)^(1/p_norm)
# L_beta_kv = sum((true_beta_kv - post_beta_kv[perm,])^p_norm)^(1/p_norm)
# 
# #record and save
# L = c(L_phi,L_beta_kdv,L_beta_kv)
# save(L,file = paste("L",toString(taskID),".rda",sep=""))