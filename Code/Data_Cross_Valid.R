#some of these are unnecssary
library(phyloseq)
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
library(randtests)

#set taskID and random seed
#assume taskID is in {1,2,...,196} because of computing restraints
#originally, values of taskID were run in parallel on Duke's Computing Cluster
taskID = 1
set.seed(1)

##########
# set up code to partition data into 4 roughly equal sections
size_vec = c(13,13,14,14)

#find a hold out set
samps = 1:54
hold_out = NULL
hold_out[[length(size_vec+1)]] = rep(0,13)
for(i in 1:length(size_vec)){
  hold_out[[i]] = sample(samps,size_vec[i],replace = FALSE)
  for(j in 1:length(hold_out[[i]])){
    samps = samps[-which(samps == hold_out[[i]][j])]
  }
}
hold_out[[length(size_vec)+1]] = NULL


##########
# Setup indexing for converting taskID to K,C,partition
M = 196 #max possible taskID
ind_mat = matrix(0,nrow=M,ncol = 3)

#Can change values of K and C
for(k in 1:7){ #let K vary in {2,3,...,8}
  for(c in 1:7){ #let C vary in {1,2,...,7}
    for(i in 1:4){ #let i vary in and control which set is held out
      ind = 4*7*(k-1) + 4*(c-1) + i
      ind_mat[ind,1] = k + 1
      ind_mat[ind,2] = c
      ind_mat[ind,3] = i
    }
  }
}

# #Once C has been chosen, run K over a large range.  
# #The upper bound is the size of the training set
# 
# M = 168
# ind_mat = matrix(0,nrow=M,ncol = 3)
# c = 8
# for(k in 1:42){
#   for(i in 1:4){
#     ind = 4*(k-1) + i
#     ind_mat[ind,1] = k + 1
#     ind_mat[ind,2] = c
#     ind_mat[ind,3] = i
#   }
# }

set.seed(taskID)

#load dataset from Dethlefsen and Relman 2011
#current version taken from Kris Sankaran's github
load('../Data/abt.rda')
#rename
ps = abt
print(ps)

#subset to patient F
ps = subset_samples(ps,ind == "F")
print(ps)


#retrieve taxonomic table and merge taxa at finest known level
tax = tax_table(ps)
for(i in 1:nrow(tax)){
  name = tax[i,1]
  for(j in 2:ncol(tax)){
    if(tax[i,j] %in% c("","uncultured","Incertae Sedis")){
      tax[i,j] = name
    } else {
      name = tax[i,j]
    }
  }
}
colnames(tax)
tax_table(ps) = tax



#merge sequencing reads together at finest known level
ps_merge  = tax_glom(ps,taxrank = "Taxon_8",NArm = TRUE)
print(ps_merge)

#Exclude taxa whic do not occur at least 100 times across all 54 samples
ps_filter = filter_taxa(ps_merge,function(x) sum(x) > 100, TRUE)
print(ps_filter)

#retrive the counts matrix
dtm = otu_table(ps_filter)

#####
# retreive parameters and set hyperparameters
D = ncol(dtm)
K = ind_mat[taskID,1] #set this according to the taskID
alpha = 1

#retrive the phylogenetic tree. 
tree = phy_tree(ps_filter)

#rename the data structure because of legacy code
toy_tree = tree.edge = tree$edge

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


#################################
# Split data into training/test #
#################################


##############
# test documents
# use the taskID to select which of the four partitions is held out
D_t = length(hold_out[[ind_mat[taskID,3]]])
#the test set is used to find local parameters values on the held-out set
test_dtm = dtm[,hold_out[[ind_mat[taskID,3]]]]
#the extra set is used to evaluate perplexity given these estimates
extra_dtm = matrix(0,nrow=V,ncol=D_t)

#partition these sets roughly in half
#some reads which occur small numbers of times are not split
#these will not affect the perplexity estimate since they occur less than 3 times
for(d in 1:D_t){
  for(v in 1:V){
    total = test_dtm[v,d]
    if( total > 3){
      extra_dtm[v,d] = sample( (round(total/2) - 1):(round(total/2) + 1),1)
      test_dtm[v,d] = test_dtm[v,d]-extra_dtm[v,d]
    }
    
  }
}

#unfold the test count matrix into token vectors
test_word_dn = NULL
test_word_dn[[D_t+1]] = rep(0,10000)
for (d in 1:D_t){ #for each document
  word = NULL
  
  for(v in 1:V){
    word = c(word,rep(v,test_dtm[v,d]))
  }
  
  word = word[sample(1:length(word),replace = FALSE)]
  test_word_dn[[d]] = word
}
test_word_dn[[D_t+1]] = NULL

#rename token vectors for legacy code
test_docs = NULL
test_docs_C = NULL
for (d in 1:D_t){
  test_docs[[d]] = test_word_dn[[d]]
  test_docs_C[[d]] = test_docs[[d]]-1
}

###############
# extra documents

#unfold the extra count matrix into token vectors
extra_word_dn = NULL
extra_word_dn[[D_t+1]] = rep(0,10000)
for (d in 1:D_t){ #for each document
  word = NULL
  
  for(v in 1:V){
    word = c(word,rep(v,extra_dtm[v,d]))
  }
  
  word = word[sample(1:length(word),replace = FALSE)]
  extra_word_dn[[d]] = word
}
extra_word_dn[[D_t+1]] = NULL

#rename the token vectors for legacy code
extra_docs = NULL
for (d in 1:D_t){
  extra_docs[[d]] = extra_word_dn[[d]]
}

################
# reestablish counts matrix without the held-out set
dtm = dtm[,-hold_out[[ind_mat[taskID,3]]]]
#rename parameters
D = ncol(dtm)

###############################
# Convert dtm to words matrix #
###############################
#unfold counts matrix to token vectors
words = NULL
words[[D+1]] = rep(0,10000)
for (d in 1:D){ #for each document
  word = NULL
  
  for(v in 1:V){
    word = c(word,rep(v,dtm[v,d]))
  }
  
  word = word[sample(1:length(word),replace = FALSE)]
  words[[d]] = word
}
words[[D+1]] = NULL


###################
# Initiliazation! #
###################
#rename token vectors for legacy code purpsoes
docs = NULL
docs_C = NULL
for (d in 1:D){
  docs[[d]] = words[[d]]
  docs_C[[d]] = docs[[d]]-1
}

# Initialize subcommunity assignment vectors according to discrete-uniform
# find the token-subcommunity matrix wt
#wt[k,v] is the number of sequencing reads v assigned to subcommunity k
wt = matrix(0, K, V) # initialize token-sub count matrix
ta = lapply(docs, function(x) rep(0, length(x))) # initialize sub assignment list
ta_C = lapply(docs, function(x) rep(0, length(x))) 
for(d in 1:length(docs)){ # for each sample
  for(w in 1:length(docs[[d]])){ # for each token in sample d
    ta[[d]][w] = sample(1:K, 1) # randomly assign sub to token w.
    ta_C[[d]][w] = ta[[d]][w]-1
    
    ti = ta[[d]][w] # find the subcommunity assignment of token w in sample d
    wi = docs[[d]][w] # find the token w in sample d
    wt[ti,wi] = wt[ti,wi]+1 # update token-subcommunity count matrix by 1 in appropriate place  
  }
}

#find the sample-subcommunity count matrix dt
#dt[d,k] is the number of tokens assigned to subcommunity k in sample d
dt = matrix(0, length(docs), K)
for(d in 1:length(docs)){ # for each sample d
  for(t in 1:K){ # for each sub t
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
      ti = k # find the assigned subcommunity
      wi = ta.temp[w] # find the relavent read
      wc_dwt[d,wi,ti] = wc_dwt[d,wi,ti]+1 # update count matrix     
    }
  }
}

#find counts of read assigned to subcommunity descended from specific nodes
#nc_dnt[d,a,k] is the number of read descened from node a in sample d assigned to subcommunity k
nc_dnt =  array(0,dim=c(D,max(tree.edge),K)) 
nc_dnt_prop = array(0,dim=c(D,max(tree.edge),K))  #for memory purposes in legacy code
for(d in 1:length(docs)){ # for each sample
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

#initiliaze ig parameters
# C is chosen by taskID
# a1 = b = 10, a2 = 10^4
gam_shape_p = rep(0,p)
gam_rate_p = rep(0,p)
for(a in 1:p){
  num_leaves = length(c(left_leaves[[internal_nodes[a]]], right_leaves[[internal_nodes[a]]]))
  #modelled shrinkage
  if(num_leaves<ind_mat[taskID,2]){
    gam_shape_p[a] = 10
  } else {
    gam_shape_p[a] = 10^4
  }
  gam_rate_p[a] = 10
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

#Initailize hyperparameter covraince matrix
Lambda = 1*diag(p)
#Phi is a holdover from legacy code, is as input from C function, not used
Phi = diag(p)

#compute Lamba Inverse
Lambda_inv = solve(Lambda)

#Initiealize chains
warmup = 10000 #number of warmup iterations 
iterations = 1000 #number of iterations recorded after warmup
thin = 10 #amount by which we thin recorded iteratoins.  Run thin*iterations after warmup

#pre-allocate chains
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

#source C function containing Gibbs Sampler
sourceCpp("PG_C_Sampler_Strict_Sparse_Tree_Shrink.cpp")

#this is here to run the pgdraw function in the C code --- probably a better way of doing this
f_pg = function(b,c){
  return(pgdraw(b,c))
}

#run the Gibbs Sampler!
begin = Sys.time()
results = LTN(results, f_pg, Sigma_ppk, W_ppk, mu_pk, v_pdk, psi_pdk, kappa_pdk, theta_kda, beta_kdv, Lambda_inv, Phi, gam_shape_p, gam_rate_p, chain_phi_dki, psi_chain_k_ipd, mu_chain_k_ip, Sigma_chain_k_ipp, nc_dnt, dt, descendants_mat_C, ta_C, docs_C, ancestors_C, internal_nodes_C, leaf_success_C, leaf_failures_C, K, p, D, V, alpha, iterations, warmup, thin)
end = Sys.time()

#record output as different data strucures
nc_dnt = results[[5]]
chain_phi_dki = results[[4]]
psi_chain_k_ipd = results[[3]]
mu_chain_k_ip = results[[2]]
Sigma_chain_k_ipp = results[[1]]

#find posterior mean psi
post_psi_pdk = array(0,dim=c(p,D,K))
for(k in 1:K){
  post_psi_pd = matrix(0,nrow=p,ncol=D)
  for(d in 1:D){
    post_psi_pd[,d] = apply(psi_chain_k_ipd[[k]][1:iterations,,d],2,mean)
  }
  post_psi_pdk[,,k] = post_psi_pd
}

#find posterior mean mu
post_mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  for(a in 1:p){
    post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
  }
}

#find posterior mean Sigma
post_Sigma_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  temp = matrix(0,nrow=p,ncol=p)
  for(i in 1:iterations){
    temp = temp + Sigma_chain_k_ipp[[k]][i,,]
  }
  post_Sigma_ppk[,,k] = temp/iterations
}

#find posterior mean phi
post_phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

##########################################
# Fit modified Gibbs Sampler on test set #
##########################################

# Initialize subcommunity assignment vectors according to discrete-uniform
# find the token-subcommunity matrix wt
#wt[k,v] is the number of sequencing reads v assigned to subcommunity k
wt = matrix(0, K, V) # initialize count matrix
ta = lapply(test_docs, function(x) rep(0, length(x))) # initialize sub assignment list
ta_C = lapply(test_docs, function(x) rep(0, length(x))) 
for(d in 1:length(test_docs)){ # for each sample
  for(w in 1:length(test_docs[[d]])){ # for each token in sample d
    ta[[d]][w] = sample(1:K, 1) # randomly assign sub to token w.
    ta_C[[d]][w] = ta[[d]][w]-1
    
    ti = ta[[d]][w] # find the subcommunity assignment of token w in sample d
    wi = test_docs[[d]][w] # find the token w in sample d
    wt[ti,wi] <- wt[ti,wi]+1 # update token-subcommunity count matrix by 1 in appropriate place  
  }
}

#find the sample-subcommunity count matrix dt
#dt[d,k] is the number of tokens assigned to subcommunity k in sample d
dt = matrix(0, length(test_docs), K)
for(d in 1:length(test_docs)){ # for each sample d
  for(t in 1:K){ # for each sub t
    dt[d,t] <- sum(ta[[d]]==t) # count tokens in sample d assigned to subcommunity t   
  }
}

#find counts of reads by sample and subcommunity
#wc_dwt[d,w,k] contains the number of times read w has subcommunity k in sample d
wc_dwt = array(0,dim=c(D_t,V,K)) 
for(d in 1:length(test_docs)){ # for each sample
  for (k in 1:K){
    ta.temp = test_docs[[d]][which(ta[[d]]==k)]
    for(w in 1:length(ta.temp)){ # for each word
      ti = k # find the assigned topic
      wi = ta.temp[w] #find the relavent read
      wc_dwt[d,wi,ti] = wc_dwt[d,wi,ti]+1 # update count matrix      
    }
  }
}


#find counts of read assigned to subcommunity descended from specific nodes
#nc_dnt[d,a,k] is the number of read descened from node a in sample d assigned to subcommunity k
nc_dnt =  array(0,dim=c(D_t,max(tree.edge),K)) 
nc_dnt_prop = array(0,dim=c(D_t,max(tree.edge),K)) #needed for memory purposes in legacy code
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

#initialize phi according to a Dir(1) dist
phi_dk = matrix(0,nrow=D_t,ncol=K)
for(d in 1:D_t){
  for(k in 1:K){
    phi_dk[d,] = rgamma(K,1,1)
  }
  phi_dk[d,] = phi_dk[d,]/sum(phi_dk[d,])
}

#initialize psi according to a N(0,I) distribution
psi_pdk = array(0,dim=c(p,D_t,K))
for(k in 1:K){
  for(d in 1:D_t){
    psi_pdk[,d,k] = chol(diag(p)) %*% matrix(rnorm(p,0,1),nrow=p)
  }
}


#find theta_kdA values
theta_kda = array(0,dim=c(K,D_t,A))
for (k in 1:K){
  for (d in 1:D_t){
    for (a in 1:p){
      theta_kda[k,d,internal_nodes[a]] = exp(psi_pdk[a,d,k])/(1+exp(psi_pdk[a,d,k]))
    }
  }
}

#find the implicit beta-distributions
beta_kdv = array(0,dim=c(K,D_t,V))
for (d in 1:D_t){
  for (k in 1:K){
    for (leaf in leaves){ #for each leaf
      beta_kdv[k,d,leaf] = prod(theta_kda[k,d,leaf_success[[leaf]]])*prod((1-theta_kda[k,d,leaf_failures[[leaf]]]))
    }
  }
}

#initialize v from a PG dist according to nc_dnt and kappa
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

#initialize Sigma to posterior mean esitmates - not updated in modified gibbs sampler
Sigma_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  Sigma_ppk[,,k] = post_Sigma_ppk[,,k]
}
#find the inverses
W_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  W_ppk[,,k] = solve(Sigma_ppk[,,k])
}

#initialize mu according to posterior mean estimates - not updated in modified gibbs sampler
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

#pre-allocae results
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

#unpack the results
nc_dnt = results[[5]]
chain_phi_dki = results[[4]]
psi_chain_k_ipd = results[[3]]
mu_chain_k_ip = results[[2]]
Sigma_chain_k_ipp = results[[1]]

#find posterior mean phi estimates for the test set
test_post_phi_dk = matrix(0,nrow=D_t,ncol=K)
for(d in 1:D_t){
  for(k in 1:K){
    test_post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

#find theta and beta estimates for the test set at each iteration
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
LDA_out = LDA(t(dtm),k = K, method = "Gibbs", control = list(iter = 20000))
gibbs_perp = perplexity(LDA_out, newdata = t(extra_dtm))

LDA_out = LDA(t(dtm),k = K, method = "VEM")
vem_perp = perplexity(LDA_out, newdata = t(extra_dtm))

#save perplexity results
perp = c(ltn_perp,gibbs_perp,vem_perp)

save(perp,file = paste("../Results/Data/Perplexity/Data_perp_",toString(taskID),".rda",sep=""))
