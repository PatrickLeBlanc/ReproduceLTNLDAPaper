time1 = Sys.time()

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

set.seed(1)

#Define Parameters
true_K = 2
D = 50
N = 1000
alpha = 1


#let's hard code the tree into an edge matrix
edge_list= matrix(0,nrow = 14, ncol = 2)
edge_list[1,] = c(9,10)
edge_list[2,] = c(9,1)



edge_list[3,] = c(11,12)
edge_list[4,] = c(12,2)
edge_list[5,] = c(12,3)
edge_list[6,] = c(11,13)
edge_list[7,] = c(13,4)
edge_list[8,] = c(13,5)

edge_list[9,] = c(10,14)
edge_list[10,] = c(14,15)
edge_list[11,] = c(15,7)
edge_list[12,] = c(15,8)
edge_list[13,] = c(14,6)
edge_list[14,] = c(10,11)


toy_tree = tree.edge = edge_list

root = setdiff(toy_tree[,1],toy_tree[,2])
internal_nodes = unique(toy_tree[,1])
internal_nodes = sort(internal_nodes)
internal_nodes_C = internal_nodes - 1
A = max(internal_nodes)
leaves = setdiff(toy_tree[,2],toy_tree[,1])
V = length(leaves)

#Generate some tree data structures

#Find the descendants
descendants = NULL
descendants_mat = matrix(0,ncol=2,nrow=max(tree.edge))
for(i in 1:max(tree.edge)){
  if (sum(which(tree.edge[,1]==i))>0){
    descendants[[i]] = tree.edge[which(tree.edge[,1] == i),2]
    descendants_mat[i,] = descendants[[i]]
  }
}
descendants_mat_C = descendants_mat - 1

#descendants[[i]] contains the two immediate descendants of node i

#Find the parents
parents = NULL
for (i in 1:max(tree.edge)){
  if (sum(which(tree.edge[,2]==i))>0){
    parents[i] = tree.edge[which(tree.edge[,2]==i),1]
  }
}
#parents[i] contains the parent of node i

#Find the ancestors
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
#ancestors[[i]] contains all of the ancestors of node i

#Cut the tree into layers
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
#layers is set up such that the i^th entry i the list contains the nodes in 
#the i^th layer of the list

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
#left_leaves[[node]] contains the leaves left-descended from node
#right_leaves[[node]] contains the leaves right-descended from node


#need to find, for each leaf
# the nodes which have to succeed
# the nodes which have to fail
#in order for the leaf to be selected

#successes first 
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
#leaf_success[[leaf]] contains the nodes which have to be successful for leaf to appear

#and now failures
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



#come up with a mapping to 
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

Phi = 0.01*diag(p)
Lambda = 1*diag(p)

#generate tree shrinkage
true_C = 4
true_gam_shape_p = rep(0,p)
gam_rate_p = rep(0,p)
for(a in 1:p){
  num_leaves = length(c(left_leaves[[internal_nodes[a]]], right_leaves[[internal_nodes[a]]]))
  
  #true shrinkage
  if(num_leaves<true_C){
    true_gam_shape_p[a] = 10
  } else {
    true_gam_shape_p[a] = 10^4
  }
  gam_rate_p[a] = 10
}

#generate true Sigma_k
true_Sigma_ppk = array(0,dim=c(p,p,true_K))
for(k in 1:true_K){
  for(a in 1:p){
    true_Sigma_ppk[a,a,k] = 1/rgamma(1,shape = true_gam_shape_p[a],rate = gam_rate_p[a])
  }
}

true_W_ppk =  array(0,dim=c(p,p,true_K))
for(k in 1:true_K){
  true_W_ppk[,,k] = solve(true_Sigma_ppk[,,k])
}

#generate true mu_k
true_mu_pk = matrix(0,nrow=p,ncol=true_K)
for(k in 1:true_K){
  true_mu_pk[,k] = Lambda %*% matrix(rnorm(p,0,1),nrow=p)
}
# true_mu_pk[,1] = 2
true_mu_pk[1,1] = 2
true_mu_pk[2,1] = 2
true_mu_pk[3,1] = 0

true_mu_pk[1,2] = 2
true_mu_pk[2,2] = -2
true_mu_pk[3,2] = 0


#generate true phi_d  
true_phi_dk = matrix(0,nrow=D,ncol=true_K)
for (d in 1:D){
  for (k in 1:true_K){
    true_phi_dk[d,] = rgamma(true_K,1,1)
  }
  true_phi_dk[d,] = true_phi_dk[d,]/sum(true_phi_dk[d,])
}

#genreate true psi_dk
true_psi_pdk = array(0,dim=c(p,D,true_K))
for(k in 1:true_K){
  for(d in 1:D){
    true_psi_pdk[,d,k] = chol(true_Sigma_ppk[,,k]) %*% matrix(rnorm(p,0,1),nrow=p) + true_mu_pk[,k]
  }
}
# for(d in 1:D){
#   for(k in 1:true_K){
#     true_psi_pdk[,d,k] = true_mu_pk[,k]
#   }
# }

#draw theta_kdA values
true_theta_kda = array(0,dim=c(true_K,D,A))
for (k in 1:true_K){
  for (d in 1:D){
    for (a in 1:p){
      true_theta_kda[k,d,internal_nodes[a]] = exp(true_psi_pdk[a,d,k])/(1+exp(true_psi_pdk[a,d,k]))
    }
  }
}

#find the implicit beta-distributions
true_beta_kdv = array(0,dim=c(true_K,D,V))
for (d in 1:D){
  for (k in 1:true_K){
    for (leaf in leaves){ #for each leaf
      true_beta_kdv[k,d,leaf] = prod(true_theta_kda[k,d,leaf_success[[leaf]]])*prod((1-true_theta_kda[k,d,leaf_failures[[leaf]]]))
    }
  }
}

# par(mfrow=c(1,true_K))
# for(k in 1:true_K){
#   matplot(true_beta_kdv[k,,],type="l")
# }
# par(mfrow=c(1,1))

#find true topic assignments
true_ta = NULL
for (d in 1:D){
  true_ta[[d]] = sample(1:true_K,N,replace=TRUE,prob=true_phi_dk[d,])
}

##################
# Generate words #
##################
words = matrix(0,nrow=D,ncol=N)
for (d in 1:D){ #for each document
  for (n in 1:N){ #for each word
    ta = true_ta[[d]][n] #find the topic assignment
    
    prob_vec = rep(0,length(leaves))
    for (leaf in leaves){ #for each leaf
      prob_vec[leaf] = prod(true_theta_kda[ta,d,leaf_success[[leaf]]])*prod((1-true_theta_kda[ta,d,leaf_failures[[leaf]]]))
    }
    
    words[d,n] = sample(leaves,1,prob = prob_vec)
  }
}

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

D_t = D
N_t = N/2
test_phi_d = matrix(0,nrow=D_t,ncol=true_K)
for (d in 1:D_t){
  test_phi_d[d,] = rgamma(true_K,shape=alpha,rate=1)
  test_phi_d[d,] = test_phi_d[d,]/sum(test_phi_d[d,])
}

test_true_z_dn = matrix(0,nrow=D_t,ncol=N_t)
for (d in 1:D_t){
  test_true_z_dn[d,] = sample(1:true_K,N_t, prob = test_phi_d[d,],replace=TRUE)
}

#genreate true psi_dk
test_psi_pdk = array(0,dim=c(p,D_t,true_K))
for(k in 1:true_K){
  for(d in 1:D_t){
    test_psi_pdk[,d,k] = chol(true_Sigma_ppk[,,k]) %*% matrix(rnorm(p,0,1),nrow=p) + true_mu_pk[,k]
  }
}

#draw theta_kdA values
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

test_word_dn = matrix(0,nrow=D_t,ncol=N_t)
for(d in 1:D_t){
  for (n in 1:ncol(test_true_z_dn)){
    test_word_dn[d,n] = sample(1:V,1,prob=test_beta_kdv[test_true_z_dn[d,n],d,])
  }
}

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

extra_true_z_dn = matrix(0,nrow=D_t,ncol=N_t)
for (d in 1:D_t){
  extra_true_z_dn[d,] = sample(1:true_K,N_t, prob = test_phi_d[d,],replace=TRUE)
}

extra_word_dn = matrix(0,nrow=D_t,ncol=N_t)
for(d in 1:D_t){
  for (n in 1:ncol(extra_word_dn)){
    extra_word_dn[d,n] = sample(1:V,1,prob=test_beta_kdv[extra_true_z_dn[d,n],d,])
  }
}

extra_docs = NULL
for (d in 1:D_t){
  extra_docs[[d]] = extra_word_dn[d,]
}

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
K= true_K#set number of topics the model expects to find
#transform into docs
docs = NULL
docs_C = NULL
docs_mat_C = matrix(0,nrow=D,ncol=N)
for (d in 1:nrow(words)){
  docs[[d]] = words[d,]
  docs_C[[d]] = docs[[d]]-1
  docs_mat_C[d,] = docs_C[[d]]
}

#some code below is stolen from 

## 1. Initialize topic assignments to words in each doc.  
#2. Generate word-topic count matrix.
wt <- matrix(0, K, V) # initialize word-topic count matrix
ta <- lapply(docs, function(x) rep(0, length(x))) # initialize topic assignment list
ta_C = lapply(docs, function(x) rep(0, length(x))) 
ta_mat_C_dn = matrix(0,nrow=D,ncol=N)
for(d in 1:length(docs)){ # for each document
  for(w in 1:length(docs[[d]])){ # for each token in document d
    ta[[d]][w] <- sample(1:K, 1) # randomly assign topic to token w.
    
    #####################
    # SET TO TRUTH HERE #
    #####################
    # ta[[d]][w] = true_ta[[d]][w]
    ta_C[[d]][w] = ta[[d]][w]-1
    ta_mat_C_dn[d,w] = ta_C[[d]][w]
    ti <- ta[[d]][w] # topic index
    wi <- docs[[d]][w] # wordID for token w
    wt[ti,wi] <- wt[ti,wi]+1 # update word-topic count matrix     
  }
}


dt <- matrix(0, length(docs), K)
for(d in 1:length(docs)){ # for each document d
  for(t in 1:K){ # for each topic t
    dt[d,t] <- sum(ta[[d]]==t) # count tokens in document d assigned to topic t   
  }
}

wc_dwt = array(0,dim=c(D,V,K)) #word-count-by-document list
for(d in 1:length(docs)){ # for each document
  for (k in 1:K){
    ta.temp = docs[[d]][which(ta[[d]]==k)]
    for(w in 1:length(ta.temp)){ # for each word
      ti = k # topic index
      wi = ta.temp[w] # wordID for token w
      wc_dwt[d,wi,ti] = wc_dwt[d,wi,ti]+1 # update word-topic count matrix     
    }
  }
}

#wcbdbt lists word-counts-by-document-by-topic
#wcbdbt[d,, contains a word-topic matrix for document d
#wcbdbdt[d,w,k] contains the number of times word w has topic k in doc d

#Need to transform into a list of node counts-by-topic-by-doc
nc_dnt =  array(0,dim=c(D,max(tree.edge),K)) #node-count-by-document-by-topic
nc_dnt_prop = array(0,dim=c(D,max(tree.edge),K))
for(d in 1:length(docs)){ # for each document
  for (k in 1:K){
    nc_dnt[d,1:max(leaves),k]=wc_dwt[d,,k] #word-count,doc d,top k on the leaves
    nc_dnt_prop[d,1:max(leaves),k]=wc_dwt[d,,k] #word-count,doc d,top k on the leaves
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
#this is set up to be doc - node - topic
#might want to redo at some point? to be more readable

#Initialize kappa
kappa_pdk = array(0,dim=c(p,D,K))
for(k in 1:K){
  for(a in 1:p){
    for(d in 1:D){
      kappa_pdk[a,d,k] = nc_dnt[d,descendants[[internal_nodes[a]]][1],k] - nc_dnt[d,internal_nodes[a],k]/2
    }
  }
}

#initialize phi
phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    phi_dk[d,] = rgamma(K,1,1)
  }
  phi_dk[d,] = phi_dk[d,]/sum(phi_dk[d,])
  # phi_dk[d,] = true_phi_dk[d,]
}

#initialize psi
psi_pdk = array(0,dim=c(p,D,K))
for(k in 1:K){
  for(d in 1:D){
    psi_pdk[,d,k] = chol(diag(p)) %*% matrix(rnorm(p,0,1),nrow=p)
    # psi_pdk[,d,k] = true_psi_pdk[,d,k]
  }
}


#draw theta_kdA values
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

#initialize v
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


#Set moelled C
C = true_C
#initialize modelled shrinkage
gam_shape_p = rep(0,p)
for(a in 1:p){
  num_leaves = length(c(left_leaves[[internal_nodes[a]]], right_leaves[[internal_nodes[a]]]))
  #modelled shrinkage
  if(num_leaves<C){
    gam_shape_p[a] = 10
  } else {
    gam_shape_p[a] = 10^4
  }
  
}

#initialize Sigma
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

#initialize mu
mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  mu_pk[,k] = diag(p) %*% matrix(rnorm(p,0,1),nrow=p)
  # mu_pk[,k] = true_mu_pk[,k]
}

#compute Lamba Inverse
Lambda_inv = solve(Lambda)


#Initiealize chains

warmup = 1000
iterations = 1000
thin = 1
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

results = NULL
results[[5]] = nc_dnt
results[[4]] = chain_phi_dki
results[[3]] =   psi_chain_k_ipd
results[[2]] = mu_chain_k_ip
results[[1]] = Sigma_chain_k_ipp

sourceCpp("PG_C_Sampler_Strict_Sparse_Tree_Shrink.cpp")

f_pg = function(b,c){
  return(pgdraw(b,c))
}

#do LTN-LDA analysis
begin = Sys.time()
results = LTN(results, f_pg, Sigma_ppk, W_ppk, mu_pk, v_pdk, psi_pdk, kappa_pdk, theta_kda, beta_kdv, Lambda_inv, Phi, gam_shape_p, gam_rate_p, chain_phi_dki, psi_chain_k_ipd, mu_chain_k_ip, Sigma_chain_k_ipp, nc_dnt, dt, descendants_mat_C, ta_C, docs_C, ancestors_C, internal_nodes_C, leaf_success_C, leaf_failures_C, K, p, D, V, alpha, iterations, warmup, thin)
end = Sys.time()

#unpack the results
nc_dnt = results[[5]]
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
    # post_psi_pd[,d] = mean(psi_chain_k_ipd[[k]][1:iterations,,d])
  }
  post_psi_pdk[,,k] = post_psi_pd
}

#find theta_kdA values
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


post_mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  for(a in 1:p){
    post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
  }
}

#draw theta_kdA values
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

#psoterior Sigma means
post_Sigma_ppk = array(0,dim=c(p,p,K))
for(k in 1:K){
  temp = matrix(0,nrow=p,ncol=p)
  for(i in 1:iterations){
    temp = temp + Sigma_chain_k_ipp[[k]][i,,]
  }
  post_Sigma_ppk[,,k] = temp/iterations
}

#posterior phi means
post_phi_dk = matrix(0,nrow=D,ncol=K)
for(d in 1:D){
  for(k in 1:K){
    post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
  }
}

model = list(Mean_Post_Phi_d = post_phi_dk,
             Mean_Post_Beta_kd = post_beta_kdv,
             Mean_Post_Beta_k = post_beta_kv,
             Chain_Phi = chain_phi_dki,
             Chain_Psi = psi_chain_k_ipd,
             Chain_Mu = mu_chain_k_ip,
             Chain_Sigma = Sigma_chain_k_ipp
)

plot(true_phi_dk,post_phi_dk,
     xlim = c(0,1),ylim = c(0,1),
     xlab = "True Phi", ylab = "Mean Posterior Phi"
)
