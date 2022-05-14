library(phyloseq)
library(LTNLDA)
library(topicmodels)

#used for setting random seed and saving results
#set using slurm ID
taskID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

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

#merge sequencing reads together at taxon_4
ps  = tax_glom(ps,taxrank = "Taxon_4",NArm = TRUE)
print(ps)

#take tree from phyloseq object and use to construct simulated dataset
tree = phy_tree(ps)
tree.edge = tree$edge
#rename the data structure because of legacy code
toy_tree = tree.edge

#set random seed
set.seed(2)

#################################
# Simulate otu table from model #
#################################


#Define Parameters
true_K = 4
D = 50
N = 10000
alpha = 1
C = 9

####################################################
# Convert Tree Edge Matrix into useable structures #
####################################################

#find characteristics of tree
#find the root node
root = setdiff(tree.edge[,1],tree.edge[,2]) 

#find the internal nodes and reorder
internal_nodes = unique(tree.edge[,1])
internal_nodes = sort(internal_nodes)
internal_nodes_C = internal_nodes - 1

#find the maximum of the internal nodes
A = max(internal_nodes)

#find the set of leaves
leaves = setdiff(tree.edge[,2],tree.edge[,1])

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


#find the number of leaves descended from each node
#num_leaves[a] is the number of leaves descended from node a
num_leaves = rep(0,p)
for(x in 1:p){
  num_leaves[x] = length(c(left_leaves[[internal_nodes[x]]], right_leaves[[internal_nodes[x]]]))
}

#make data structure recording which nodes belong in the upper part of the tree and which belong in the lower part
#U_nodes list the nodes in the upper matrix; p_U is the number of such nodes
#L_nodes list the nodes in the upper matrix; p_L is the number of such nodes
U_nodes = which(num_leaves >= C)
U_nodes_C = U_nodes - 1
p_U = length(U_nodes)
L_nodes = which(num_leaves < C)
L_nodes_C = L_nodes - 1
p_L = length(L_nodes)

#define scale matrices for Inverse-Wishart covariance priors for upper and lower matrices
#initialized to diagonal matrices
Phi_U = diag(p_U)
Phi_L = diag(p_L)
#also define Lambda
Lambda = diag(p)

#set hyperparameters for covariance priors
a_L = 5
b_L = 5*(p_L+2)
a_U = 10^4
b_U = 10

######################
# Generate true data #
######################


###################
#                 #
###################
# Data Generation #
###################
#                 #
###################

#generate true Sigma_k^U
true_Sigma_U_ppk = array(0,dim=c(p_U,p_U,true_K))
for(k in 1:true_K){
  for(a in 1:p_U){
    true_Sigma_U_ppk[a,a,k] = 1/rgamma(1,shape = a_U,rate = b_U)
  }
}
true_W_U_ppk =  array(0,dim=c(p_U,p_U,true_K))
for(k in 1:true_K){
  true_W_U_ppk[,,k] = solve(true_Sigma_U_ppk[,,k])
}

#generate true Sigma_k^L

#G-matrix
g_prior = 0.25
true_G_L_ppk = array(0,dim=c(p_L,p_L,true_K))
#make data structures to store coordinates of existing/nonexisting edges
num_upptri_L = p_L*(p_L-1)/2
upper_coord_L = matrix(0,nrow = num_upptri_L,ncol=2)
upper_coord_L_C = matrix(0,nrow = num_upptri_L,ncol=2)
for(k in 1:true_K){
  coord_ct = 1
  for(i in 1:(p_L-1)){
    for(j in (i+1):p_L){
      #initialize G
      true_G_L_ppk[i,j,k] = sample(0:1,size = 1,prob = c(1-g_prior,g_prior))
      
      #Create  a numbering on the edges
      upper_coord_L[coord_ct,1] = i
      upper_coord_L[coord_ct,2] = j
      
      #iterate counter
      coord_ct = coord_ct + 1
    }
  }
  diag(true_G_L_ppk[,,k]) = 1
}
upper_coord_L_C = upper_coord_L - 1

true_W_L_ppk =  array(0,dim=c(p_L,p_L,true_K))
for(k in 1:true_K){
  true_W_L_ppk[,,k] = f_gwish(true_G_L_ppk[,,k],a_L*(p_L+2),b_L*Phi_L)
}
true_Sigma_L_ppk = array(0,dim=c(p_L,p_L,true_K))
for(k in 1:true_K){
  true_Sigma_L_ppk[,,k] = solve(true_W_L_ppk[,,k])
}


#combine Sigma_k^L and Sigma_k_^U into Sigma_k
true_Sigma_ppk = array(0,dim = c(p,p,true_K))
for(k in 1:true_K){
  true_Sigma_ppk[U_nodes,U_nodes,k] = true_Sigma_U_ppk[,,k]
  true_Sigma_ppk[L_nodes,L_nodes,k] = true_Sigma_L_ppk[,,k]
}
true_W_ppk = array(0,dim=c(p,p,true_K))
for(k in 1:true_K){
  true_W_ppk[,,k] = solve(true_Sigma_ppk[,,k])
}

#generate true mu_k
true_mu_pk = matrix(0,nrow=p,ncol=true_K)
for(k in 1:true_K){
  true_mu_pk[,k] = Lambda %*% matrix(rnorm(p,0,1),nrow=p)
}
#establish four subcommunities
true_mu_pk[1,1] = 3
true_mu_pk[2,1] = 2
true_mu_pk[3,1] = 2
true_mu_pk[4,1] = 2

true_mu_pk[1,2] = 3
true_mu_pk[2,2] = 2
true_mu_pk[3,2] = 2
true_mu_pk[4,2] = -2

true_mu_pk[1,3] = 3
true_mu_pk[2,3] = -2

true_mu_pk[1,4] = -3


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

par(mfrow=c(1,true_K))
for(k in 1:true_K){
  matplot(true_beta_kdv[k,,],type="l")
}
par(mfrow=c(1,1))

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
dtm = t(dtm) #transpose to fit requirements

rownames(dtm) = tree$tip.label

#construct phyloseq object
OTU = otu_table(dtm,taxa_are_rows = TRUE)
ps_sim = phyloseq(OTU,tree)

##################
# Make test set #
#################
#####################
# Generate test set #
#####################

#this set is used to estimate the sample-specific parameters for the test set
#this is partitioned into a test set and an extra set because of legacy code

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

#convert token vectors into counts matrix
test_dtm = matrix(0,nrow=D_t,ncol=V)
for(d in 1:D_t){
  for(v in 1:V){
    test_dtm[d,v] = length(which(test_docs[[d]]==v))
  }
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

#combine these two sets into one 
dtm = test_dtm + extra_dtm
dtm = t(dtm)#transpose to fit requirements
rownames(dtm) = tree$tip.label

#construct phyloseq object
OTU = otu_table(dtm,taxa_are_rows = TRUE)
ps_sim_test = phyloseq(OTU,tree)


#################################
# Run model and find perplexity #
#################################

#set the number of subcommunities
K = 2
#set C to truth
C = 9

#run the LTNLDA model - default values correspond to true values except for K and C
model = LTNLDA_cov(ps_sim,K = K, C = C)

#find the perplexity on a test set
perp = LTNLDA_Perplexity(model, test_ps, iterations, burnin, thin)

ltn_perp = perp$Perplexity

###############################
# run lda and find perpelxity #
###############################

##########
# topicmodels
# use the topic models package to estimate perplexity for LDA
# same number of iterations as LTN-LDA
# use a Gibbs sampler and a variational EM algorithm

dtm = t(otu_table(ps_sim))
LDA_out = LDA(dtm,k = K, method = "Gibbs", control = list(iter = 20000))
dtm = t(otu_table(ps_sim_test))
gibbs_perp = perplexity(LDA_out, newdata = dtm)

dtm = t(otu_table(ps_sim))
LDA_out = LDA(dtm,k = K, method = "VEM")
dtm = t(otu_table(ps_sim_test))
vem_perp = perplexity(LDA_out, newdata = extra_dtm)

perp = c(ltn_perp,gibbs_perp,vem_perp)

#save the perplexity results, name according to taskID
save(perp,file = paste("../Results/Simulations/Perplexity/perp",toString(taskID),".rda",sep=""))
