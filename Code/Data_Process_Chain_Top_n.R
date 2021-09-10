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
#here taskID controls the number of modelled subcommunities K
taskID = 2
set.seed(1)

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

#retrieve the taxonomix table
tax = tax_table(ps_filter)

#####
# retreive parameters and set hyperparameters
D = ncol(dtm)
K = taskID #set this according to the taskID
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

#############
# load data #
#############
load("../Results/Data/Markov Chains/Sank_C_8_K_ 7 _data_results.rda")

#unpack results
nc_dnt = results[[5]]
K = dim(nc_dnt)[3]
chain_phi_dki = results[[4]]
iterations = dim(chain_phi_dki)[3]
psi_chain_k_ipd = results[[3]]
mu_chain_k_ip = results[[2]]
Sigma_chain_k_ipp = results[[1]]

#calculate posterior mean psi
post_psi_pdk = array(0,dim=c(p,D,K))
for(k in 1:K){
  post_psi_pd = matrix(0,nrow=p,ncol=D)
  for(d in 1:D){
    post_psi_pd[,d] = apply(psi_chain_k_ipd[[k]][1:iterations,,d],2,mean)
  }
  post_psi_pdk[,,k] = post_psi_pd
}

#find implied theta values
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

#find posterior mean amu
post_mu_pk = matrix(0,nrow=p,ncol=K)
for(k in 1:K){
  for(a in 1:p){
    post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
  }
}

#find implied theta values
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

#do LDA analysis using VEM algorithm - results can change slightly on each run
lda_mod = LDA(t(dtm),k = K, method = "VEM")
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

##################################################
# Make a plot to show top_n for both LTN and LDA #
##################################################

#find posterior mean values for average ASV-sub dist
terms_post_ltn = post_beta_kv
terms_post_lda = lda_post$terms[dat_perm,]

#find ordering on these dists in decreasing proportions
order_mat_ltn = matrix(0,nrow = K, ncol = V)
order_mat_lda = matrix(0,nrow = K, ncol = V)
for(k in 1:K){
  order_mat_ltn[k,] = order(terms_post_ltn[k,], decreasing = TRUE)
  order_mat_lda[k,] = order(terms_post_lda[k,], decreasing = TRUE)
}

#make long form datastructures
names = NULL #contains names of ASVs
props = NULL #contains proportions of ASVs
topic = NULL #contains subcommunity IDs
model = NULL #contains model ID

#find the top n names for each model
top_n = 5
for(k in 1:K){
  names = c(names,tax[order_mat_ltn[k,top_n:1],dim(tax)[2]])
  props = c(props,terms_post_ltn[k,order_mat_ltn[k,top_n:1]])
  topic = c(topic,rep(k,top_n))
  model = c(model,rep("LTN LDA",top_n))
  
  names = c(names,tax[order_mat_lda[k,top_n:1],dim(tax)[2]])
  props = c(props,terms_post_lda[k,order_mat_lda[k,top_n:1]])
  topic = c(topic,rep(k,top_n))
  model = c(model,rep("LDA",top_n))
}

#make a dataframe with all of the data
Data_topn_df = data.frame(names,props,topic,model)

#reorder topics to match other plots
top_perm = c(1,7,6,4,2,5,3)
for(i in 1:nrow(Data_topn_df)){
  Data_topn_df$topic[i] = top_perm[Data_topn_df$topic[i]]
}

save(Data_topn_df,file="../Results/Data/Markov Chains/Data_topn_df.rda")
