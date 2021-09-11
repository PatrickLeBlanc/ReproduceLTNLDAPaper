# ReproduceLTNLDAPaper

This repository contains the code to reproduce the numerical examples and data analysis in LeBlanc and Ma, as well as some results.  The code used to produce the results in the paper was written before the LTNLDA package and so does not make use of the package.  In addition, it contains a number of places where the code could be improved to make it more efficient; in the interests of remaining close to the original code, I have not modified these.  

Included in this repository are three main directories: Code, Data, and Results.  The Code repository contains the code necessary to reproduce the results in the paper; the Data repository contains a tree structure and a phyloseq object, both retrieved from (Dethlefsen and Relman 2011) which were used in the paper; Results contains summaries of results for both numerical simulations and data analysis.

Much of this code was originally designed to run on the Duke Computing Cluster (DCC).  Multiple instances of the same code would be run in parallel on different cores, with the only difference being the taskID.  Thus, in much of the code the taskID parameter controls many other parameters.  

We now present a description of each file and directory in the repository:

# Code

Contains R and Rcpp scripts used in the analysis.  Some of these have been modified to run in general R environments instead of on the Duke Computing Cluster.

## Data_Chain.R

Data_Chain.R contains the code used to load the dataset from (Dethlefsen and Relman 2011), subset the dataset, and run an LTN-LDA Gibbs sampler on the entire dataset.  The parameter $K$ is controlled by the taskID parameter, and $C$ is set equal to 8 as suggested by perplexity analysis.  The output is a list where each entry contains Markov Chains for different parameters.  All parameters are included except for the Polya-Gamma auxiliary variables.  Note that these Markov chains are not included in the repository because they exceed github's maximum file size.  

## Data_Cross_Valid.R

Data_Cross_Valid.R contains the code to perform cross-validation on the dataset from (Dethlefsen and Relman 2011).  The dataset is partitioned into four approximately equal datasets.  According to the taskID parameter, we select the values for $K$ and $C$ as well as select which of these partitions to treat as the hold-out set.  It fits a Gibbs sampler for LTN-LDA on the training set, and uses the output of this to estimate the perplexity on the test set using a modified Gibbs sampler.  The output is a vector containing three three perplexity estimates: the LTN-LDA Gibbs sampler, the LDA Gibbs sampler, and the LDA VEM algorithm.

This code also contains an option, obtained via commentating out and uncommenting relevant parts of the code, to obtain perplexity estimates for a fixed value of $C$ and varying levels of $K$ controlled by the taskID.  

## Data_Perplexity_Plots.R

Data_Perplexity_Plots.R loads the perplexity results obtained from running many different iterations of the Data_Cross_Valid.R code.  It then produces a plot of perplexity curves for LTN-LDA for varying levels of $K$ as we vary $C$ as well as a plot of perplexity curves for LTN-LDA and LDA as $K$ varies while $C = 8$.  

## Data_Plot_Chains.R

Data_Plot_Chains.R loads the three dataframes produced by Data_Process_Chain_Dist.R and Data_Process_Chain_Top_n.R.  It uses the dataframes produced by Data_Process_Chain_Dist.R to produce a side-by-side plot where the left plot displays posterior mean subcommunity abundance estimates for LTN-LDA and LDA as $K$ varies in $\{3,4,7\}$ and the right plot displays subcommunity compositions for LTN-LDA and LDA.  The dataframe produced by Data_Process_Chain_Top_n.R is used to produce a plot displaying the top $n$ ASVs per subcommunity for LDA and LTN-LDA as well as their prevalence.  

## Data_Process_Chain_Dist.R

Data_Process_Chain_Dist.R loads Markov Chains produced by the Data_Chain.R script for $K \in \{3,4,7\}$ and $C = 8$.  It processes these into two dataframes: Data_abun_df and Data_dist_df.  Data_abun_df contains posterior mean estimaets of $\phi_d$ for LTN-LDA and LDA for all values of $K$.  For each value, it includes a sample ID, a permutation to correct for label switching as $K$ changes, a model ID, and a subcommunity ID.  Data_dist_df contains posterior mean estimates of $\beta_{d,k}$ distributions for LTN-LDA, posterior mean estimates of $\beta_k$ for LTN-LDA, and posterior mean estimates of $\beta_k$ for LDA.  For each value, we record the corresponding the corresponding leaf number, sample number, distribution ID, and subcommunity ID.  The subcommunities for Data_dist_df are permuted so that they correspond to the subcommunities in Data_abun_df.

## Data_Process_Chain_Top_n.R

Data_Process_Chain_Top_n.R loads a Markov chain produced by the Data_Chain.R script for $K = 7$ and $C = 8$.  It processes this into a dataframe: Dist_topn_df.  This contains, for the top $n = 5$ ASVs in each subcommunity according to posterior mean $\beta_k$ estimates for LTN-LDA and LDA, ASV names, ASV prevalence, subcommunity IDs, and model IDs.  The subcommunities are permuted so that they correspond to the same subcommunities in Data_dist_df and Data_abun_df.

## PG_C_Sampler_Strict_Sparse_Tree_Shrink_Perp.cpp

An Rcpp implementation of a modified LTN-LDA Gibbs sampler used to recover $\psi_{d,k}$ and $\phi_d$ on a test set. 

## PG_C_Sampler_Strict_Sparse_Tree_Shrink.cpp

An Rcpp implementation of an LTN-LDA Gibbs sampler.

## Simulation_Chain.R

Simulation_Chain.R generates a training set of $D = 50$ samples and $N_d = 10,000$ reads per sample.  We set $K = 4$, $C = 10$, $\alpha = 1$, $\Lambda = I$, and $\mu = 0$.  The phylogenetic tree is the tree structure stored in tree_mat.rda.  It generates training and test sets according to these parameters.  It then runs an LTN-LDA Gibbs sampler on the training set.  The parameter $K$ is controlled by the taskID parameter, and $C$ is set equal to 10.  The output is a list where each entry contains Markov Chains for different parameters.  All parameters are included except for the Polya-Gamma auxiliary variables.  Note that these Markov chains are not included in the repository because they exceed github's maximum file size.  

## Simulation_Perplexity.R

Simulation_Perplexity.R generates a training set of $D = 50$ samples and $N_d = 10,000$ reads per sample.  We set $K = 4$, $C = 10$, $\alpha = 1$, $\Lambda = I$, and $\mu = 0$.  The phylogenetic tree is the tree structure stored in tree_mat.rda.  It generates training and test sets according to these parameters.  We choose to set the value of the pair $(K,C)$ manually.  It fits a Gibbs sampler for LTN-LDA on the training set, and uses the output of this to estimate the perplexity on the test set using a modified Gibbs sampler.  The output is a vector containing three three perplexity estimates: the LTN-LDA Gibbs sampler, the LDA Gibbs sampler, and the LDA VEM algorithm.

To obtain the perplexity results in the paper, we can run this code over a large number of $(K,C)$ pairs.  

This script also contains optional, commented out, code to find the $L_p$ distances between posterior mean estimates and truth for $\phi_d$, $\beta_k$, and $\beta_{d,k}$.  

## Simulation_Perplexity_Plots.R

Simulation_Perplexity_Plots.R loads the perplexity results obtained from running many different iterations of the Simulation_Cross_Valid.R code.  It then produces three plots.  The first plots perplexity curves for LDA and LTN-LDA as we fix $C$ to truth and vary $K$.  The second plots a perplexity curve for LTN-LDA as we fix $K$ to truth and vary $C$.  The third plots the curves of $L_2$ distances between truth and posterior mean estimates for the $\phi_d$, $\beta_k$, and $\beta_{d,k}$ for LTN-LDA as we fix $K$ to truth and vary $C$.

## Simulation_Plot_Chains.R

Simulation_Plot_Chains.R loads the two dataframes produced by Simulation_Process_Chain_Dist.R.  It uses the dataframes produced by Simulation_Process_Chain_Dist.R to produce a side-by-side plot where the left plot displays posterior mean subcommunity abundance estimates for LTN-LDA and LDA as $K$ varies in $\{4,5,7,10\}$ and the right plot displays subcommunity compositions for LTN-LDA and LDA.  

## Simulation_Process_Chain_Dist.R

Simulation_Process_Chain_Dist.R loads Markov Chains produced by the Data_Chain.R script for $K \in \{4,5,7,10\}$ and $C = 10$.  It processes these into two dataframes: Simulation_abun_df and Simulation_dist_df.  Simulation_abun_df contains posterior mean estimaets of $\phi_d$ for LTN-LDA and LDA for all values of $K$.  For each value, it includes a sample ID, a permutation to correct for label switching as $K$ changes, a model ID, and a subcommunity ID.  Simulation_dist_df contains posterior mean estimates of $\beta_{d,k}$ distributions for LTN-LDA, posterior mean estimates of $\beta_k$ for LTN-LDA, and posterior mean estimates of $\beta_k$ for LDA.  For each value, we record the corresponding the corresponding leaf number, sample number, distribution ID, and subcommunity ID.  The subcommunities for Data_dist_df are permuted so that they correspond to the subcommunities in Simulation_abun_df.

# Data

Stores the two data objects used in the paper.

## abt.rda

abt.rda is a phyloseq object containing the dataset of (Dethlefsen and Relman 2011).  This version was obtained from Kris Sankaran's github.

## tree_mat.rda

A phylogenetic tree structure obtained from modifying the abt.rda dataset.

# Results

Stores results.

## Data

Stores results obtained from analyzing the data of (Dethlefsen and Relman 2011).

### Markov Chains

Stores the dataframes obtained by processing raw Markov chain output.  The Markov chains, unfortuneately, exceeded githubs maximum file size.

#### Data_abun_df

Data_abun_df contains posterior mean estimaets of $\phi_d$ for LTN-LDA and LDA for all values of $K$.  For each value, it includes a sample ID, a permutation to correct for label switching as $K$ changes, a model ID, and a subcommunity ID.  

#### Data_dist_df

Data_dist_df contains posterior mean estimates of $\beta_{d,k}$ distributions for LTN-LDA, posterior mean estimates of $\beta_k$ for LTN-LDA, and posterior mean estimates of $\beta_k$ for LDA.  For each value, we record the corresponding the corresponding leaf number, sample number, distribution ID, and subcommunity ID.  The subcommunities for Data_dist_df are permuted so that they correspond to the subcommunities in Data_abun_df.

#### Data_topn_df

This contains, for the top $n = 5$ ASVs in each subcommunity according to posterior mean $\beta_k$ estimates for LTN-LDA and LDA, ASV names, ASV prevalence, subcommunity IDs, and model IDs.  The subcommunities are permuted so that they correspond to the same subcommunities in Data_dist_df and Data_abun_df.

### Perplexity

Contains files storing perplexity results for the data.

#### LDA_Perp.rda

LDA_Perp.rda contains cross-validation perplexity scores on the data as $K$ varies.

#### Sank_perp_C_8_K.csv

Contains cross-validation perplexity results for LTN-LDA on the data with $C = 8$ and $K$ varying.

#### Sank_perp_results_mat_1.csv

Contains cross-validation perplexity results for LTN-LDA as $K$ varies in $\{2,3,\dots,8\}$ and $C$ varies in $\{1,2,\dots,7\}$.

#### Sank_perp_results_mat_2.csv

Contains cross-validation perplexity results for LTN-LDA as $K$ varies in $\{2,3,\dots,8\}$ and $C$ varies in $\{8,9,\dots,14\}$.

#### Sank_perp_results_mat_3.csv

Contains cross-validation perplexity results for LTN-LDA as $K$ varies in $\{2,3,\dots,8\}$ and $C$ varies in $\{15,16,\dots,21\}$.

## Simulations

Stores results obtained from running numerical simulations.

### Markov Chains

Stores the dataframes obtained by processing raw Markov chain output.  The Markov chains, unfortunately, exceeded githubs maximum file size and were not included.

#### Sim_abun_df

Sim_abun_df contains posterior mean estimaets of $\phi_d$ for LTN-LDA and LDA for all values of $K$.  For each value, it includes a sample ID, a permutation to correct for label switching as $K$ changes, a model ID, and a subcommunity ID.  

#### Sim_dist_df

Sim_dist_df contains posterior mean estimates of $\beta_{d,k}$ distributions for LTN-LDA, posterior mean estimates of $\beta_k$ for LTN-LDA, and posterior mean estimates of $\beta_k$ for LDA.  For each value, we record the corresponding the corresponding leaf number, sample number, distribution ID, and subcommunity ID.  The subcommunities for Data_dist_df are permuted so that they correspond to the subcommunities in Sim_abun_df.

### Perplexity

This directory is blank --- the means and standard errors of the perplexity results are instead recorded in Simulation_Perplexity_Plots.  This directory still exists to record any saved simulation perplexity results.

# References

* Les Dethlefsen and David A. Relman. Incomplete recovery and individualized responses of the human distal gut microbiota to repeated antibiotic perturbation. Proceedings of the National Academy of the Sciences of the United States of America. 18(Supplement 1): 4554-4561, 2011. 
* LeBlanc and Ma

Additionally, I inspiration for coding a collapsed LDA Gibbs sampler from:

* Brooks, Andrew.  "Latent Dirichlet Allocation – under the hood".  Web blog post. data science side projects, thoughts, & experiments.  github.io.  January 17, 2015.  Web.  September 11, 2021.