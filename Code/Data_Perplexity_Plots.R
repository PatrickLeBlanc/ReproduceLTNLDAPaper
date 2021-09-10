#load perplexity cross valdiation results
mat_kc_1 = read.csv("../Results/Data/Perplexity/Sank_perp_results_mat_1.csv")
mat_kc_2 = read.csv("../Results/Data/Perplexity/Sank_perp_results_mat_2.csv")
mat_kc_3 = read.csv("../Results/Data/Perplexity/Sank_perp_results_mat_3.csv")

#merge the three datasets into one
mat_kc = cbind(mat_kc_1,mat_kc_2,mat_kc_3)
mat_kc = mat_kc[,- which(colnames(mat_kc) == "X")]
rownames(mat_kc) = 2:8
colnames(mat_kc) = 1:21

#plot perplexity curves for different values of K as C changes
#each curve corresponds to a different avlue of K
matplot(t(mat_kc),type="l",
        xlab = "C",
        ylab = "Perplexity")

#load perplexity results for LDA as K varies
load("../Results/Data/Perplexity/LDA_perp.rda")
perp_k = perp_k[-1]
#load perplexity results for LTN-LDA fixing C at 8 and varying K
LTN_perp_k = read.csv("../Results/Data/Perplexity/Sank_perp_C_8_K.csv")
LTN_perp_k = LTN_perp_k[,-1]

#plot these two curves next to each other
plot(2:42,perp_k[1:41],type="l",col="red",
     ylim = c(min(perp_k,LTN_perp_k),max(perp_k,LTN_perp_k)),
     xlab = "Number of Subcommunities",
     ylab = "Perplexity")
lines(2:42,LTN_perp_k[1:41],type="l",col="blue")
legend(x = 30,y = 26, legend=c("LDA", "LTN LDA"),
       col=c("red", "blue"),lty = 1)
