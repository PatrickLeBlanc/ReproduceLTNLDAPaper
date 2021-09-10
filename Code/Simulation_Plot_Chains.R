library(ggplot2)
library(gridExtra)

#load abundance dataframe
load("../Results/Simulations/Markov Chains/Sim_abun_df.rda")

#make first plot, Posterior Mean Estimates of ASV-Subcommunity Distributions
p1 = ggplot(Sim_abun_df, aes(x=samps, y=props)) +
  geom_line(aes(col = Subcommunities)) +
  scale_color_manual( values = c("firebrick","darkgreen","blue","black")) +
  facet_grid(tops ~ model, scale = "free_x") +
  xlab("") +
  ylab("") +
  ggtitle("Posterior Mean Estimates of ASV-Subcommunity Distributions") + 
  theme_classic()

#load distribution dataframe
load("../Results/Simulations/Markov Chains/Sim_dist_df.rda")

#make plots of sample-to-sample heterogeneity in distributions
p2 = ggplot(Sim_dist_df, aes(x=input, y=vals, group=samp, col = Distribution)) +
  geom_line(aes(size = Distribution,col = Distribution)) +
  scale_color_manual( values = c("black","firebrick","skyblue1")) +
  scale_size_manual( values = c(1,1,0.1)) +
  facet_wrap(Sim_dist_df$tops,ncol = 1) +
  xlab("") +
  ylab("") +
  ggtitle("Posterior Mean Estimates of ASV-Subcommunity Distributions") + 
  theme_classic()

#plot both plots side-by-side
grid.arrange(p1,p2,widths = c(2,1.35))