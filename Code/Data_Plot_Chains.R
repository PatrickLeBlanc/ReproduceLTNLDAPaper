library(ggplot2)
library(gridExtra)

#load top n dataframe
load("../Results/Data/Markov Chains/Data_topn_df.rda")

#plot top n
ggplot(data=Data_topn_df, aes(x=names, y=props, fill = factor(topic))) +
  geom_col(show.legend = FALSE) +
  facet_grid(topic~model,scale = "free_y") +
  coord_flip() +
  xlab("") +
  ylab("") +
  ggtitle("Most Prevalent ASVs in each Community")

#load abundance dataframe
load("../Results/Data/Markov Chains/Data_abun_df.rda")

#make first plot, Posterior Mean Estimates of ASV-Subcommunity Distributions
p1 = ggplot(Data_abun_df, aes(x=samps, y=props)) +
  geom_rect(xmin = 12, xmax = 23, ymin=0,ymax=1,
            fill="grey", size=0.1, alpha=0.1) +
  geom_rect(xmin = 41, xmax = 51, ymin=0,ymax=1,
            fill="grey", size=0.1, alpha=0.1) +
  geom_line(aes(col = Subcommunities)) +
  scale_color_manual( values = c("firebrick","blue","black")) +
  facet_grid(tops ~ model, scale = "free_x") +
  xlab("") +
  ylab("") +
  ggtitle("Posterior Mean Estimates of Subcommunity-Sample Distributions") + 
  theme_classic()

#load distribution dataframe
load("../Results/Data/Markov Chains/Data_dist_df.rda")

#make plots of sample-to-sample heterogeneity in distributions
p2 = ggplot(Data_dist_df, aes(x=input, y=vals, group=samp, col = Distribution)) +
  geom_line(aes(size = Distribution,col = Distribution)) +
  scale_color_manual( values = c("black","firebrick","skyblue1")) +
  scale_size_manual( values = c(1,1,0.1)) +
  facet_wrap(Data_dist_df$tops, ncol = 1) +
  xlab("") +
  ylab("") +
  ggtitle("Posterior Mean Estimates of ASV-Subcommunity Distributions") + 
  theme_classic()

#plot abundance and distributions side-by-side
grid.arrange(p1,p2,widths = c(2,1.35))
