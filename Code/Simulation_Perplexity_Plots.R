######
# Set K = 4, C = 10
#Vary modelled K
# Find Perpleixy reuslts on test set

#mean perplexity results
ltn_mean = c(18.85588, 18.57939, 18.28391, 18.28295, 18.28338, 18.28417,  18.27827,18.28633,18.28561,18.28602,18.27237,18.28669,18.28677,18.28720,18.27380)
gibbs_mean =  c(22.02326, 21.08663, 20.18087, 19.91968, 19.67414, 19.48571, 19.32166,19.19323,19.06823,18.97301,18.86739,18.80594,18.73669,18.67695,18.60863)
vem_mean = c(22.10172, 21.17243, 20.32376, 20.03573, 19.82287, 19.67847, 19.52686,19.42411,19.33230,19.26611,19.18195,19.1915190,19.10025,19.08471,19.03745)
diff_mean = c(-3.167385, -2.507233, -1.896959, -1.636729, -1.390758, -1.201541, -1.043386,-0.9069052,-0.7826182,-0.6869926,-0.590123,-0.5192471,-0.4499262,-0.3897482,-0.3348277)

#standard errors
ltn_se = c(0.2425307, 0.2394357, 0.2366763, 0.2368156, 0.236882, 0.2368877, 0.2394086,0.236728,0.2368847,0.236911,0.2385101,0.2369245,0.2369199,0.2369105,0.2385367)
gibbs_se = c(0.3000582, 0.2878427, 0.2793704, 0.2740694, 0.2693248, 0.2676298, 0.269996,0.2656374,0.2649179,0.2635284,0.2630384,0.2600536,0.259374,0.2558966,0.2553406)
vem_se = c(0.3033448, 0.2917799, 0.2855275, 0.2773037, 0.2718938, 0.2715212, 0.2721278,0.267465,0.2665721,0.2659874,0.265572,0.2639662,0.2620169,0.2609801,0.2612784)
diff_se = c(0.07454353, 0.06378899, 0.05707894, 0.04858733, 0.04331082, 0.04086327, 0.03914101,0.03810839,0.03511972,0.0334383,0.03076392,0.02844463,0.02616278,0.02316887,0.02108921)

#the range of subcommunities
tops = 2:16

#find range
bottom = min(c(ltn_mean,gibbs_mean,vem_mean) - max(c(ltn_se,gibbs_se,vem_se)))
top = max(c(ltn_mean,gibbs_mean,vem_mean) + max(c(ltn_se,gibbs_se,vem_se)))

#plot LDA Gibbs and LTN-LDA perplexity results as K changes
plot(tops, ltn_mean, type="l",
     ylim = c(bottom, top),
     xlab = "Modelled Number of Subcommunities",
     ylab = "Perplexity",
     col="blue")
points(tops, ltn_mean,pch=19,col="blue")
arrows(tops, ltn_mean - ltn_se, tops, ltn_mean + ltn_se,
       length=0.05, angle=90, code=3,col="blue")
lines(tops,gibbs_mean,col="red")
points(tops,gibbs_mean,pch=19,col="red")
arrows(tops, gibbs_mean-gibbs_se, tops, gibbs_mean + gibbs_se,
       length=0.05, angle=90, code=3,col="red")
legend(6.5, 22, legend=c("LDA", "LTN LDA"),
       col=c("red", "blue"), lty=1, cex=0.8)


######
# Set K = 4, C = 5
#Vary modelled threshold
# Find Perpleixy reuslts on test set


#mean resutls for LTN-LDA perplexity
ltn_mean = c(15.7306, 15.29598, 14.80268, 14.80798, 14.81480, 14.82844, 14.78313, 14.86270, 14.79338, 14.83189, 14.84519, 14.80821, 14.78866, 14.80821, 14.80985, 14.80985, 14.80985, 14.80985, 14.80985)
#standard errors
ltn_se = c(0.3208034, 0.3121336, 0.297144, 0.2927624, 0.2935443, 0.2920883, 0.3033424, 0.2997002, 0.3009515, 0.2946937, 0.2948669, 0.29275, 0.2967749,0.29275,0.2928447,0.2928447,0.2928447,0.2928447,0.2928447)

#set domain values for C
thresh = c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37)

#find range
bottom = min(c(ltn_mean) - max(c(ltn_se)) - 1)
top = max(c(ltn_mean) + max(c(ltn_se)) + 1)

#plot
plot(thresh, ltn_mean, type="l",
     ylim = c(bottom, top),
     xlab = "Modelled Cutoff C",
     ylab = "Perplexity",
     col="blue")
points(thresh, ltn_mean,pch=19,col="blue")
arrows(thresh, ltn_mean - ltn_se, thresh, ltn_mean + ltn_se,
       length=0.05, angle=90, code=3,col="blue")


#############
#Set K = 4, C = 10
# Vary modelled C
#find L_p distance between true and posterior mean distributions

# mean L_p distance
phi = c(1.8430024,1.3103418,0.2555306,0.3685165,0.3632270,2.2327761,rep(2.7628544,4),2.8441400)
kdv = c(4.7155191,3.4634802,0.9094481,1.0548659,1.0494419,5.0848770,rep(6.2498393,4),6.3777029)
kv = c(0.5085376,0.3592188,0.1062151,0.1230520,0.1210231,0.5973031,rep(0.6967518,4),0.7070941)

#standard error
phi_sd = c(0.2075234,0.2007645,0.08971663,0.1222815,0.1289433,0.0908768,rep(0.04787115,4),0.04344717)
kdv_sd = c(0.3255159,0.323039,0.1415926,0.1723617,0.1771282,0.1790466,rep(0.1196509,4),0.1161416)
kv_sd = c(0.0556162,0.0494964,0.01907214,0.02252207,0.02349504,0.02937618,rep(0.02050077,4),0.01876762)

#domain range for C
thresh = c(0,5,10,15,20,25,30,35,40,45,50)

#find range
bottom = min(c(phi,kdv,kv) - max(c(phi_sd,kdv_sd,kv_sd)) - 1)
bottom = 0
top = max(c(phi,kdv,kv) + max(c(phi_sd,kdv_sd,kv_sd)) + 0.5)

#plot
plot(thresh, phi, type="l",
     ylim = c(bottom, top),
     xlab = "Modelled Cutoff C",
     ylab = "L_2 Distance",
     main = "L_2 distance between posterior estimates and truth",
     col="blue")
points(thresh, phi,pch=19,col="blue")
arrows(thresh, phi-phi_sd, thresh, phi+phi_sd,
       length=0.05, angle=90, code=3,col="blue")
lines(thresh,kdv,col="red")
points(thresh, kdv,pch=19,col="red")
arrows(thresh, kdv-kdv_sd, thresh, kdv+kdv_sd,
       length=0.05, angle=90, code=3,col="red")
lines(thresh,kv,col="black")
points(thresh, kv,pch=19,col="black")
arrows(thresh, kv-kv_sd, thresh, kv+kv_sd,
       length=0.05, angle=90, code=3,col="black")
legend(5, 7, legend=c("Sub-Samp Dist", "Samp ASV-Sub Dist", "Avg ASV-Sub Dist"),
       col=c("red", "blue","black"), lty=1, cex=0.8)
