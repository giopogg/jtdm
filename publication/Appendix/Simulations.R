rm(list = ls())
library(mvtnorm)
library(ggplot2)
library(cowplot)

S = 100
n = 100

set.seed(1234)

traits = rmvnorm(S, sigma = matrix(data=c(1,0.8,0.8,1), ncol = 2))
rownames(traits) = paste0("S",1:100)
colnames(traits) = c("t1", "t2")

p1 = ggplot(as.data.frame(traits)) + geom_point(aes(x=t1,y=t2), colour = "red", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = 2, alpha = 0.5) +
  annotate("text", x = 1, y = -2, label = paste0("Correlation = ", signif(cor(traits[,1], traits[,2]), 2)))+
  theme_minimal()
#plot(traits)
#abline(0,1,col = "red")

#####Now build communities.

## Neutral model
# Set 100 communities that are made by 20 species each, drawn randomly with resampling
spp_comm = matrix(0, nrow = n, ncol = S)
colnames(spp_comm) = rownames(traits)

for(i in 1:n){
  comm = sample(colnames(spp_comm),20, replace = T)
  spp_comm[i, names(table(comm))] = table(comm)
}
# Compute community weighted means
CWM = spp_comm %*% traits /rowSums(spp_comm)
p2 = ggplot(as.data.frame(CWM)) + geom_point(aes(x=t1,y=t2), colour = "red", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = 2, alpha = 0.5) +
  annotate("text", x = 0.2, y = -0.4, label = paste0("Correlation = ", signif(cor(CWM[,1], CWM[,2]), 2)))+
  theme_minimal()
#plot(CWM)
#abline(0,1,col = "red")
cor(CWM[,1], CWM[,2])

## Environmental filtering
# The 5 communities are only made of species with traits ina certain range
ranges = c(-2, -1, -0.25, 0.25, 1, 2)
spp_comm = matrix(0, nrow = n, ncol = S)
colnames(spp_comm) = rownames(traits)

for(i in 0:(n-1)){
  comm = sample(which(traits[,1] > ranges[floor(i/20) + 1] & traits[,1] < ranges[floor(i/20) + 2] &
          traits[,2] > ranges[floor(i/20) + 1] & traits[,2] < ranges[floor(i/20) + 2] ), 20, replace = T)
  spp_comm[i + 1, paste0("S",names(table(comm)))] = table(comm)
}

# Compute community weighted means
CWM = spp_comm %*% traits /rowSums(spp_comm)
p3 = ggplot(as.data.frame(CWM)) + geom_point(aes(x=t1,y=t2), colour = "red", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = 2, alpha = 0.5) +
  annotate("text", x = 0.3, y = -1, label = paste0("Correlation = ", signif(cor(CWM[,1], CWM[,2]), 2)))+
  theme_minimal()
#plot(CWM)
#abline(0,1,col = "red")
cor(CWM[,1], CWM[,2])



## Limiting similarity
# Sample 10000 communities that are made by 20 species each, drawn randomly with resampling
# Keep only the 100 communities that maximise trait variance.
spp_comm_all = matrix(0, nrow = 10000, ncol = S)
colnames(spp_comm_all) = rownames(traits)

for(i in 1:10000){
  comm = sample(colnames(spp_comm_all),20, replace = T)
  spp_comm_all[i, names(table(comm))] = table(comm)
}

# Compute variance as the trace of the covariance of observed traits
tot_var = sapply(1:1000, function(i) sum(diag(var(traits[rep(names(spp_comm_all[i,]),spp_comm_all[i,]),]))))

spp_comm =  spp_comm_all[order(tot_var , decreasing = TRUE)[1:100],]

# Compute community weighted means
CWM = spp_comm %*% traits /rowSums(spp_comm)

p4 = ggplot(as.data.frame(CWM)) + geom_point(aes(x=t1,y=t2), colour = "red", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = 2, alpha = 0.5) +
  annotate("text", x = 0.1, y = -0.4, label = paste0("Correlation = ", signif(cor(CWM[,1], CWM[,2]), 2)))+
  theme_minimal()

#plot(CWM)
#abline(0,1,col = "red")
cor(CWM[,1], CWM[,2])

# Final plot
plot_grid(NULL, p1, NULL, p2, p3, p4, ncol = 3,
          labels = c(NA,'a) Traits of the species pool', NA,'   b) CWM traits under Neutral assembly',
                     'c) CMW traits under environmental filtering', '   d) CWM traits under limiting similarity'),
          hjust = 0, label_x = 0.1, label_size = 13)

ggsave("SimulatedCWM.pdf")


