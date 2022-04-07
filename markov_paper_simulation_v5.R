library(ggplot2)
library(expm)
library("gridExtra")
library(cowplot)
library(ggpubr)
library(gridExtra)
library(matrixStats)
# Version 5: Modeling person by person for the Markov Model paper Numerical Simulation
# By: Cassie Berns, 3 April 2022
# The aim of this code is as follows: to first run a Monte Carlo simulation to advance
#   each person in the cohort model as an individual, rather than as a cohort, through
#   four states with given transition probabilities in the matrix p.trans. Once the data
#   on state transitions are obtained for the entire cohort, the mean, variance, and
#   covariances are then calculated.
currpath <- dirname(rstudioapi::callFun("getActiveDocumentContext")$path)  
setwd(currpath)
# Defining Parameters
set.seed(123)
n.cohort <- 1000
time.end <- 30
current_state <- c(n.cohort,0,0,0)
cohort.run <- 10000
time.step <- 1
t.points <- seq(0,time.end,time.step)
time.cycles <- length(t.points)-1
#transition probability matrix: note that it is assumed that once one advances past a certain state value,
#   they cannot backtrack (eg: once someone is in state 3 they cannot move to states 1 or 2).

p <- matrix(c(0.71, 0.1, 0.05, 0.14, 0, 0.76, 0.07, 0.17, 0, 0, 0.89, 0.11, 0, 0, 0, 1), nrow=4, ncol=4, byrow=T)

states<-c(1:4)
numstates<-length(states)
p1<-c(1,0,0,0)
trace<-matrix(1,nrow=n.cohort,ncol=time.cycles+1)
cohort.sample <- array(c(rep(0,cohort.run*4*(time.cycles+1))),dim=c(cohort.run,4,time.cycles+1))
s1.mean <- array(c(rep(0,time.cycles+1)),dim=c(time.cycles+1))
s1.sd <- array(c(rep(0,time.cycles+1)),dim=c(time.cycles+1))
s2.mean <- array(c(rep(0,time.cycles+1)),dim=c(time.cycles+1))
s2.sd <- array(c(rep(0,time.cycles+1)),dim=c(time.cycles+1))
s3.mean <- array(c(rep(0,time.cycles+1)),dim=c(time.cycles+1))
s3.sd <- array(c(rep(0,time.cycles+1)),dim=c(time.cycles+1))
s4.mean <- array(c(rep(0,time.cycles+1)),dim=c(time.cycles+1))
s4.sd <- array(c(rep(0,time.cycles+1)),dim=c(time.cycles+1))
start.time <-proc.time()
for (i in 1:cohort.run){
  for(t in 2:(time.cycles+1)){
    prev.state<-trace[,t-1]
    sample.time <- proc.time()
    probability <-p[prev.state,]
    usample<-runif(n.cohort,0,1)
    sum.p <-rowCumsums(probability)
    state <- max.col(sum.p >=usample,ties.method = "first")
    trace[,t] <- state
    cohort.sample[i,1,t]=count(trace[,t]==1)
    cohort.sample[i,2,t]=count(trace[,t]==2)
    cohort.sample[i,3,t]=count(trace[,t]==3)
    cohort.sample[i,4,t]=count(trace[,t]==4)
  }
  # cohort.sample[i,2,]=sum(trace[which(trace[,t]==2),t])
  # cohort.sample[i,3,]=sum(trace[which(trace[,t]==3),t])
  # cohort.sample[i,4,]=sum(trace[which(trace[,t]==4),t])
  # s1.mean[i,t] <- mean(cohort.sample[,1,t])
  # s1.sd[i,t] <- sd(cohort.sample[,1,t])
  # 
  # s2.mean[i,t] <- mean(cohort.sample[,2,t])
  # s2.sd[i,t] <- sd(cohort.sample[,2,t])
  # 
  # s3.mean[i,t] <- mean(cohort.sample[,3,t])
  # s3.sd[i,t] <- sd(cohort.sample[,3,t])
  # 
  # s4.mean[i,t] <- mean(cohort.sample[,4,t])
  # s4.sd[i,t] <- sd(cohort.sample[,4,t])
}
elapsed.time <-proc.time()-start.time
print(elapsed.time)
for(t in 2:(time.cycles+1)){
  s1.mean[t] <- mean(cohort.sample[,1,t])
  s1.sd[t] <- sd(cohort.sample[,1,t])

  s2.mean[t] <- mean(cohort.sample[,2,t])
  s2.sd[t] <- sd(cohort.sample[,2,t])

  s3.mean[t] <- mean(cohort.sample[,3,t])
  s3.sd[t] <- sd(cohort.sample[,3,t])

  s4.mean[t] <- mean(cohort.sample[,4,t])
  s4.sd[t] <- sd(cohort.sample[,4,t])
}
s1.mean[1] <- n.cohort
s1.mean.micro <- as.data.frame(cbind(t.points,s1.mean))
s1.mean.high.micro <- as.data.frame(cbind(t.points,s1.mean+s1.sd))
s1.mean.low.micro <- as.data.frame(cbind(t.points,s1.mean-s1.sd))

s2.mean.micro <- as.data.frame(cbind(t.points,s2.mean))
s2.mean.high.micro <- as.data.frame(cbind(t.points,s2.mean+s2.sd))
s2.mean.low.micro <- as.data.frame(cbind(t.points,s2.mean-s2.sd))

s3.mean.micro <- as.data.frame(cbind(t.points,s3.mean))
s3.mean.high.micro <- as.data.frame(cbind(t.points,s3.mean+s3.sd))
s3.mean.low.micro <- as.data.frame(cbind(t.points,s3.mean-s3.sd))

s4.mean.micro <- as.data.frame(cbind(t.points,s4.mean))
s4.mean.high.micro <- as.data.frame(cbind(t.points,s4.mean+s4.sd))
s4.mean.low.micro <- as.data.frame(cbind(t.points,s4.mean-s4.sd))

colnames(s1.mean.micro) <- c("time","count")
colnames(s1.mean.high.micro) <- c("time","count")
colnames(s1.mean.low.micro) <- c("time","count")

colnames(s2.mean.micro) <- c("time","count")
colnames(s2.mean.high.micro) <- c("time","count")
colnames(s2.mean.low.micro) <- c("time","count")

colnames(s3.mean.micro) <- c("time","count")
colnames(s3.mean.high.micro) <- c("time","count")
colnames(s3.mean.low.micro) <- c("time","count")

colnames(s4.mean.micro) <- c("time","count")
colnames(s4.mean.high.micro) <- c("time","count")
colnames(s4.mean.low.micro) <- c("time","count")

multinom_means_func <- function(p) {
  multinom_means=array(rep(0,4*(time.cycles+1)),dim=c(time.cycles+1,4))
  for (i in 2:(time.cycles+1)) {
    pt = p%^%(i-1)
    multinom_means[i,] = colSums((pt*current_state))
  }
  
  return(multinom_means)
}
multinom_means=multinom_means_func(p)
multinom_means[1,] <- c(n.cohort,0,0,0)

multinom_sd_func <- function(p) {
  multinom_sd = array(rep(0,4*(time.cycles+1)),dim=c(time.cycles+1,4))
  for (k in 2:(time.cycles+1)) {
    probs <- multinom_means[k,]/n.cohort
    #nval <- final_multinom_means[k,]
    var <- n.cohort * (1-probs) * probs
    multinom_sd[k,] <- sqrt(var)
  }
  return(multinom_sd)
}

#standard deviations for model values 
multinom_sd <- multinom_sd_func(p)
multinom_sd[1,] <- c(0,0,0,0)

s1.mean.multi <- as.data.frame(cbind(t.points,multinom_means[,1]))
s1.mean.high.multi <- as.data.frame(cbind(t.points,multinom_means[,1]+multinom_sd[,1]))
s1.mean.low.multi <- as.data.frame(cbind(t.points,multinom_means[,1]-multinom_sd[,1]))

s2.mean.multi <- as.data.frame(cbind(t.points,multinom_means[,2]))
s2.mean.high.multi <- as.data.frame(cbind(t.points,multinom_means[,2]+multinom_sd[,2]))
s2.mean.low.multi <- as.data.frame(cbind(t.points,multinom_means[,2]-multinom_sd[,2]))

s3.mean.multi <- as.data.frame(cbind(t.points,multinom_means[,3]))
s3.mean.high.multi <- as.data.frame(cbind(t.points,multinom_means[,3]+multinom_sd[,3]))
s3.mean.low.multi <- as.data.frame(cbind(t.points,multinom_means[,3]-multinom_sd[,3]))

s4.mean.multi <- as.data.frame(cbind(t.points,multinom_means[,4]))
s4.mean.high.multi <- as.data.frame(cbind(t.points,multinom_means[,4]+multinom_sd[,4]))
s4.mean.low.multi <- as.data.frame(cbind(t.points,multinom_means[,4]-multinom_sd[,4]))

colnames(s1.mean.multi) <- c("time","count")
colnames(s1.mean.high.multi) <- c("time","count")
colnames(s1.mean.low.multi) <- c("time","count")

colnames(s2.mean.multi) <- c("time","count")
colnames(s2.mean.high.multi) <- c("time","count")
colnames(s2.mean.low.multi) <- c("time","count")

colnames(s3.mean.multi) <- c("time","count")
colnames(s3.mean.high.multi) <- c("time","count")
colnames(s3.mean.low.multi) <- c("time","count")

colnames(s4.mean.multi) <- c("time","count")
colnames(s4.mean.high.multi) <- c("time","count")
colnames(s4.mean.low.multi) <- c("time","count")

# c.color <- c(multinom_mean="red",
#              multinom_highsd="red",
#              multinom_lowsd="red",
#              micro_mean="#1E90FF" ,
#              micro_highsd="#1E90FF",
#              micro_lowsd="#1E90FF")

c.color <- c(multinom_mean="red",
             micro_mean="#0000FF")
x.breaks <- seq(0,time.end,by=5)
y.breaks <- seq(0,n.cohort,by=200)
ltype1 <- "solid"
ltype2 <- "dotted"

p1 <- ggplot()+   
  geom_line(data = s1.mean.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5, linetype=ltype1) +
  geom_line(data = s1.mean.high.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5,linetype=ltype1) +
  geom_line(data = s1.mean.low.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5,linetype=ltype1) +
  #geom_point(data = s1.mean.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.9, size=0.7, shape=3) +
  #geom_point(data = s1.mean.high.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.9, size=0.7,shape=3) +
  #geom_point(data = s1.mean.low.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.9, size=0.7,shape=3) +
  geom_line(data = s1.mean.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1,linetype=ltype2) +
  geom_line(data = s1.mean.high.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1, linetype=ltype2) +
  geom_line(data = s1.mean.low.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1,linetype=ltype2) +
  # geom_line(data = median_lower.df, aes(x=x.values, y=V2, color = "c3_minmax_median"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = std_upper.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = std_lower.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = median_std_upper.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=0.7,linetype=ltype) +
  # geom_line(data = median_std_lower.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=0.7,linetype=ltype) +
#geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3)+
scale_x_continuous(name = "Time (years)",breaks=x.breaks,expand=c(0,0),limits=c(0,31)) +
  scale_y_continuous(name = expression(Number~of~individuals~S[1]),breaks=y.breaks,expand=c(0,0),limits=c(0, n.cohort)) +
  scale_color_manual(values = c.color,
                     name = "Approaches")+
  #ggtitle("") +
  theme(legend.position = "none",
        axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.title = element_text(color = "black", size = 17),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(color = "black", size = 16)) +
  guides(color = guide_legend(title="CDF types"), linetype = guide_legend(title="CDF types"))
ggsave(paste(currpath,"/S1.tiff",sep="") , units="in", width=7, height=4, dpi=200, compression = 'lzw')

p2 <- ggplot()+   geom_line(data = s2.mean.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5, linetype=ltype1) +
  #geom_ribbon(aes(x=x, ymin=minmax_lower.df$V2, ymax=minmax_upper.df$V2),alpha=0.8,fill="#DCDCDC")+
  geom_line(data = s2.mean.high.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5,linetype=ltype1) +
  geom_line(data = s2.mean.low.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5,linetype=ltype1) +
  geom_line(data = s2.mean.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1,linetype=ltype2) +
  geom_line(data = s2.mean.high.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1, linetype=ltype2) +
  geom_line(data = s2.mean.low.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1,linetype=ltype2) +
  # geom_line(data = median_lower.df, aes(x=x.values, y=V2, color = "c3_minmax_median"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = std_upper.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = std_lower.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = median_std_upper.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=0.7,linetype=ltype) +
  # geom_line(data = median_std_lower.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=0.7,linetype=ltype) +
  #geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3)+
  scale_x_continuous(name = "Time (years)",breaks=x.breaks,expand=c(0,0),limits=c(0,31)) +
  scale_y_continuous(name = expression(Number~of~individuals~S[2]),breaks=y.breaks,expand=c(0,0),limits=c(0, n.cohort)) +
  scale_color_manual(values = c.color,
                     name = "Approaches")+
  #ggtitle("") +
  theme(legend.position = "none",
        axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.title = element_text(color = "black", size = 17),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(color = "black", size = 16)) +
  guides(color = guide_legend(title="CDF types"), linetype = guide_legend(title="CDF types"))
ggsave(paste(currpath,"/S2.tiff",sep="") , units="in", width=7, height=4, dpi=200, compression = 'lzw')

p3 <- ggplot()+   geom_line(data = s3.mean.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5, linetype=ltype1) +
  #geom_ribbon(aes(x=x, ymin=minmax_lower.df$V2, ymax=minmax_upper.df$V2),alpha=0.8,fill="#DCDCDC")+
  geom_line(data = s3.mean.high.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5,linetype=ltype1) +
  geom_line(data = s3.mean.low.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5,linetype=ltype1) +
  geom_line(data = s3.mean.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1,linetype=ltype2) +
  geom_line(data = s3.mean.high.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1, linetype=ltype2) +
  geom_line(data = s3.mean.low.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1,linetype=ltype2) +
  # geom_line(data = median_lower.df, aes(x=x.values, y=V2, color = "c3_minmax_median"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = std_upper.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = std_lower.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = median_std_upper.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=0.7,linetype=ltype) +
  # geom_line(data = median_std_lower.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=0.7,linetype=ltype) +
  #geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3)+
  scale_x_continuous(name = "Time (years)",breaks=x.breaks,expand=c(0,0),limits=c(0,31)) +
  scale_y_continuous(name = expression(Number~of~individuals~S[3]),breaks=y.breaks,expand=c(0,0),limits=c(0, n.cohort)) +
  scale_color_manual(values = c.color,
                     name = "Approaches")+
  #ggtitle("") +
  theme(legend.position = "none",
        axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.title = element_text(color = "black", size = 17),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(color = "black", size = 16)) +
  guides(color = guide_legend(title="CDF types"), linetype = guide_legend(title="CDF types"))
ggsave(paste(currpath,"/S3.tiff",sep="") , units="in", width=7, height=4, dpi=200, compression = 'lzw')

p4 <- ggplot()+   geom_line(data = s4.mean.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5, linetype=ltype1) +
  #geom_ribbon(aes(x=x, ymin=minmax_lower.df$V2, ymax=minmax_upper.df$V2),alpha=0.8,fill="#DCDCDC")+
  geom_line(data = s4.mean.high.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5,linetype=ltype1) +
  geom_line(data = s4.mean.low.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.7, size=0.5,linetype=ltype1) +
  geom_line(data = s4.mean.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1,linetype=ltype2) +
  geom_line(data = s4.mean.high.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1, linetype=ltype2) +
  geom_line(data = s4.mean.low.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=1, size=1,linetype=ltype2) +
  # geom_line(data = median_lower.df, aes(x=x.values, y=V2, color = "c3_minmax_median"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = std_upper.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = std_lower.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=0.7,linetype=ltype2) +
  # geom_line(data = median_std_upper.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=0.7,linetype=ltype) +
  # geom_line(data = median_std_lower.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=0.7,linetype=ltype) +
  #geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3)+
  scale_x_continuous(name = "Time (years)",breaks=x.breaks,expand=c(0,0),limits=c(0,31)) +
  scale_y_continuous(name = expression(Number~of~individuals~S[4]),breaks=y.breaks,expand=c(0,0),limits=c(0, n.cohort)) +
  scale_color_manual(values = c.color,
                     name = "Approaches")+
  #ggtitle("") +
  theme(legend.position = "none",
        axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.title = element_text(color = "black", size = 17),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(color = "black", size = 16)) +
  guides(color = guide_legend(title="CDF types"), linetype = guide_legend(title="CDF types"))
ggsave(paste(currpath,"/S4.tiff",sep="") , units="in", width=7, height=4, dpi=200, compression = 'lzw')

# 
# #making new plot to use for legend
ggpl_legend <- ggplot()+   geom_line(data = s4.mean.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=1, size=1, linetype=ltype1) +
  #geom_ribbon(aes(x=x, ymin=minmax_lower.df$V2, ymax=minmax_upper.df$V2),alpha=0.8,fill="#DCDCDC")+
  geom_line(data = s4.mean.high.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.9, size=1,linetype=ltype1) +
  geom_line(data = s4.mean.low.micro, aes(x=t.points, y=count, color = "micro_mean"), alpha=0.9, size=1,linetype=ltype1) +
  geom_line(data = s4.mean.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=0.9, size=1,linetype=ltype2) +
  geom_line(data = s4.mean.high.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=0.9, size=1, linetype=ltype2) +
  geom_line(data = s4.mean.low.multi, aes(x=t.points, y=count, color = "multinom_mean"), alpha=0.9, size=1,linetype=ltype2) +
  # geom_line(data = median_lower.df, aes(x=x.values, y=V2, color = "c3_minmax_median"), alpha=0.9, size=1,linetype=ltype2) +
  # geom_line(data = std_upper.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=1,linetype=ltype2) +
  # geom_line(data = std_lower.df, aes(x=x.values, y=V2, color = "c5_minmax_mean_std"), alpha=0.9, size=1,linetype=ltype2) +
  # geom_line(data = median_std_upper.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=1,linetype=ltype) +
  # geom_line(data = median_std_lower.df, aes(x=x.values, y=V2, color = "c6_full"), alpha=0.9, size=1,linetype=ltype) +
  #geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3)+
  scale_x_continuous(name = "Time (years)",breaks=x.breaks,expand=c(0,0),limits=c(0,31)) +
  scale_y_continuous(name = expression(Number~of~individuals~S[4]),breaks=y.breaks,expand=c(0,0),limits=c(0, n.cohort)) +
  scale_color_manual(values = c.color,
                     name = "Approaches",labels=c("Multinomial distribution","Microsimulation"))+
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"),
        legend.title = element_text(color = "black", size = 17,face="bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white"),
        legend.box.margin = margin(5),
        legend.text = element_text(color = "black", size = 16,face="bold")) +
  guides(color = guide_legend(title="Approaches"), linetype = guide_legend(title="Approaches"))
scale_fill_discrete(name="Approaches", labels=c("Microsimulation","Multinomial distribution"))

#function to extract legend
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
# #creating the new legend
shared_legend <- extract_legend(ggpl_legend)
# 
# plotlist = list(p0, p3, p1, p4, p2, p5)
# ml <- marrangeGrob(plotlist, nrow=2, ncol=3,top=NULL)
# 
# ml <-grid.arrange(arrangeGrob(p0, p1, p2, p3, p4, p5, nrow=2, ncol=3),
#                   shared_legend, nrow = 2, ncol=1, heights = c(16, 2))
# 
# ggsave(paste(currpath,"/plot/input_feb2022/multi_simple_same_shade.tiff",sep=""),ml, units="in", width=16, height=8, dpi=300, compression = 'lzw')
# 
# ml <- ggarrange(p1, p2, p3, p4,ncol=2, nrow=2,shared_legend)

ml <-grid.arrange(arrangeGrob(p1, p2, p3, p4, nrow=2, ncol=2),
                  shared_legend, nrow = 2, ncol=1, heights = c(16, 2))

ggsave(paste(currpath,"/multi_final.tiff",sep=""),ml, units="in", width=16, height=8, dpi=300, compression = 'lzw')
