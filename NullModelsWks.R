######Null models / Randomization tests code##########
# Date: 4/24/2012
# Author: E. M. Hart
# Description: Code to demonstrate and teach 
# null models and Monte Carlo methods.
# Includes plotting methods from standard R libraries
# and ggplot2 libraries
#####################################################

###Example 1
# Simple test of difference of two means
####

###Generate two random data sets###

#set sample size
n <- 15

####Create data frame

s.dat <- data.frame(cbind(c(rnorm(n,20,3),rnorm(n,27,3)),c(rep(1,n),rep(2,n))))
colnames(s.dat) <- c("Measure","Group")
####Stanard R plot
plot(density(s.dat[s.dat[,2]==1,1]),xlim=c(10,40))
lines(density(s.dat[s.dat[,2]==2,1]),col=2)
####ggplot2 code
ggplot(s.dat,aes(x=Measure,group=as.factor(Group),fill=as.factor(Group),alpha=.5))+geom_density()+xlim(0,15)

####Parametric t.test#####
tt1 <- t.test(x=s.dat[s.dat[,2]==2,1],y=s.dat[s.dat[,2]==1,1])

#####Randomization test#######
# In this case the test metric
# is if the difference between the two means
# is different from 0
############

##Set the number of randomizations
nran <- 5000
output <- vector()
for(i in 1:nran){
  ###Shuffle association via sample()
  new.ind <- sample(s.dat$G,length(s.dat$G),replace=F)
  output[i] <- mean(s.dat$M[new.ind==1]) - mean(s.dat$M[new.ind==2])
}

##Find true metric
true.val <- mean(s.dat[s.dat[,2]==1,1])-mean(s.dat[s.dat[,2]==2,1])
###Normal R plot
plot(density(output),xlim=c(-6,6))
abline(v=true.val)
####ggplot2###
ggplot(as.data.frame(output),aes(x=output))+geom_density()+xlim(-6,6)+geom_vline(xintercept=true.val)

###Check true value against quantiles
quantile(output,c(.025,.975))




#########A null model example from Gotelli 2000######
# A simple null model example



###Generate random matrix
spec <- 5
sites <- 7
rand.mat <- rbinom(spec*sites,1,.5)
dim(rand.mat) <- c(spec,sites)

###Function to calculate C-Score our default metric
# See Gotelli 2000 for full description
# or Stone and Roberts 1990
# Input is a species x site data matrix
####

c.score <- function(sp.mat){
  CU <- vector()
  sp.num <- dim(sp.mat)[1]
  count <- 1
  for (i in 1:(sp.num-1)){
    for(k in (i+1):sp.num){
      #first count the number of shared sites
      Q <- sum(apply(rbind(sp.mat[i,],sp.mat[k,]),2,sum)==2)
      CU[count] <- (sum(sp.mat[i,])-Q)*(sum(sp.mat[k,])-Q)
      count<- count+1
    }
  }
  return(mean(CU))
}

#### Randomization algorithm########
#### this randomization algorithm is based on SIM9 from Gotelli 2000
#### In works by keeping row and column totals fixed
#### Inputs are: nsim: the number of simulations to perform
#### sp.mat:  The presence-absence species x sits matrix

sim9 <- function(nsim,sp.mat){
  
  #get constraints
  r.cons <- apply(sp.mat,1,sum)
  c.cons <- apply(sp.mat,2,sum)
  ####Loop through and fill the initial matrix
  outmat <- matrix(0,nrow=length(r.cons),ncol=length(c.cons))
  for(i in 1:length(r.cons)){
    s.ind <- sample(1:length(c.cons),r.cons[i],replace=F)
    outmat[i,s.ind] <- 1
     }
  ### Now row constraints should be the same, but the column counts need to be shuffled
  col.offset<- apply(outmat,2,sum) - c.cons
  for(i in 1:length(col.offset)){
      which()
    
  }
  
}




























