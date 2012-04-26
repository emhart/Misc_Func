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
#### Inputs are:
#### sp.mat:  The presence-absence species x sits matrix
#### This is a brute force algorithm that uses a while loop
#### it has a fail safe for the loop, but may output bad results
### although it hasn't yet.  Also note that it can handle matrix saturation
### up to about 40% (40% of all values are 1's), after that it begins to seriously lag

sim9<-function(sp.mat){
r.size <- dim(sp.mat)[1]
c.size <- dim(sp.mat)[2]
##Define output matrix
outmat <- matrix(0,nrow=r.size,ncol=c.size)  
#get constraints
r.cons <- apply(sp.mat,1,sum)
c.cons <- apply(sp.mat,2,sum)
# create index list from constraints
r.i <- rep(1:r.size,r.cons)
c.i <- rep(1:c.size,c.cons)
# Randomize our vectors
index <- 1:length(r.i)
r.ind <- sample(1:length(r.i),length(r.i),replace=F)
c.ind <- sample(1:length(c.i),length(c.i),replace=F)
# create our first loop
nmat <- cbind(r.i[r.ind],c.i[c.ind])

# figure out how many duplicates we have
dup.i <- index[duplicated(nmat)]

count <- 1
# loop until all duplicates are removed
while(length(dup.i > 0) || count > 100){

# loop through each duplicate swapping out
# rows and columns
for(i in 1:length(dup.i)){
oldval.r <- nmat[dup.i[i],1]
oldval.c <- nmat[dup.i[i],2]

new.ind.r <- sample(r.ind,1)
new.ind.c <- sample(c.ind,1)

newval.r <- nmat[new.ind.r,1]
newval.c <- nmat[new.ind.c,2]


nmat[new.ind.r,1] <- oldval.r
nmat[new.ind.c,2] <- oldval.c


nmat[dup.i[i],1]<-newval.r
nmat[dup.i[i],2]<-newval.c

}
# recalculate number of duplicates
dup.i <- index[duplicated(nmat)]
# adjust failsafe counter
count <- count + 1
}
# set new matrix
outmat[nmat] <- 1

return(outmat)}


#####Just a bit of code to verify that the
#### randomization is working
### Uncomment to test
#r.c <- apply(rand.mat,1,sum)
#c.c <- apply(rand.mat,2,sum)
#test <- sim9(rand.mat)
#r.t <- apply(test,1,sum)
#c.t <- apply(test,2,sum)
#sum(c(r.c-r.t,c.c-c.t))

####randomization

###Generate random matrix
spec <- 10
sites <- 15
rand.mat <- rbinom(spec*sites,1,.35)
dim(rand.mat) <- c(spec,sites)

# set the number of simulations
# 5000 is a good compromise of speed and sample size
n.sim <- 5000
###calculate true value
t.val <- c.score(rand.mat)

prog <- txtProgressBar(min=0, max=n.sim, char="*", style=3)
null.vec <- vector()
for(i in 1:n.sim){
  tmp <-try(c.score(sim9(rand.mat)),silent=T)
  if(is.numeric(tmp)){ null.vec[i] <- tmp }
   ###Code for a progress bar
  setTxtProgressBar(prog, i)
  
}

hist(null.vec)
abline(v=t.val,col=2,lwd=2)
##Check the quantiles
quantile(null.vec,c(.025,.975),na.rm=T)

###Now use real data from Gotelli 2000, Virginia Ants
##Get vaAnts.csv from the github file
ants <- read.csv("vaAnts.csv")

###calculate true value
t.val <- c.score(ants)

ant.vec <- vector()

prog <- txtProgressBar(min=0, max=n.sim, char="*", style=3)
null.vec <- vector()
for(i in 1:n.sim){
    tmp <-try(c.score(sim9(ants)),silent=T)
   if(is.numeric(tmp)){ ant.vec[i] <- tmp }
  ###Code for a progress bar
  
  setTxtProgressBar(prog, i)
  
}

hist(ant.vec)
abline(v=t.val,col=2,lwd=2)
quantile(ant.vec,c(.025,.975),na.rm=T)

###########Abundance models coming soon#########
