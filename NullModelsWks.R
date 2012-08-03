######Null models / Randomization tests code##########
# Date: 4/24/2012
# Author: E. M. Hart
# Description: Code to demonstrate and teach 
# null models and Monte Carlo methods.
# Includes plotting methods from standard R libraries
# and ggplot2 libraries
# requires two data files
# available at https://github.com/emhart/Misc_Func
# they are vaAnts.csv or aussieAnts.csv
#####################################################

###Example 1
# Simple test of difference of two means
####

#############Note I do some plots in regular and again in ggplot, so install ggplot if you want to try it.
#install.packages("ggplot2")
library(ggplot2)
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
ggplot(s.dat,aes(x=Measure,group=as.factor(Group),fill=as.factor(Group),alpha=.5))+geom_density()

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
plot(density(output),xlim=c(-10,10))
abline(v=true.val,col=2)
####ggplot2###
ggplot(as.data.frame(output),aes(x=output))+geom_density()+xlim(-10,10)+geom_vline(xintercept=true.val)

###Check true value against quantiles
quantile(output,c(.025,.975))



######### Boot strap linear regression to accont for non-independence
######### Here we test for signifigance of the slope
###Simulate regression data

x <- rnorm(80,10,3)
###You can adjust the parameter in the term below
### rnorm(80,0,4.5), increase the sd to change the effect
### of sample size.  As noise decreases, the importance of sample
### size also decreases
y <- cbind(rep(1,80),x)%*%c(.1,.6) + rnorm(80,0,4.5)
### the number of independent sites is 20
### now randomly sort them into groups
g.index <- vector()
for(i in 1:4){
  g.index <- c(g.index,sample(1:20,round(runif(1,10,20)),replace=F))
}
x <- x[1:length(g.index)]
y <- y[1:length(g.index)]
plot(x,y)
###Now check the full model without the bootstrap
summary(lm(y~x))
index<- 1:length(g.index)
###Create matrix to hold data
outmat <- matrix(0,ncol=2,nrow=n.sim)
###Code for progress bar
prog <- txtProgressBar(min=0, max=n.sim, char="*", style=3)

for(i in 1:n.sim){
    #define new resampled index
    rsamp.i <- vector()
    ##Loop through and randomly sample 1 instance from each site
    for(k in 1:length(unique(g.index))){
      ###Control for behavior of the sample function
      if(sum(g.index==k)>1){rsamp.i<- c(rsamp.i,sample(which(g.index==k),1))}
      if(sum(g.index==k)==1){rsamp.i <- c(rsamp.i,which(g.index==k))}
      
    }
    setTxtProgressBar(prog, i)
  t.mod <- lm(y[rsamp.i]~x[rsamp.i])
  outmat[i,] <- coef(t.mod)
    
}
hist(outmat[,2])
quantile(outmat[,2],c(.025,.975))

###Plot the full slope and then the boot strap slope in red

plot(x,y)
abline(lm(y~x))
abline(apply(outmat,2,mean),col=2)



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

########### Abundance null models ########
### First select a metric ####
### In this case we use U ####

U.calc <- function(sp.mat){
  V <- var(apply(sp.mat,1,sum))
  W <- sum(apply(sp.mat,2,var))
  return(V/W)
}

###Function to draw from the log normal
r.Ni <- function(n,a){
  return(exp((rnorm(n)/(2*a))))
}

### Converts a matrix to a vector ###

find.matrix.index <- function(x.coord,y.coord,max.x){
  return.index <- vector()
  return.index <- ((y.coord -1)*max.x) + x.coord
  return(return.index)
}

###generate random matrix following Gotelli and Ulrich 2010

sites <- round(runif(1,5,50))
species <- round(runif(1,10,200))

sp.totals <- sort(r.Ni(species,.5),decreasing=T)
### Mimic veil line
Smax <- sample((species/2):species,1)
sp.totals <- sp.totals[1:Smax]

rand.mat <- matrix(0,ncol=sites,nrow=Smax)
ccA <- runif(sites,0,1)
for(i in 1:Smax){
rand.mat[i,] <- t(rmultinom(1,exp(sp.totals[i]),prob=ccA))
}



####Randomization IT, randomization 5 from Gotelli and Ulrich 2010

IT <- function(sp.mat){
  #get the total number of individuals
  i.t <- sum(sp.mat)
  ###calculate marginal probabilities for r and c
  r.p <- apply(sp.mat,1,sum)/i.t
  c.p <-  apply(sp.mat,2,sum)/i.t
  
  ####draw rows and colums
  new.row <- sample(1:length(r.p),i.t,replace=T,prob=r.p)
  new.col <- sample(1:length(c.p),i.t,replace=T,prob=c.p)
  ###Get vector positions for each now 
  vec.pos <- find.matrix.index(new.col,new.row,length(c.p))
  ###now sum all the counts
  my.counts <- table(vec.pos)
  #expand the count vector back into the right size for the matrix
  outmat <- rep(0,prod(dim(sp.mat)))
  outmat[as.numeric(names(my.counts))] <- my.counts
  dim(outmat) <- c(dim(sp.mat)[2],dim(sp.mat)[1])
  
   return(t(outmat))
}


t.val <- U.calc(rand.mat)
n.sim <- 5000
abun.rand<- vector()

prog <- txtProgressBar(min=0, max=n.sim, char="*", style=3)
null.vec <- vector()
for(i in 1:n.sim){
  tmp <-try(U.calc(IT(rand.mat)),silent=T)
  if(is.numeric(tmp)){ abun.rand[i] <- tmp }
  ###Code for a progress bar
  setTxtProgressBar(prog, i)
  }

hist(abun.rand)
abline(v=t.val,col=2,lwd=2)
quantile(abun.rand,c(.025,.975),na.rm=T)

#########Now using real data from 
## Gotelli and Ulrich ecological archives
ausAnts <- read.csv("aussieAnts.csv")
t.val <- U.calc(ausAnts)
n.sim <- 5000
ants.vec<- vector()

prog <- txtProgressBar(min=0, max=n.sim, char="*", style=3)
for(i in 1:n.sim){
  tmp <-try(U.calc(IT(ausAnts)),silent=T)
  if(is.numeric(tmp)){ ants.vec[i] <- tmp }
  ###Code for a progress bar
  setTxtProgressBar(prog, i)
}

hist(ants.vec,xlim=c(t.val,max(ants.vec)))
abline(v=t.val,col=2,lwd=2)
quantile(ants.vec,c(.025,.975),na.rm=T)


####Another example with the fish parasites data set
fishPar <- read.csv("fishParasites.csv")
t.val <- U.calc(fishPar)
n.sim <- 5000
fish.vec<- vector()

prog <- txtProgressBar(min=0, max=n.sim, char="*", style=3)
for(i in 1:n.sim){
  tmp <-try(U.calc(IT(fishPar)),silent=T)
  if(is.numeric(tmp)){ fish.vec[i] <- tmp }
  ###Code for a progress bar
  setTxtProgressBar(prog, i)
}

hist(fish.vec,xlim=c(t.val,max(fish.vec)))
abline(v=t.val,col=2,lwd=2)
quantile(ants.vec,c(.025,.975),na.rm=T)


####Integration example###########
x <- seq(-4,4,length=1000) 
y <- dnorm(x,0,1)
plot(x,y,type='l',ylim=c(0,.5))
###Set the lower limits of interest
ll <- -1
ul <- 1
######Hit and miss MC integration#####
###First pick domain
x.dom <- seq(-4,4,length=1000)
y.dom <- seq(0,.5,length=1000)
prog <- txtProgressBar(min=0, max=n.sim, char="*", style=3)
n.sim <- 10000
x.r <- y.r <- vector()
under.curve <- rep(0,n.sim)
for(i in 1:n.sim){
  ###Draw random x and y
  x.r[i] <- sample(x.dom,1)
  y.r[i] <- sample(y.dom,1)
  if(y.r[i] < dnorm(x.r[i],0,1)){
    ###Set further constraints on the x axis
    if(x.r[i] >= ll && x.r[i] <= ul){under.curve[i]<-1}}
  setTxtProgressBar(prog, i)
  
}

plot(-100,-100,ylim=range(y.dom),xlim=range(x.dom),col=2,lwd=2)
points(x.r,y.r,col=under.curve+1,pch=19)
lines(x,y,col=2,lwd=4)

f.h <- sum(under.curve)/n.sim
V.samp <- diff(range(x.dom))*diff(range(y.dom))
area <- f.h*V.samp









