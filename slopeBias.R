#####Generate noisy data
###Create design matrix
X<- cbind(rep(1,100),rnorm(100,0,4),rnorm(100,0,4),rnorm(100,0,4),rnorm(100,0,4))
###Set effects sizes
effs <- c(1,1,1,1,1)
Y<- X%*%effs

####Full model shows perfect correlation
full.mod <- lm(Y~X)
###Choose random covariate, 2 in this case
m1 <- lm(Y~X[,3])
my.y <- Y
my.x <- X[,3]


####Now create some output matrices to store the data####
##Store the r^2 of the original model, the r^2 (correlation) of the between the true 
## x and the measured x, and the sd of the newly fit model
sdvars <- seq(.001, 20,length=1000)
xerr.output <- matrix(0,ncol=5,nrow=1000)
colnames(xerr.output) <- c("alpha","beta","rsq","cor","sd")

for(i in 1:1000){
  ####create new variables with measurement error
  xstar <- my.x + rnorm(100,0,sd=sdvars[i])
  tmp.m1 <- lm(my.y~xstar)
  cormod <- lm(xstar~my.x)  
  ####Collect output
  xerr.output[i,] <- c(coef(tmp.m1),summary(tmp.m1)$r.sq,summary(cormod)$r.sq,sd(resid(tmp.m1)))
  
}

yerr.output <- matrix(0,ncol=5,nrow=1000)
colnames(yerr.output) <- c("alpha","beta","rsq","cor","sd")

for(i in 1:1000){
  ####create new variables with measurement error
  ystar <- my.y + rnorm(100,0,sd=sdvars[i])
  tmp.m1 <- lm(ystar~my.x)
  cormod <- lm(ystar~my.y)  
  ####Collect output
  yerr.output[i,] <- c(coef(tmp.m1),summary(tmp.m1)$r.sq,summary(cormod)$r.sq,sd(resid(tmp.m1)))
  
}

####Set up data frames for plotting

yerr.output <- data.frame(yerr.output)
xerr.output <- data.frame(xerr.output)

####Code to make first blog figure#####
plot(my.y~my.x, ylab="Stupid Response",xlab="Silly Question")
for(i in 1:100){abline(xerr.output[i,1:2],lwd=.5,col=2)}
for(i in 1:100){abline(yerr.output[i,1:2],lwd=.5,col=4)}
abline(m1,lty=2,lwd=3)

comb.out <- rbind(yerr.output,xerr.output)
comb.out$errt <- c(rep("y",1000),rep("x",1000))
###Code for slope figures
library(ggplot2)
#####
ggplot(comb.out,aes(x=sqrt(cor),y=beta,colour=errt))+geom_point()+scale_colour_discrete(name="Variable \n with error")+xlab("Correlation with true value")+ylab("Slope")+geom_abline(intercept=1.1638,slope=0)


###Code for R^2 figure
ggplot(comb.out,aes(x=sqrt(cor),y=rsq,colour=errt))+geom_point()+scale_colour_discrete(name="Variable \n with error")+xlab("Correlation with true value")+ylab("R-Sqared")+geom_abline(intercept=0.2617,slope=0)




