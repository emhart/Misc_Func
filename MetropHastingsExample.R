
data <- rnorm(20, 30,5)

########## This set of functions creates a metropolis-hastings
########## sampler for a vector of normally distributed data
########## It returns a data file with three elements, the raw samples
########## of mu, sigma, and then the acceptance rate

sample.sigmu <-function(x,dx,dsd,sig,mu){
sdstar <- runif(1,sig-dsd,sig+dsd) 
xstar = runif(1,mu-dx,mu+dx)
while(sdstar < 0){
  	sdstar <- runif(1,sig-dsd,sig+dsd)
	}

logalpha = sum(dnorm(x,mean=xstar,sd=sdstar,log=T) - dnorm(x,mu,sig,log=T))- log(sdstar/sig)
logu = log(runif(1,0,1))
acc = (logu < logalpha)
my.mu = acc*xstar + (1-acc)*mu
my.sig = acc*sdstar + (1-acc)*sig
return(list(mu = my.mu,sig=my.sig, acc=acc))
}

sample.mh<-function(n=10000,x,mu=3,sig=4,dx=2,dsd=2)
{

accab  <- 0
mus<-rep(NaN,n)
sigs <- rep(NaN,n)
for(i in 1:n){
z = sample.sigmu(x,dx,dsd,sig,mu)

mus[i] = mu = z$mu
sigs[i] = sig = z$sig

accab = accab + z$acc


}
return(invisible(list(mean=mus, sig=sigs, accab=accab/n)))
}

mh.exam <- sample.mh(x=data)
######Now compare to our actual data
rbind(c(mean(data),mean(mh.exam$m)),c(sd(data),mean(mh.exam$sig)))
####We can plot our functions to see if they were sticky or well mixed
par(mfrow=c(1,2))
plot(mh.exam$m)
plot(mh.exam$s)




#using similar code as above, setting up a simple linear regression.
x <- c(44.61646,39.24717, 43.77985, 40.40638, 43.59542, 38.44058, 37.92276, 41.35069, 43.58636, 37.81024)

y <- c(48.49584, 44.69342, 48.62739, 44.46996, 48.37165, 40.59726, 43.14914, 48.73189, 54.60816, 41.23829)

####Function to center the data
fn <- function(x,a=0,b=1){
	
	centered<- a+b*(x-mean(x))
	return(centered)
	}

sample.b <-function(x,y,alpha,beta,sig,db){

	bstar <- abs(runif(1,beta-db,beta+db))


logalpha = sum(dnorm(y,mean=fn(x,alpha,bstar),sd=sig,log=T) - dnorm(y,mean=fn(x,alpha,beta),sd=sig,log=T))

logu = log(runif(1,0,1))
acc = (logu < logalpha)


my.beta = acc*bstar + (1-acc)*beta

return(list(beta=my.beta, accb=acc))
}
sample.sig <-function(x,y,alpha,beta,sig,dsd){

sdstar <- abs(runif(1,sig-dsd,sig+dsd))




logalpha = sum(dnorm(y,mean=fn(x,alpha,beta),sd=sdstar,log=T) - dnorm(y,fn(x,alpha,beta),sd=sig,log=T)) - log(sdstar/sig)
logu = log(runif(1,0,1))
accs = (logu < logalpha)

my.sig = accs*sdstar + (1-accs)*sig

return(list(sig = my.sig,accs = accs))}

sample.a <-function(x,y,alpha,beta,sig,da){
	astar = abs(runif(1,alpha-da,alpha+da))
	logalpha = sum(dnorm(y,mean=fn(x,astar,beta),sd=sig,log=T) - dnorm(y,mean=fn(x,alpha,beta),sd=sig,log=T))

logu = log(runif(1,0,1))
acc = (logu < logalpha)

my.alpha = acc*astar + (1-acc)*alpha
	return(list(alpha = my.alpha,acca = acc))}
	

#############
### Here we put it all together
#######



sample.mh.linreg<-function(n=1000,x,y,alpha=-2,beta=3,sig=5,dsd=2,da=2,db=2)
{

acca  <- 0
accb <- 0
accs <- 0
alphas<-rep(NaN,n)
betas<-rep(NaN,n)
sigs <- rep(NaN,n)
for(i in 1:n){
z = sample.a(x,y,alpha,beta,sig,da)
#cat("alpha and beta to be passed to q ",z$alpha," ",z$beta,"\n")
q <- sample.sig(x,y,alpha,beta,sig,dsd)
v <- sample.b(x,y,alpha,beta,sig,db)

betas[i] = beta = v$beta
alpha = z$alpha
alphas[i] <- mean(y)-mean(x)*beta
sigs[i] = sig = q$sig
acca = acca + z$acca
accs = accs + q$accs
accb = accb + v$accb

}
return(invisible(list(sig = sigs,alpha = alphas,beta =betas, accb=accb/n,acca=acca/n,accs=accs/n)))
}

my.n <- 5000
lreg.test <- sample.mh.linreg(n=my.n,x=x,y=y)

par(mfrow=c(3,1))
#####I'm allowing for a 1000 step burnin'
plot(lreg.test$sig[1000:my.n])
plot(lreg.test$alpha[1000:my.n])
plot(lreg.test$beta[1000:my.n])


####Compare to normal linear regression estimates
my.est <- lm(y~x)
my.table <- rbind(c(mean(lreg.test$alpha[1000:my.n]),my.est$coef[1]),c(mean(lreg.test$beta[1000:my.n]),my.est$coef[2]),c(mean(lreg.test$sig[1000:my.n]),sd(resid(my.est))))
colnames(my.table)<-c("MH Est","LS Est")
rownames(my.table)<-c("Intercep","Slope","Std Dev")
my.table


