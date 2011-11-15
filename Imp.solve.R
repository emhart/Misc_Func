# 
# 
# Author: Edmund Hart
# imp.solve - A routine to solve an implicit function
# and put it in an x-y coordinate system for easy graphing
#
###############################################################################



###########Returns a pair of xy coordinates for an implicit function###
#### x is the x data
#### y is the y data
#### my.fun is the function to solve, using the form: my.fun(x,y)
#####function works by calculating a length(x) x length(y) grid of values
####and then takes the minimized value for each x given y
####The size of the grid used is determined by the length of the x and y
#### vectors you supply
imp.solve <- function(x,y,my.fun){
	###Suppress warnings
	options("warn"=-1)
	z<-outer(x,y,my.fun)
####Temporary x to hold solutions

	x1 <- vector()
####Real vector to build on
 my.x <- vector()
 my.y <- vector()
	index <- 1:length(x)
	for(i in 1:length(x)){
	####Now it gets tricky because polynomial functions will
	#have multiple solutions that = 0 
	#  By searching over the grid we need to find all those points where cross overs occur 
	
		#positive solutions
		index2 <- index[z[,i]>0]
		####negative solutions
		index3 <- index[z[,i]<0]
		####now we need to find the number of solutions
		#This works by looking at discontinuities in the index
		#I  have to look in both and make sure the solution counts
		#agree
		n.sol.p <- which(diff(index2)>1)
		n.sol.n <- which(diff(index3)>1)
		
		####now check and see if the same number of solutions exists
		if(length(n.sol.p) > length(n.sol.n)){
		x1 <- x[index2[c(n.sol.p,n.sol.p+1)]]
		}
		
		if(length(n.sol.p) < length(n.sol.n)){
			x1 <- x[index2[c(n.sol.n,n.sol.n+1)]]
		}
		
		if(length(n.sol.p) == length(n.sol.n)){
			
			if(min(index2)==1){x1 <- x[max(index2)]}
			if(max(index2)==length(x)){x1 <- x[min(index2)]}
			
			###now just check if there was a solution
			if(length(index2)==length(x) || length(index2)==0){
			x1 <- NA}
		}
		
	####Now we  bind the solutions together
	my.x <- c(my.x,x1)
	my.y <- c(my.y,rep(y[i],length(x1)))
		
	}
	
	###Now clean it up...
	my.y <- my.y[!is.na(my.x)]
	my.x <- my.x[!is.na(my.x)]
	#t.y <- my.y[order(my.x)]
	#t.x <- sort(my.x)
	
	t.x <- my.x[order(my.y)]
	t.y <- sort(my.y)
	####Create a data fram for easy plotting in ggplot
	tdat <- data.frame(cbind(my.x,my.y))
	colnames(tdat) <- c("x","y")
	return(tdat)
}

####Example###

###Define a function######
#Here is a trivial example
my.fun <- function(x,y){y^2+x^2 -1}
####Set a range to plot over#####

y<-seq(-1.5,1.5,length=1000)
x<-seq(-1.5,1.5,length=1000)

#####collect coordinates
to.plot <- imp.solve(x,y,my.fun)
#Plot with ggplot
ggplot(to.plot,aes(x=x,y=y))+geom_point()+xlim(range(x))+ylim(range(y))



####Compare to a typical use without the imp.solve using a similar method
##But produces ugly plots
x<-seq(-1.5,1.5,length=1000)
y<-seq(-1.5,1.5,length=1000)
z<-outer(x,y,my.fun)
contour(x,y,z,level=0)



###############Here is an ecological example#########
###Here you can recreate the surface to determine transitions
###of a modified Ricker equation
###Recreating the figure from Aviles 1999
#AvilŽs, L. 1999. Cooperation and non-linear dynamics:
#An ecological perspective on the evolution of sociality. 
#Evolutionary Ecology Research 1: 459-477
#See the appendix for the derivation of the implicit equation

##########First set up the parameters to vary over########
y<-seq(0.01,.9,length=1000)
x<-seq(-2,8,length=1000)
###################Define the function############
my.fun <- function(x,y){.01*exp((y-x)/y)-y }
to.plot <- imp.solve(x,y,my.fun)

ggplot(to.plot,aes(x=x,y=y))+geom_line()

