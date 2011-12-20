###############Food web plotting function################
### Date: 12/13/2011
### Author: Edmund Hart (edmund.m.hart@gmail.com)
### Description: A function to create grahps of trophic networks using ggplot2
### Plots food webs in a circular graph
### requires ggplot2 
#########################################################
library(ggplot2)

########## support function create.xy
########## returns regularly spaced circular coordinates for the size 
########## of your web

create.xy <- function(po){
  degs <- seq(0,2*pi,by=(2*pi/(po)))
	return(cbind(cos(degs),sin(degs)))
}




############Plotting function##############
# plot.webgg() requires two arguments:
# ARGUMENTS: web - first is a square S x S matrix of 0's and 1's where S is the total
# number of species in the web
# columns are consumers and rows are species that are consumed
# therefore a 1 in row 5 and column 8 means that species 8 eats species 5
#  Use the function t() if your matrices are oppositely arranged
#labels - is a vector of length S containing string values for the color of each consumer
# species link.
#
# RETURNED VALUES:  The function creates a list with two values.  The first is
# a ggplot2 object accessed as object$plot  
# The next is a dataframe of raw x-y coordinates to connect accessed as object$rawdat
plot.webgg <- function(web,labels){
  	xy <- create.xy(dim(web)[1])
		xy <-data.frame(xy)
		cons <- vector()
		my.plot<- ggplot(xy,aes(x=X1,y=X2))+geom_point()+xlab("")+ylab("")
  	###Suppress warnings
    options(warn=-1)
    ####Create output frame
    p.df <- data.frame(matrix(NA,ncol=4,nrow=0))
		
    ###Create label index
    label.in <- vector()
    
    for(i in 1:dim(web)[1]){
			
      #select the consumer links
			cons <- which(web[,i]==1)
      #####Now I need to have an if statement if there are no links
  if(length(cons)>0){
      #Now create repeat the one set of points that is being drawn
      #ordering with odd numbers
			tmp1 <-data.frame(cbind(rep(xy[i,1],length(cons)),rep(xy[i,2],length(cons)),seq(1,2*length(cons),by=2)))
      #Next select the nodes to connect to
      tmp2 <- data.frame(cbind(xy[cons,]))
      #add an index of sequenced even numbers 
      tmp2$X3 <- seq(2,2*length(cons),by=2)
      tmp <- rbind(tmp1,tmp2)
      #now sort them to get the order correct
      tmp <- tmp[order(tmp$X3),]
      #Finally add a grouping variable
      tmp$color <- rep(paste(labels[i],i),dim(tmp)[1])
      p.df <- rbind(p.df,tmp)
     # A tricky thing can be that your colors will be all wrong if 
     # If I don't adjust the label vector, shortening if for columns
     # that are all 0's so instead I build up a label index
    label.in <- c(label.in,i)
      }
	
    
  }
      my.plot<-my.plot+geom_path(data=p.df,aes(x=X1,y=X2,group=color,colour=color))+scale_colour_manual(values=labels[label.in],legend=F)
		my.plot <-   my.plot+opts(axis.ticks = theme_blank(), axis.title.y = theme_blank(), axis.text.y =  theme_blank(), axis.title.x = theme_blank(), axis.text.x =  theme_blank()) 
	return(list(plot=my.plot,rawdat=p.df))
}



##############Sample Data Generation######
# Data inputs:  In this formulation the species that is consumed are the rows
# and the consumers are the columns
# Data example: Here I use a random network example


######Large random matrix
rand.mat <- matrix(rbinom(1600,1,.2),ncol=40,nrow=40)
diag(rand.mat)<- 0
lab <- rep("Black",40)
test<- plot.webgg(rand.mat,lab)
test$plot

####Small random matrix with colored links
rand.mat <- matrix(rbinom(100,1,.3),ncol=10,nrow=10)
diag(rand.mat)<- 0
lab <- c(rep("blue",3),rep("red",7))
test<- plot.webgg(rand.mat,lab)
test$plot
##############Note that overlapping colors
# can cause other colors (in this example it would be purple)
# this will normally not happen unless you have reciprocal links
# if anyone really wants that fixed I could figure it out 
# lastly here is an example using a real food web from 
# Bascompte et al. 10.1073/pnas.0501562102 at http://www.pnas.org/cgi/content/full/0501562102/DC1
# found at http://knb.ecoinformatics.org/knb/metacat?action=read&qformat=nceas&docid=bowdish.272 
#
basc.web <- as.matrix(read.csv("basc_web.csv",header=F))
####Note that I cleaned it up by removing the row and column
#labels in the original files and saved it as a CSV

#now replace values greater than 0 but less than 1 with 1's
for(i in 1:249){basc.web[which(basc.web[,i]>0),i]<-1 }
labs <- rep("Black",249)
basc.plot <- plot.webgg(basc.web,labs)
basc.plot$plot

##########Finally an example generated with a niche-foodweb simulation
# from Williamns and Martinez: http://www.nature.com/nature/journal/v404/n6774/abs/404180a0.html
#  I have a brief undocumented function below to generate the web, it 
# requires S (# of species) and C (connectivity)

niche.model <- function(S,C){
  new.mat <- matrix(0,nrow=S,ncol=S)
	ci <- vector()
	niche <- runif(S,0,1)
	r <- rbeta(S,1,((1/(2*C))-1)) * niche
	for(i in 1:S){ci[i]<-runif(1,r[i]/2,niche[i])}
	
	#now set the smallest species niche value to have an n of 0
	r[which(niche==min(niche))] <- .00000001
  for(i in 1:S){
		for(j in 1:S){
		   if(niche[j] > (ci[i]-(.5*r[i])) && niche[j]< (ci[i]+.5*r[i])){new.mat[j,i]<-1}
			}
    }
new.mat <- new.mat[,order(apply(new.mat,2,sum))]
return(new.mat)
}
###Plot
niche.web <- niche.model(20,.25)
labs <- rep("Black",20)
niche.plot <- plot.webgg(niche.web,labs)
niche.plot$plot


  

  