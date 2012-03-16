

#this is simply the log 5 rule function to calculate the point estimate of the Pomeroy rankings
l5rule <- function(ratings){
a <- ratings[1]
b <- ratings[2]
wpct <- (a - (a*b))/ ((a + b) - (2*a*b))
return(wpct)}


#this function is used in my prior calculation function
#It will take a probability of winning and a confidence ranking and convert it into a beta prior
#parameters
find.betap <- function(p){
#you can adjust the 2 to change how tight the distribution should get.  So by adjusting this parameter
#you can add more weight to the priors
a<- 3
bp <- a / p
b <- bp-a
return(c(a,b))}


#prior.db's should be in two columns C1 is the odds of winning a region (or the champhionship), C2 the
#cofidence and C3 the team index

#######Calculates beta distribution values based on power rankings
###### It assumes a 0 distribution with an calculated sd and then calculates
###### then integrates the point spread with a normal distribution
######  I add that difference to .5 to get a probability of victory
####pvec - a vector of the two power rankings
#### rankvals - a vector of the mean and variance of the power rankings
#### flat - if True just returns a 1,1
calc.priors <- function(pvec,rankvals,flat=F){
#if this is flagged, it will automatically return flat priors
if(flat==T){return(c(1,1))}
#This asks what is the probablity that team.index[1] beats team.index[2]

p1 <- pvec[1]
p2 <- pvec[2]
p.win <- -1*diff(pnorm(c(p1,p2),rankvals[1],rankvals[2]))

#####Conf is a measure of your confidence in the prior.  Larger values
##### mean more confidence and a tighter beta distribution
conf <-  2
p.win <- .5+p.win

if(p.win > 1){p.win <- .99}
if(p.win < 0){p.win <- .01}

#finally I'll solve for the beta parameters
return(find.betap(p.win))}

#a quick function to reassign position numbers
posassign <- function(vec){
out <- vector()
new.length <- length(vec)/2
for(x in 1:new.length){
out <- c(out,rep(x,2))}
return(out)}


                
                

#Setting up my data
setwd("H:/NCAAModel/")
ncaa.dat <- read.csv("ncaa2012.csv")



#The total number of simulations of the entire tournament
n.sim <- 10000
#The total number of games to simulate to make the fake data based of the probability derived by the log5 rule
#I chose 34 because thats the number of games in a regular season, but we can change this in sensitivity tests
n.games <- 20

######Set the priors####
###I have two sources of priors information####
####Priors 1 http://www.usatoday.com/sports/sagarin/bkt1011.htm
#### Priors 2 http://sonnymoorepowerratings.com/m-basket.htm
### I will calculate each beta value and then take the mean values
### to create average prior info
priors1 <- c(mean(ncaa.dat$powerscore1),sd(ncaa.dat$powerscore1))
priors2 <- c(mean(ncaa.dat$powerscore2),sd(ncaa.dat$powerscore2))


#This function will simulate a division champion
region.sim <- function(bracket.dat,priors) {

   #an output matrix that holds the wins as 1 and the losses as 0's
       #if a team wins the first round it gets a 1 in column 2, etc...
                
out.mat <- cbind(bracket.dat$index,matrix(0,nrow=16,ncol=4))

#First I set up the loop to cover all the rounds of each division
#in this case 4 rounds to establish a champion

for(i in 1:4){
      #next I can get the number of games from the position code in the file
      div.games <- max(unique(bracket.dat$position)) #to start this is 8
      

      for(j in 1:div.games){
             
                 #first I need to figure out the index of the two teams playing
                 
                #Now I simulate 34 games between the first team
                 opponents.ratings <- bracket.dat[bracket.dat[,7]==j,2]
                 index <- bracket.dat[bracket.dat[,7]==j,8]
                #first I get my point estimate of theta
                 theta <- l5rule(opponents.ratings)
                 #now I simulate the odds of team A beating team B
                 game.sim <- rbinom(n.games,1,theta)
                 
                #now I randomly draw from the posterior distribution of this game given my simulated data
                #and my priors (set to 1,1)
                
                #I calculate the number of wins and losses
                s <- length(which(game.sim==1))
                f <- length(which(game.sim==0))
                #Then I add my prior information
                # I have two sources of prior info detailed above
                # I will calculate each one and then take the average
                p1 <-  bracket.dat[bracket.dat[,7]==j,4]
                beta1 <-  calc.priors(p1,priors1)
                p2 <-  bracket.dat[bracket.dat[,7]==j,5]

                beta2 <- calc.priors(p2,priors2)
                 a <- mean(c(beta1[1],beta2[1]))
                 b <- mean(c(beta1[2],beta2[2]))

                
                #now I'll take a single random draw from  my posterior beta distribution
                
                outcome <- rbinom(1,1,rbeta(1,a+s,b+f))
                
                #record the outcome in the output matrix
                if(outcome==1){out.mat[which(out.mat[,1]==index[1]),i+1]<-1}
                if(outcome==0){out.mat[which(out.mat[,1]==index[2]),i+1]<-1}
                
            
      }
      
       #now I need to readjust the bracket positions  for the first three rounds
          if(i < 4){
             #subset out the winners
             winners <- out.mat[out.mat[,i+1]>0,1]
             bracket.dat <- bracket.dat[bracket.dat[,8]%in%winners,]
             #reorder the bracket by position and then re assign position numbers
             bracket.dat <- bracket.dat[order(bracket.dat$position),]
             bracket.dat$position <- posassign(bracket.dat$position)
             }
      
}

return(out.mat)}



champion.sim <- function(bracket.dat,priors)  {

#this function works the same as a region but the only difference are the priors 
# and the size of the output matrix


   #an output matrix that holds the wins as 1 and the losses as 0's
       #if a team wins the first round it gets a 1 in column 2, etc...
                
out.mat <- cbind(bracket.dat$index,matrix(0,nrow=4,ncol=2))

#First I set up the loop to cover all the rounds of each division
#in this case 4 rounds to establish a champion

for(i in 1:2){
      #next I can get the number of games from the position code in the file
      div.games <- max(unique(bracket.dat$position)) #to start this is 8
      

      for(j in 1:div.games){
             
                 #first I need to figure out the index of the two teams playing
                 
                #Now I simulate 34 games between the first team
                 opponents.ratings <- bracket.dat[bracket.dat[,7]==j,2]
                 index <- bracket.dat[bracket.dat[,7]==j,8]
                #first I get my point estimate of theta
                 theta <- l5rule(opponents.ratings)
                 #now I simulate the odds of team A beating team B
                 game.sim <- rbinom(n.games,1,theta)
                 
                #now I randomly draw from the posterior distribution of this game given my simulated data
                #and my priors (set to 1,1)
                
                #I calculate the number of wins and losses
                s <- length(which(game.sim==1))
                f <- length(which(game.sim==0))
                #Then I add my prior information
              p1 <-  bracket.dat[bracket.dat[,7]==j,4]
                beta1 <-  calc.priors(p1,priors1)
                p2 <-  bracket.dat[bracket.dat[,7]==j,5]

                beta2 <- calc.priors(p2,priors2)
                 a <- mean(c(beta1[1],beta2[1]))
                 b <- mean(c(beta1[2],beta2[2]))
                #now I'll take a single random draw from  my posterior beta distribution
                
                outcome <- rbinom(1,1,rbeta(1,a+s,b+f))
                
                #record the outcome in the output matrix
                if(outcome==1){out.mat[which(out.mat[,1]==index[1]),i+1]<-1}
                if(outcome==0){out.mat[which(out.mat[,1]==index[2]),i+1]<-1}
                
            
      }
      
       #now I need to readjust the bracket positions  for the first three rounds
          if(i < 2){
             #subset out the winners
             winners <- out.mat[out.mat[,i+1]>0,1]
             bracket.dat <- bracket.dat[bracket.dat[,8]%in%winners,]
             #reorder the bracket by position and then re assign position numbers
             bracket.dat <- bracket.dat[order(bracket.dat$position),]
             bracket.dat$position <- posassign(bracket.dat$position)
             }
      
}

return(out.mat)}



#this is the master data array that holds the index value of a team and its record

master.out<- array(0,c(n.sim,64,7))
#this is the array to hold all the region tournament results
region.out <- array(0,c(4,n.sim,16,5))
#I set my priors here

for(x in 1:n.sim){

      for(i in 1:4){
            #first we simulate each region
            #we subset the original data set
            bracket.dat <- subset(ncaa.dat,ncaa.dat$region==i)
            region.out[i,x,,] <- region.sim(bracket.dat,regprior)
      }
   
 #next I put the regional data into the master data set and create the championship bracket
 
 master.out[x,1:16,1:5] <- region.out[1,x,,]
  master.out[x,17:32,1:5] <- region.out[2,x,,]
   master.out[x,33:48,1:5] <- region.out[3,x,,]
    master.out[x,49:64,1:5] <- region.out[4,x,,]
    
#now I will create the championship bracket
champ.index <- master.out[x,master.out[x,,5]==1,1]
champ.dat <- ncaa.dat[ncaa.dat[,8]%in%champ.index,]
champ.dat <- champ.dat[order(champ.dat$region),]
champ.dat$position <- posassign(champ.dat$position)
champ.out <- champion.sim(champ.dat,kapchampprior)
cat(x,"/n")
for(i in 1:4){
master.out[x,which(champ.out[i,1] == master.out[x,,1]),6:7] <- champ.out[i,2:3]}


}


#now I can processes the data and create a table with it

summed.mat <- matrix(0,nrow=64,ncol=6)
for(x in 1:64){

  for(i in 1:6){
  
          summed.mat[x,i] <- sum(master.out[,x,i+1])
  }
  }
  
  p.mat <- (summed.mat/10000)*100

final.output <- cbind(ncaa.dat[1],p.mat)

flat.output <- final.output
combined.output <- final.output








