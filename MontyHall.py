
'''
Python version of the Monty Hall problem
by EM Hart 2/20/2012
Change scenario by changing the code 
in strat_dictionary (0,1,2)

'''

#from matplotlib import pyplot
import numpy.random as np
import numpy
from scipy import *
from matplotlib import pyplot

#create an array for Wins
pwins = zeros(1000)
#####Create an an array of all possible values
potential = array([1,2,3])
#####Change this to change your strategy
strat_dict = ["Stay","Switch","Random"]
master_strat = strat_dict[2]
for k in range(1000):
    wins = zeros(100)

    for i in range(100):

####Assign a prize value
        prize = np.random_integers(1,3,1)
###Now make a guess
        guess = np.random_integers(1,3,1)
   
    
####Now we need to figure out which of the doors are revealed
        if prize == guess:
            reveals = potential[where(prize!=potential)]
            reveals = reveals[np.random_integers(0,1,1)]
####Here is where I might have used which in R
        if prize != guess:
            reveals = potential[where(prize!=potential)]
            reveals = reveals[where(guess != reveals)]
####This formulation allows me to have Random strategy
        if master_strat == "Random":
            strat = strat_dict[np.random_integers(0,1,1)]
        if master_strat == "Stay":
            strat = "Stay"
        if master_strat == "Switch":
            strat = "Switch"
        
#Now its simple if we just stay       
        if strat == "Stay":
            guess = guess
####Switch is a bit more complicated, this is a very inelegant solution compared to R
        if strat == "Switch":
            switch = concatenate((guess,reveals))
            for j in range(3):
                exc = potential[j] in switch
                if exc==False:
                    guess = potential[j]

        
    
        if guess == prize:
            wins[i]=1

        pwins[k]= wins.sum()
        

pwins = pwins/1000
pyplot.hist(pwins,100)

pyplot.show()

