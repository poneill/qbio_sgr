"""
This module provies exact stochastic simulation via the Gillespie
algorithm.

Data structures:

A simulation is a list of reactions + an initial state

A reaction is a list of products, a list of reactants and a reaction
rate constant.

Output:

A list of lists of the form [[time, species1, species2,...,speciesN]]

Must optionally incorporate forcing function, i.e. species
concentration is known at all times

Usage:

Reactions:
3A + 2B <-> 1C (k1)
1B + 1C <-> 1D (k2)
1A + 1C <-> 1B (k(t)) <- forcing function for reaction rate

sim = Simulation([[[-3,-2,1,0],k1],...]) #choose(A,3)*choose(B,2)
                                         #should be inferred from mass action
"""

from utils import transpose,inverse_cdf_sample,product,normalize,zipWith,choose
from random import random,expovariate
from matplotlib import pyplot as plt
rexp = expovariate


class SSA(object):
    def __init__(self,reactions,init_state):
        self.stoich_vectors,self.rate_constants = transpose(reactions)
        self.state = init_state
        self.time = 0
        self.history = [(self.time,self.state)]
        
    def propensity(self,stoich_vector,rate_constant):
        choices = [choose(x_j,-v_j)
                   for x_j,v_j
                   in zip(self.state,stoich_vector)
                   if v_j < 0]
        #print choices
        propensity = rate_constant*product(choices)
        #print "state:",self.state,"stoich:",stoich_vector,"rate const:",rate_constant,"prop:",propensity
        return propensity

    def update(self):
        # Compute propensities
        propensities = [self.propensity(v,k) for v,k in 
                           zip(self.stoich_vectors,self.rate_constants)]
        p = sum(propensities) #rate of sum of random variables
        # determine time of next reaction
        dt = rexp(p)
        #optimize next line later
        # determine which reaction
        v = inverse_cdf_sample(self.stoich_vectors,normalize(propensities))
        # update state vector
        self.state = zipWith(lambda x,y:x+y,self.state,v)
        self.time += dt
        self.history.append((self.time,self.state))
        #print self.state

    def run(self,final_time):
        while self.time < final_time:
            self.update()

    def plot_trajectory(self):
        times, trajs = transpose(self.history)
        for traj in transpose(trajs):
            plt.plot(times,traj)
    
print "loaded"
