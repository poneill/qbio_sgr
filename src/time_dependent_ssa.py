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
    def __init__(self,reactions,init_state,species_names):
        self.stoich_vectors,self.propensities = transpose(reactions)
        self.state = init_state
        self.time = 0
        self.history = [(self.time,self.state)]
        self.verbose = False
        self.species_names = species_names
        self.reactions_performed = 0
        self.finished_run = False

    def update(self,final_time):
        # Compute propensities
        propensities = [propensity(self.state,self.time)
                        for propensity in self.propensities]
        self.logging(propensities)
        p = sum(propensities) #rate of sum of random variables
        # determine time of next reaction
        if p == 0:
            raise Exception("No possible reactions")
        dt = rexp(p)
        #optimize next line later
        # determine which reaction
        new_time = self.time + dt
        if new_time < final_time:
            v = inverse_cdf_sample(self.stoich_vectors,normalize(propensities))
    # update state vector
            self.state = zipWith(lambda x,y:x+y,self.state,v)
            self.time = new_time
            self.logging(str(self.time) + " " + str(self.state))
        
            self.history.append((self.time,self.state))
        #print self.state
            self.reactions_performed += 1
        else:
            self.finished_run = True
        
    def run(self,final_time):
        self.finished_run = False
        while not self.finished_run:
            self.update(final_time) # sets self.finished_run to true
        # except Exception as e:
        #     if e.message == "No possible reactions":
        #         return
        #     else:
        #         raise e

    def state_at(self,desired_time):
        key = lambda(time,state):time if time < desired_time else 0
        best_time,best_state = max(self.history,
                                   key=key)
        return best_state

    def plot_trajectory(self):
        times, states = transpose(self.history)
        trajs = transpose(states)
        for traj,name in zip(trajs,self.species_names):
            plt.plot(times,traj,label=name)

    def logging(self,x):
        if self.verbose:
            print x

        
print "loaded"

