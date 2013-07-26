from ssa import SSA
from scipy.stats import binom
rbinom = binom.rvs

k1,k2,g1,g2,k3 = 10,10,1,0,1000
gnxp_reactions = [([1,0,0],k1), #mRNA production
                  ([-1,0,0],g1), #mRNA degradation
                  ([-1,0,1],k2), #protein complex formation
                  ([1,1,-1],k3), #protein formation we do things this way
                  #to pack everything into mass action kinetics
                  ([0,-1,0],g2), #protein degradation
              ]
replication_time = 10
division_prob = 0.5 # prob given molecule goes to given daughter cell
init_state = [0,0,0]

def cell_cycle_sim(time):
    exp_time = 0
    sim = SSA(gnxp_reactions,init_state)
    while exp_time < time:
        exp_time += replication_time
        sim.run(exp_time) # run until exp_time
        new_state = [rbinom(n,division_prob) if n > 0 else 0
                     for n in sim.state]
        sim.state = new_state
        print exp_time
    return sim

