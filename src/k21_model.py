from ssa import SSA
from math import *
from scipy.stats import norm
dnorm = norm.rvs

gamma_m = log(2)/15
gamma_p = log(2)/15 #???
sample_k_r2 = lambda :dnorm(0.5,0.1)
sample_k_r3 = lambda :dnorm(2,0.4)
sample_forward_kinase_independent_rate = lambda : dnorm(1,sqrt(0.2)) # mean, sd
sample_backward_kinase_independent_rate = lambda : dnorm(2,sqrt(0.4)) # mean, sd
sample_forward_alpha = lambda :dnorm(-1,sqrt(0.2))
sample_backward_alpha = lambda :dnorm(1,sqrt(0.2))
sample_forward_beta = lambda :dnorm(2,sqrt(0.5))
sample_backward_beta = lambda :dnorm(-2,sqrt(0.5))
k_t = 10
FAST = 10**6
# State vector is as follows: S1, S2, S3, M, P, C2, C3, CM

kinase = 1
alpha = sample_forward_alpha()
beta = sample_forward_beta()
k12 = sample_forward_kinase_independent_rate()
k21 = max(0,alpha + beta * kinase)
k23 = sample_forward_kinase_independent_rate()
k32 = sample_backward_kinase_independent_rate()
k_r2 = sample_k_r2()
k_r3 = sample_k_r3()


init_state = (1,0,0,0,0,0,0,0)
species_names = "S1,S2,S3,M,P,C2,C3,Cm".split(",")

model1_reactions = [((-1, 1, 0, 0, 0, 0, 0, 0), k12),  # S1 -> S2
                    (( 1,-1, 0, 0, 0, 0, 0, 0), k21),  # S2 -> S1
                    (( 0,-1, 1, 0, 0, 0, 0, 0), k23),  # S2 -> S3
                    (( 0, 1,-1, 0, 0, 0, 0, 0), k32),  # S3 -> S2
                    (( 0,-1, 0, 0, 0, 1, 0, 0), k_r2),  # S2 -> C2 (begin transcription from S2)
                    (( 0, 1, 0, 1, 0,-1, 0, 0), FAST),  # C2 -> S2 + M (finish transcription from s2)
                    (( 0, 0,-1, 0, 0, 0, 1, 0), k_r3),  # S3 -> C3 (begin transcription from S3)
                    (( 0, 0, 1, 1, 0, 0, -1, 0), FAST),  # C3 -> S3 + M (finish transcription from s2)
                    (( 0, 0, 0,-1, 0, 0, 0, 1), k_t),  # M -> Cm (begin transcription from S3)\
                    (( 0, 0, 0, 1, 1, 0, 0,-1), FAST),  # Cm -> S3 + M (finish transcription from s2)
                    (( 0, 0, 0,-1, 0, 0, 0, 0), gamma_m),  # M -> {}
                    (( 0, 0, 0, 0,-1, 0, 0, 0), gamma_p),  # P -> {}
                  ]

k21_sims = [SSA(model1_reactions,init_state,species_names=species_names) for i in range(10)]
print "loaded"
for i,sim in enumerate(k21_sims):
    print i
    try:
        sim.run(60)
    except ZeroDivisionError:
        continue


