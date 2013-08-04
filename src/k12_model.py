from ssa import SSA
from math import *
from scipy.stats import norm
dnorm = norm.rvs
from utils import *
from time_indep_params import *

# State vector is as follows: S1, S2, S3, M, P, C2, C3, CM

def generate_model(model_class):
    kinase = 1
    k12 = sample_forward_kinase_independent_rate()
    k23 = sample_forward_kinase_independent_rate()
    k21 = sample_backward_kinase_independent_rate()
    k32 = sample_backward_kinase_independent_rate()
    # make the appropriate rate kinase dependent
    if model_class == 12 or model_class == 23:
        eval_string = ("k%s = sample_forward_kinase_dependent_rate()" %
                       model_class)
    else:
        eval_string = ("k%s = sample_backward_kinase_dependent_rate()" %
                       model_class)
    exec(eval_string)
    k_r2 = sample_k_r2()
    k_r3 = sample_k_r3()
    k_t = 0 # for now

    init_state = (1,0,0,0,0,0,0,0)
    species_names = "S1,S2,S3,M,P,C2,C3,Cm".split(",")

    model_reactions = [((-1, 1, 0, 0, 0, 0, 0, 0), k12),     # S1 -> S2
                        (( 1,-1, 0, 0, 0, 0, 0, 0), k21),     # S2 -> S1
                        (( 0,-1, 1, 0, 0, 0, 0, 0), k23),     # S2 -> S3
                        (( 0, 1,-1, 0, 0, 0, 0, 0), k32),     # S3 -> S2
                        # S2 -> C2 (begin transcription from S2)
                        (( 0,-1, 0, 0, 0, 1, 0, 0), k_r2),
                        # C2 -> S2 + M (finish transcription from s2)
                        (( 0, 1, 0, 1, 0,-1, 0, 0), FAST),
                        # S3 -> C3 (begin transcription from S3)
                        (( 0, 0,-1, 0, 0, 0, 1, 0), k_r3),
                        # C3 -> S3 + M (finish transcription from S3)
                        (( 0, 0, 1, 1, 0, 0, -1, 0), FAST),
                        # M -> Cm (begin translation)
                        (( 0, 0, 0,-1, 0, 0, 0, 1), k_t),
                        # Cm -> S3 + M (finish translation)
                        (( 0, 0, 0, 1, 1, 0, 0,-1), FAST),    
                        (( 0, 0, 0,-1, 0, 0, 0, 0), gamma_m), # M -> {}
                        (( 0, 0, 0, 0,-1, 0, 0, 0), gamma_p), # P -> {}
                  ]
    return SSA(model_reactions,init_state,species_names=species_names)

def model_class_inference_experiment():
    for model_class in [12,21,23,32]:
        exec("k%s_sims = [generate_model(%s) for i in range(1000)]"
             % (model_class,model_class))
        exec("sims = k%s_sims" % model_class)
        for sim in verbose_gen(sims):
            sim.run(120)
