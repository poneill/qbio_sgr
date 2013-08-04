from time_dependent_ssa import SSA
from math import *
from scipy.stats import norm
dnorm = norm.rvs
from utils import *
from collections import defaultdict
import mpmath

gamma_m_rate = log(2)/15
gamma_p_rate = 0 #log(2)/20 #???
sample_kr2_rate = lambda :max(0,dnorm(0.5,0.1))
sample_kr3_rate = lambda :max(0,dnorm(2,0.4))
sample_forward_kinase_independent_rate = lambda : max(0,dnorm(1,(0.2))) # mean, sd
sample_backward_kinase_independent_rate = lambda : max(0,dnorm(2,(0.4))) # mean, sd
sample_forward_alpha = lambda :dnorm(-1,0.2)
sample_backward_alpha = lambda :dnorm(1,0.2)
sample_forward_beta = lambda :dnorm(2,0.5)
sample_backward_beta = lambda :dnorm(-2,0.5)
k_t_rate = 0 # for now
FAST = 10**6

def kinase(t):
    return 0 if t < 120 else 10

def kinase2(t):
    if t < 240:
        if t < 120:
            return 0
        else: 
            return 0.5
    else:
        if t < 360:
            return 1
        else:
            return 2

def f_k12_indep(k):
    #k = sample_forward_kinase_independent_rate()
    return lambda (s1,s2,s3,m),t:k*s1

def f_k12_dep(alpha,beta,kinase):
    # alpha = sample_forward_alpha()
    # beta = sample_forward_beta()
    return lambda (s1,s2,s3,m),t:max(0,alpha + beta*kinase(t))*s1

def f_k21_indep(k):
    return lambda (s1,s2,s3,m),t:k*s2

def f_k21_dep(alpha,beta,kinase):
    return lambda (s1,s2,s3,m),t:max(0,alpha + beta*kinase(t))*s2

def f_k23_indep(k):
    return lambda (s1,s2,s3,m),t:k*s2

def f_k23_dep(alpha,beta,kinase):
    return lambda (s1,s2,s3,m),t:max(0,alpha + beta*kinase(t))*s2

def f_k32_indep(k):
    return lambda (s1,s2,s3,m),t:k*s3

def f_k32_dep(alpha,beta,kinase):
    return lambda (s1,s2,s3,m),t:max(0,alpha + beta*kinase(t))*s3

def f_kr2(kr2):
    return lambda (s1,s2,s3,m),t:kr2*s2

def f_kr3(kr3):
    return lambda (s1,s2,s3,m),t:kr3*s3

def f_gamma_m():
    return lambda (s1,s2,s3,m),t:m*gamma_m_rate

# State vector is as follows: S1, S2, S3, M, P, C2, C3, CM
model_classes = [12,21,23,32]
forwards = [12,23]
backwards = [21,32]

def generate_model(model_class,kinase,init_run_time=0):
    k12_rate = sample_forward_kinase_independent_rate()
    k12 = f_k12_indep(k12_rate)
    k21_rate = sample_backward_kinase_independent_rate()
    k21 = f_k21_indep(k21_rate)
    k23_rate = sample_forward_kinase_independent_rate()
    k23 = f_k23_indep(k23_rate)
    k32_rate = sample_backward_kinase_independent_rate()
    k32 = f_k32_indep(k32_rate)
    for mc in [12,21,23,32]:
        if mc in forwards:
            alpha_rate = sample_forward_alpha()
            beta_rate = sample_forward_beta()
        elif mc in backwards:
            alpha_rate = sample_backward_alpha()
            beta_rate = sample_backward_beta()
        else:
            raise Exception("invalid model class")
        if mc == model_class:
            fmt_str = (mc,mc, alpha_rate,beta_rate)
            exec("k{0} = f_k{1}_dep({2},{3},kinase)".format(*fmt_str))
    kr2_rate = sample_kr2_rate()
    kr3_rate = sample_kr3_rate()
    kr2 = f_kr2(kr2_rate)
    kr3 = f_kr3(kr3_rate)
    gamma_m = f_gamma_m()

    init_state = (1,0,0,0)
    species_names = "S1 S2 S3 M P C2 C3 Cm".split()
    
    model_reactions = [((-1, 1, 0, 0), k12),     # S1 -> S2
                        (( 1,-1, 0, 0), k21),     # S2 -> S1
                        (( 0,-1, 1, 0), k23),     # S2 -> S3
                        (( 0, 1,-1, 0), k32),     # S3 -> S2
                        # S2 -> C2 (transcribe from S2)
                        (( 0,0, 0, 1), kr2),
                        # S3 -> C3 (transcribe from S3)
                        (( 0, 0,0, 1), kr3),
                        # M -> Cm (translate)
                        (( 0, 0, 0,-1), gamma_m), # M -> {}
                       ((0,0,0,0),lambda (s1,s2,s3,m),t:1*s1)
                  ]
    sim = SSA(model_reactions,init_state,species_names=species_names)
    sim.model_class = model_class
    for param in "k12 k21 k23 k32 kr2 kr3 gamma_m alpha beta".split():
        exec_string = "sim.{0} = {1}_rate".format(param,param)
        exec(exec_string)
    sim.parameters = (sim.k12,
                      sim.k21,
                      sim.k23,
                      sim.k32,
                      sim.kr2,
                      sim.kr3,
                      sim.alpha,
                      sim.beta)
    sim.run(init_run_time)
    return sim

def log_likelihood(p_sample,q_sample):
    """Compute ll of ps, given qs"""
    pseudo_qs = list(q_sample) + range(100)
    q_n = float(len(pseudo_qs))
    qs = {k:v/q_n for k,v in Counter(pseudo_qs).items()} # pseudocount
    return sum(log(qs[p]) for p in p_sample)

def model_class_inference_experiment():
    reference_sims = {}
    num_sims = 1000
    for model_class in model_classes:
        print "generating samples for:",model_class
        sims = [generate_model(model_class,kinase,240) for i in range(num_sims)]
        reference_sims[model_class] = sims
    correct_guesses = 0
    trials = 100
    for trial in range(trials):
        mc = random.choice(model_classes)
        print "beginning trial:",trial
        #trial_sims = [generate_model(mc,kinase,240) for i in range(num_sims)]
        trial_sims = [generate_model(mc,kinase,120) for i in range(num_sims)]
        trial_mrnas_120 = [sim.state_at(120)[3] for sim in trial_sims]
        trial_mrnas_240 = [sim.state_at(240)[3] for sim in trial_sims]
        lls = {}
        for model_class in model_classes:
            print "comparing", mc," to: ",model_class
            ref_sims = reference_sims[model_class]
            ref_mrnas_120 = [sim.state_at(120)[3] for sim in ref_sims]
            ref_mrnas_240 = [sim.state_at(240)[3] for sim in ref_sims]
            ll_120 = log_likelihood(ref_mrnas_120, trial_mrnas_120)
            ll_240 = log_likelihood(ref_mrnas_240, trial_mrnas_240)
            print "partial lls:",ll_120,ll_240
            ll = ll_120 + ll_240
            print "total ll:",ll
            lls[model_class] = ll
        best_class = max(lls.keys(),key=lambda k:lls[k])
        if best_class == mc:
            print "Guessed correctly:",lls[best_class], "vs:",[lls[c] for c in model_classes if not c == best_class]
            correct_guesses += 1 
        else:
            print "Guessed incorrectly:",lls[best_class],"vs:",lls[mc]
    print "guessed correctly:",correct_guesses/float(trials)
print "loaded"

def likelihood_experiment(mrnas120,mrnas240):
    mc = 21 # by previous experiment
    ref_sims = [generate_model(mc,kinase) for i in range(1000)]
    all_mrnas120 = []
    all_mrnas240 = []
    for j,ref_sim in enumerate(ref_sims):
        print "ref_sim:",j
        mrnas120 = []
        mrnas240 = []
        for i in range(100):
            print i
            ref_sim.history = []
            ref_sim.time = 0
            ref_sim.run(240)
            mrnas120.append(ref_sim.state_at(120)[3])
            mrnas240.append(ref_sim.state_at(240)[3])
        all_mrnas120.append(mrnas120)
        all_mrnas240.append(mrnas240)
    return all_mrnas120,all_mrnas240,ref_sims

param_estimate = [1.0199423579837563, 2.621820189733821, 1.0811274643284017, 2.621820189733821, 0.5214392227617346, 1.6005073371928782, 0.9140125458538542, -2.3908331183609306]

def likelihood_experiment2(timepoints,data,kinase,samples,runs_per_sample):
    """Given timepoints,data at said time points, and kinase function,
    return a list of parameter values with associated likelihoods"""
    mc = 21
    ref_sims = [generate_model(mc,kinase2) for i in range(samples)]
    output = []
    for j,ref_sim in enumerate(ref_sims):
        print "ref_sim:",j
        mrna_dict = defaultdict(list)
        for i in range(runs_per_sample):
            print i
            ref_sim.history = []
            ref_sim.time = 0
            ref_sim.run(480)
            for tp in timepoints:
                mrna_dict[tp].append(ref_sim.state_at(tp)[3])
        ll = sum(log_likelihood(mrna_dict[tp],data[i])
                 for i,tp in enumerate(timepoints))
        print (ref_sim.parameters,ll)
        output.append((ref_sim.parameters,ll))
    return output

def samples_per_trial_experiment():
    """How many samples are required to reliably estimate the ll?"""
    sim = generate_model(21,kinase2)
    lls = []
    for j in range(10):
        mrna_dict = defaultdict(list)
        for i in range(runs_per_sample):
            print i
            ref_sim.history = []
            ref_sim.time = 0
            ref_sim.run(480)
            for tp in timepoints:
                mrna_dict[tp].append(ref_sim.state_at(tp)[3])
            ll = sum(log_likelihood(mrna_dict[tp],data[i])
                     for i,tp in enumerate(timepoints))
        lls.append(ll)
        print "ll:",ll
    return lls

def mean_from_output(output):
    """Given output of the form [params,ll], return mean parameter
    vector"""
    param_vectors,lls = transpose(output)
    z = sum(mpmath.exp(ll) for ll in lls)
    return map(sum,transpose([[(mpmath.exp(ll)/z)*p for p in param_vector]
                               for param_vector,ll in output]))


