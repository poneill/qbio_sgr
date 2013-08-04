kinase = 10
gamma_m = log(2)/15
gamma_p = 0 #log(2)/20 #???
sample_k_r2 = lambda :max(0,dnorm(0.5,0.1))
sample_k_r3 = lambda :max(0,dnorm(2,0.4))
sample_forward_kinase_independent_rate = lambda : max(0,dnorm(1,(0.2))) # mean, sd
sample_backward_kinase_independent_rate = lambda : max(0,dnorm(2,(0.4))) # mean, sd
sample_forward_alpha = lambda :dnorm(-1,0.2)
sample_backward_alpha = lambda :dnorm(1,0.2)
sample_forward_beta = lambda :dnorm(2,0.5)
sample_backward_beta = lambda :dnorm(-2,0.5)
sample_forward_kinase_dependent_rate = lambda: max(0,sample_forward_alpha() +
                                                   sample_forward_beta() * kinase)
sample_backward_kinase_dependent_rate = lambda: max(0,sample_backward_alpha() +
                                                   sample_backward_beta() * kinase)
k_t = 0 # for now
FAST = 10**6
