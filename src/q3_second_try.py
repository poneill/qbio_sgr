from ssa import SSA
FAST = 10**2
kr = 1/100.0
gr = 1
gp = FAST
kp = 500 * gp

k3 = 10**10
gnxp_reactions = [((1,0,0),kr), #mRNA production
                  ((-1,0,0),gr), #mRNA degradation
                  ((-1,0,1),kp), #protein complex formation
                  ((1,1,-1),k3), #protein formation we do things this way
                  #to pack everything into mass action kinetics
                  ((0,-1,0),gp), #protein degradation
                  ]
init_state = (0,0,0)
gnxp_sim = SSA(gnxp_reactions,init_state)
gnxp_sim.verbose = True
gnxp_sim.run(100)

print "loaded"
