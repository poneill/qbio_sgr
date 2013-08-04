from ssa import SSA

kr = 1/100.0
gr = 10**4
kp = 500 * gr
gp = 6.79
k3 = 10**6
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
gnxp_sim.run(1000)

print "loaded"
