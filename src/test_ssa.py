from ssa import SSA

k1,k2,g1,g2,k3 = 10,10,1,1,1000
gnxp_reactions = [([1,0,0],k1), #mRNA production
                  ([-1,0,0],g1), #mRNA degradation
                  ([-1,0,1],k2), #protein complex formation
                  ([1,1,-1],k3), #protein formation we do things this way
                  #to pack everything into mass action kinetics
                  ([0,-1,0],g2), #protein degradation
              ]
init_state = [0,0,0]
gnxp_sim = SSA(gnxp_reactions,init_state)
gnxp_sim.run(1000)
