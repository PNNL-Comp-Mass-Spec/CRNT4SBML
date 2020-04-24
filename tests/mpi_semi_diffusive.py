import numpy
import pandas
import sympy
import time
import sys

sys.path.insert(0, "..")
import crnt4sbml

# 1.
#network = crnt4sbml.CRNT("../sbml_files/open_fig5B.xml")

# 2.
# network = crnt4sbml.CRNT("../sbml_files/open_fig5A.xml")

network = crnt4sbml.CRNT("../sbml_files/Fig1Cii.xml")

# optimization approach
opt = network.get_semi_diffusive_approach()

# the decision vector
#opt.get_decision_vector()

# overwriting specie concentration bounds for s4. Concentrations are in pM.
# opt.get_concentration_bounds_species()

# 1.
# bounds = [(1e-3, 1e2)]*17
# iters = 50
# response = "s5"
# signal = 're2'

# 2.
bounds = opt.get_optimization_bounds()
iters = 10 #500
# response = "s4"
# signal = 're3'

params_for_global_min, obj_fun_val_for_params, my_rank = opt.run_mpi_optimization(bounds=bounds, iterations=iters, confidence_level_flag=False)

if my_rank == 0:
    numpy.save('params.npy', params_for_global_min)

#opt.generate_report()

#sys.exit()

# The reponse-related specie should be picked based on CellDesigner IDs. In our case phoshorylated A is s2.
# How to pick continuation parameter? In our case it is the amount of A protein, thus the conservation law 3.
#print(opt.get_conservation_laws())

# multistable_param_ind, sample_points, plot_specifications = opt.run_mpi_continuity_analysis(species=response, parameters=params_for_global_min,
#                                                                          auto_parameters={'PrincipalContinuationParameter': signal,
#                                                                                           'RL0': 0.1, 'RL1': 30, 'A0': 0.0, 'A1': 10000},
#                                                                          dir_path="./num_cont_graphs_parallel")

# multistable_param_ind, sample_point, plot_specifications = opt.run_mpi_greedy_continuity_analysis(species=response, parameters=params_for_global_min,
#                                                            auto_parameters={'PrincipalContinuationParameter': signal},
#                                                            print_lbls_flag=False, dir_path="./num_cont_graphs_parallel")

opt.generate_report()
