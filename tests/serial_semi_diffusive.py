import numpy
import pandas
import sympy
import time
import sys

sys.path.insert(0, "..")
import crnt4sbml

# 1.
network = crnt4sbml.CRNT("../sbml_files/open_fig5B.xml")

# optimization approach
opt = network.get_semi_diffusive_approach()

# the decision vector
#opt.get_decision_vector()

# overwriting specie concentration bounds for s4. Concentrations are in pM.
# opt.get_concentration_bounds_species()

# 1.
bounds = [(1e-3, 1e2)]*17
iters = 50
response = "s5"
signal = 're2'


params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=iters, confidence_level_flag=False)


#print(params_for_global_min)
print(obj_fun_val_for_params)
print(len(obj_fun_val_for_params))

#opt.generate_report()

# multistable_param_ind, plot_specifications = opt.run_continuity_analysis(species=response, parameters=params_for_global_min,
#                                                                          auto_parameters={'PrincipalContinuationParameter': signal,
#                                                                                           'RL0': 0.1, 'RL1': 30, 'A0': 0.0, 'A1': 10000},
#                                                                          dir_path="./num_cont_graphs_serial")

multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': signal},
                                                           print_lbls_flag=False, dir_path="./num_cont_graphs_serial")


opt.generate_report()
