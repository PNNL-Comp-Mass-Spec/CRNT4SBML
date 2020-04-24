import numpy
#import pandas
import sympy
import time
import sys

sys.path.insert(0, "..")
import crnt4sbml


network = crnt4sbml.CRNT("../sbml_files/DoublePhos.xml")

# network.basic_report()
#
# network.print_c_graph()
#
# ldt = network.get_low_deficiency_approach()
# ldt.report_deficiency_zero_theorem()
# ldt.report_deficiency_one_theorem()
#
opt = network.get_mass_conservation_approach()
#
# print("Decision Vector:")
# print(opt.get_decision_vector())
# print("")
#
# print("Species for concentration bounds:")
# print(opt.get_concentration_bounds_species())

print(opt.get_conservation_laws())

sys.exit()

bounds, concentration_bounds = opt.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = opt.run_mpi_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds,
                                                                     iterations=10)

opt.generate_report()


multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s4", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 'C2'})

opt.generate_report()










sys.exit()

# network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")
#
# opt = network.get_mass_conservation_approach()

# bounds, concentration_bounds = opt.get_optimization_bounds()
#
# params_for_global_min, obj_fun_val_for_params, my_rank = opt.run_mpi_optimization(bounds=bounds,
#                                                                                   concentration_bounds=concentration_bounds)
#
# if my_rank == 0:
#     print(params_for_global_min)
#
# multistable_param_ind, sample_points, plot_specifications = opt.run_mpi_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
#                                                                                                    auto_parameters={'PrincipalContinuationParameter': 'C3'})
#
# opt.generate_report()
#
#
# sys.exit()

# 1.
# network = crnt4sbml.CRNT("../sbml_files/DoublePhos.xml")
# signal = "C2"
# response = "s4"
# iters = 100

# 2.
# network = crnt4sbml.CRNT("../sbml_files/Fig4C_closed.xml")

# 3.
# network = crnt4sbml.CRNT("../sbml_files/closed_fig5A.xml")
# signal = "C2"
# response = "s9"

# 4.
# network = crnt4sbml.CRNT("../sbml_files/irene2014.xml")
# signal = "C1"
# response = "s1"

# 5.
# network = crnt4sbml.CRNT("../sbml_files/irene2009.xml")

# 6.
# network = crnt4sbml.CRNT("../sbml_files/hervagault_canu.xml")

# 7.
# network = crnt4sbml.CRNT("../sbml_files/conradi2007.xml")

# 8.
# network = crnt4sbml.CRNT("../sbml_files/double_insulin_binding.xml")
# signal = "C2"
# response = "s5"

# 10.
# network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")
# signal = 'C3'
# response = 's15'

# optimization approach
opt = network.get_mass_conservation_approach()

# the decision vector
#opt.get_decision_vector()

# overwriting specie concentration bounds for s4. Concentrations are in pM.
# opt.get_concentration_bounds_species()

# 1.
# bounds, concentration_bounds = opt.get_optimization_bounds()

# 2.
# bounds, concentration_bounds = opt.get_optimization_bounds()
# iters = 10

# 3.
# bounds = [(1e-2, 1e2)]*12
# concentration_bounds = [(1e-2, 1e2)]*6
# iters = 100

# 4.
bounds = [(1e-2, 1e2)]*12
concentration_bounds = [(1e-2, 1e2)]*5
iters = 100

# 5.
# bounds = [(1e-2, 1e2)]*10
# concentration_bounds = [(1e-2, 1e2)]*4
# iters = 100

# 6.
# bounds = [(1e-2, 1e2)]*13
# concentration_bounds = [(1e-2, 1e2)]*4
# iters = 100

# 7.
# bounds = [(1e-2, 1e2)]*20
# concentration_bounds = [(1e-2, 1e2)]*7
# iters = 100

# 8.
# bounds, concentration_bounds = opt.get_optimization_bounds()
# iters = 100

# params_for_global_min, obj_fun_val_for_params, my_rank = opt.run_mpi_optimization(bounds=bounds, iterations=iters, confidence_level_flag=False,
#                                                                                   concentration_bounds=concentration_bounds)

# params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=iters, confidence_level_flag=False,
#                                                                      concentration_bounds=concentration_bounds)

if my_rank == 0:
    numpy.save('params.npy', params_for_global_min)

params_for_global_min = numpy.load('params.npy')

#opt.generate_report()

#sys.exit()

# The reponse-related specie should be picked based on CellDesigner IDs. In our case phoshorylated A is s2.
# How to pick continuation parameter? In our case it is the amount of A protein, thus the conservation law 3.
#print(opt.get_conservation_laws())

# multistable_param_ind, sample_poriton, plot_specifications = opt.run_mpi_greedy_continuity_analysis(species=response, parameters=params_for_global_min,
#                                                            auto_parameters={'PrincipalContinuationParameter': signal}, dir_path="./num_cont_graphs_parallel",
#                                                            print_lbls_flag=False)

multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min[[58]],
                                                           auto_parameters={'PrincipalContinuationParameter': signal}, dir_path="./num_cont_graphs_parallel",
                                                           print_lbls_flag=True)

# multistable_param_ind, sample_poriton, plot_specifications = opt.run_mpi_continuity_analysis(species=response, parameters=params_for_global_min,
#                                                            auto_parameters={'PrincipalContinuationParameter': signal,
#                                                                             'RL0': 1e2, 'RL1': 1e6, 'A0': 0.0, 'A1': 5e6,
#                                                                             'DSMAX': 1e3}, dir_path="./num_cont_graphs_parallel",
#                                                            print_lbls_flag=False)

# multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s4", parameters=params_for_global_min,
#                                                            auto_parameters={'PrincipalContinuationParameter': 'C2'},
#                                                            print_lbls_flag=False)

opt.generate_report()
