import sys
sys.path.insert(0, "..")
import crnt4sbml
import numpy
import sympy

network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_1.xml")

# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_4c.xml")

signal = "C1"
response = "s11"

network.basic_report()

GA = network.get_general_approach()  # (signal=signal, response=response, fix_reactions=False)

print(GA.get_conservation_laws())

GA.initialize_general_approach(signal=signal, response=response)

# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species()) - 1) + [(0.0, 1e-5)]

bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(0.0, 100.0)]*(len(network.get_c_graph().get_species()))

# sympy.pprint(network.get_c_graph().get_ode_system())

# basinhopping
# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=10, seed=0, print_flag=True,
#                                                           feasible_iters=2, main_iters=10, confidence_level_flag=True,
#                                                           method='basin')

# dual anealing
params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=10, seed=0, print_flag=True,
                                                          feasible_iters=100, main_iters=100, confidence_level_flag=True,
                                                          method='dual')



GA.generate_report()

numpy.save('./num_cont_direct/params.npy', params_for_global_min)

# params_for_global_min = numpy.load('./num_cont_lagrangian/params.npy')


# sys.exit()

print(network.get_c_graph().get_species())
print(params_for_global_min)

multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal},
                                                                               dir_path='./num_cont_lagrangian')

# multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
#                                                                         auto_parameters={'PrincipalContinuationParameter': signal, 'ISW': -1, 'ISP': 0},
#                                                                         dir_path='./num_cont_lagrangian')
