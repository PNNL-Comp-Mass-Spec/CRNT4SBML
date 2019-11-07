import crnt4sbml
import numpy

network = crnt4sbml.CRNT("../sbml_files/p85-p110-PTEN.xml")

network.basic_report()

network.print_c_graph()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

print("")

approach = network.get_mass_conservation_approach()

print("Decision Vector:")
print(approach.get_decision_vector())
print("")

print("Species for concentration bounds:")
print(approach.get_concentration_bounds_species())

bounds, concentration_bounds = approach.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds, iterations=5000,
                                                                          concentration_bounds=concentration_bounds)

approach.generate_report()

#numpy.save('def_params.npy', params_for_global_min)
#params_for_global_min = numpy.load('def_params.npy')

multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s37", parameters=params_for_global_min,
                                                                                     auto_parameters={'PrincipalContinuationParameter': "C4"})

# multistable_param_ind, plot_specifications = approach.run_continuity_analysis(species="s37", parameters=params_for_global_min,
#                                                                               auto_parameters={'DSMAX': 1e3,
#                                                                                                'PrincipalContinuationParameter': "C4",
#                                                                                                'RL0': 5e-2, 'RL1': 5e6, 'A0': 1e-6, 'A1': 1e6},
#                                                                               print_lbls_flag=True)

approach.generate_report()
