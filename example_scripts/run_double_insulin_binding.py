import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/double_insulin_binding.xml")

network.basic_report()

network.print_c_graph()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

approach = network.get_mass_conservation_approach()

print("Decision Vector:")
print(approach.get_decision_vector())
print("")

print("Species for concentration bounds:")
print(approach.get_concentration_bounds_species())

bounds, concentration_bounds = approach.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds,
                                                                          concentration_bounds=concentration_bounds,
                                                                          iterations=100)

multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s5", parameters=params_for_global_min,
                                                                                     auto_parameters={'PrincipalContinuationParameter': 'C2'})

approach.generate_report()
