import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/double_insulin_binding.xml")

network.basic_report()

network.print_c_graph()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

opt = network.get_mass_conservation_approach()

print("Decision Vector:")
print(opt.get_decision_vector())
print("")

print("Species for concentration bounds:")
print(opt.get_concentration_bounds_species())

# bounds = [(1e-8, 1e-4), (1e-5, 1e-3), (1e-8, 1e-4), (1e-5, 1e-3), (1e-8, 1e-4), (1e-5, 1e-3), (1e-3, 1.0), (1e-8, 1e-4),
#           (1e-5, 1e-3), (1e-3, 1.0), (1e-3, 1.0), (5e-1, 5e5), (5e-1, 5e5), (5e-1, 5e5)]
# concentration_bounds = [(5e-1, 5e5)]*5

bounds, concentration_bounds = opt.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds,
                                                                     iterations=100)

opt.generate_report()


multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s5", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 'C2'},
                                                           print_lbls_flag=True)

# multistable_param_ind, plot_specifications = opt.run_continuity_analysis(species="s5",
#                                                                          parameters=params_for_global_min,
#                                                                          auto_parameters={'DSMAX': 1e3,
#                                                                                           'PrincipalContinuationParameter': "C2",
#                                                                                           'RL0': 1e2, 'RL1': 1e6,
#                                                                                           'A0': 0, 'A1': 5e6})

opt.generate_report()
