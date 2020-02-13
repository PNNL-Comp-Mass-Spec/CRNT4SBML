import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/DoublePhos.xml")

network.basic_report()

network.print_c_graph()

network.get_network_graphml()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

opt = network.get_mass_conservation_approach()

print("Decision Vector:")
print(opt.get_decision_vector())
print("")

print("Species for concentration bounds:")
print(opt.get_concentration_bounds_species())

bounds, concentration_bounds = opt.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds,
                                                                     iterations=100)

opt.generate_report()


multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s4", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 'C2'})

opt.generate_report()
