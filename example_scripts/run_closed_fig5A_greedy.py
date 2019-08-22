import crnt4sbml_test

network = crnt4sbml_test.CRNT("../sbml_files/closed_fig5A.xml")

network.basic_report()

network.print_c_graph()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

print("")

opt = network.get_mass_conservation_approach()

print("Decision Vector:")
print(opt.get_decision_vector())

bounds = [(1e-2, 1e2)]*12
concentration_bounds = [(1e-2, 1e2)]*9

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds,
                                                                     iterations=100)

multistable_param_ind = opt.run_greedy_continuity_analysis(species="s9", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 'C2'})

opt.generate_report()
