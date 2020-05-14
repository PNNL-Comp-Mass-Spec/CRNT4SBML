import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")

approach = network.get_mass_conservation_approach()

bounds, concentration_bounds = approach.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds,
                                                                          concentration_bounds=concentration_bounds)

multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
                                                                                     auto_parameters={'PrincipalContinuationParameter': 'C3'})

approach.generate_report()
