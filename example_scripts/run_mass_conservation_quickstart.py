import crnt4sbml_test

network = crnt4sbml_test.CRNT("../sbml_files/Fig1Ci.xml")

opt = network.get_mass_conservation_approach()

bounds = [(1e-2, 1e2)]*12
concentration_bounds = [(1e-2, 1e2)]*4

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds,
                                                                     iterations=15)

multistable_param_ind = opt.run_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 'C3'})

opt.generate_report()
