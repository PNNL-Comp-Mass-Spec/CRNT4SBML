import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Cii.xml")

approach = network.get_semi_diffusive_approach()

bounds = approach.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds)

multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s7", parameters=params_for_global_min,
                                                                                     auto_parameters={'PrincipalContinuationParameter': 're17'})

approach.generate_report()
