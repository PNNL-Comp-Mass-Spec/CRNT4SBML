import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Cii.xml")

opt = network.get_semi_diffusive_approach()

bounds = [(1e-3, 1e2)]*12

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds)

multistable_param_ind = opt.run_greedy_continuity_analysis(species="s7", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 're17'})

opt.generate_report()
