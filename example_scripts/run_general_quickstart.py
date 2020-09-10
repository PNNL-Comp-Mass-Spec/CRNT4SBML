import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")

approach = network.get_general_approach()

bnds = approach.get_optimization_bounds()

approach.initialize_general_approach(signal="C3", response="s15", fix_reactions=True)

params_for_global_min, obj_fun_vals = approach.run_optimization(bounds=bnds, dual_annealing_iters=100)

multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
                                                                                     auto_parameters={'PrincipalContinuationParameter': "C3"})

approach.generate_report()
