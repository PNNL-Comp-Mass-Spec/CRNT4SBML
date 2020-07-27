import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")

network.print_c_graph()

approach = network.get_mass_conservation_approach()

print("species for full system")
print(network.get_c_graph().get_species())
print("")
print("full system")
print(network.get_c_graph().get_ode_system())
print("")


bounds, concentration_bounds = approach.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds,
                                                                          concentration_bounds=concentration_bounds)

multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
                                                                                     auto_parameters={'PrincipalContinuationParameter': 'C3'})

print("independent species")
print(approach.get_independent_species())
print("")
print("independent system")
print(approach.get_independent_odes())

approach.generate_report()
