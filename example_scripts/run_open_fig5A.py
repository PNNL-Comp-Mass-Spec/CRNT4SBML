import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/open_fig5A.xml")

network.basic_report()

network.print_c_graph()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

approach = network.get_semi_diffusive_approach()
print("")

approach.print_decision_vector()

print("Key species:")
print(approach.get_key_species())
print("")

print("Non key species:")
print(approach.get_non_key_species())

print("")
print("Boundary species:")
print(approach.get_boundary_species())

bounds = approach.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds, iterations=500)

multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s4", parameters=params_for_global_min,
                                                                                     auto_parameters={'PrincipalContinuationParameter': 're3'})

approach.generate_report()
