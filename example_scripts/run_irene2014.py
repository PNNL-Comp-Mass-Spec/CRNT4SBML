import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/irene2014.xml")

network.basic_report()

network.print_c_graph()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

print("")

opt = network.get_mass_conservation_approach()

print("Decision Vector:")
print(opt.get_decision_vector())
print("")

print("Species for concentration bounds:")
print(opt.get_concentration_bounds_species())

bounds = [(1e-2, 1e2)]*12
concentration_bounds = [(1e-2, 1e2)]*5

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds,
                                                                     iterations=100)

multistable_param_ind = opt.run_continuity_analysis(species="s1", parameters=params_for_global_min,
                                                    auto_parameters={'PrincipalContinuationParameter': 'C1'})

opt.generate_report()
