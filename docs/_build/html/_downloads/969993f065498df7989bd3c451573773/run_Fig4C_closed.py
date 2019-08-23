import crnt4sbml_test

network = crnt4sbml_test.CRNT("../sbml_files/Fig4C_closed.xml")

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

bounds = [(1e-2, 1e2)]*10
concentration_bounds = [(1e-2, 1e2)]*3

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds,
                                                                     iterations=10000)

opt.generate_report()
