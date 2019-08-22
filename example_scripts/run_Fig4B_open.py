import crnt4sbml_test

network = crnt4sbml_test.CRNT("../sbml_files/Fig4B_open.xml")

network.basic_report()

network.print_c_graph()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

opt = network.get_semi_diffusive_approach()
print("")

opt.print_decision_vector()

print("Key species:")
print(opt.get_key_species())
print("")

print("Non key species:")
print(opt.get_non_key_species())

print("")
print("Boundary species:")
print(opt.get_boundary_species())

bounds = [(1e-2, 1e2)]*11

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=10000)

opt.generate_report()
