import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/open_fig5A.xml")

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

bounds = opt.get_optimization_bounds() #[(1e-4, 1e2)]*12

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=500)

import numpy
numpy.save('params.npy', params_for_global_min)


multistable_param_ind = opt.run_greedy_continuity_analysis(species="s4", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 're3'})

opt.generate_report()
