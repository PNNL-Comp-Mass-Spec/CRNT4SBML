import crnt4sbml
import numpy

c = crnt4sbml.CRNT("../sbml_files/Fig1Cii.xml")

c.basic_report()

ldt = c.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

opt = c.get_semi_diffusive_approach()
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

bounds = [(1e-3, 1e2)]*12

#params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=100)

#numpy.save('inj_params.npy', params_for_global_min)
params_for_global_min = numpy.load('inj_params.npy')

#multistable_param_ind, plot_specifications = opt.run_continuity_analysis(species="s7", parameters=params_for_global_min,
#                                                                         auto_parameters={'PrincipalContinuationParameter': 're17',
#                                                                                          'RL0': 0.1, 'RL1': 100, 'A0': 0.0, 'A1': 10000})

multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s7", parameters=params_for_global_min,
                                                                                auto_parameters={'PrincipalContinuationParameter': 're17'})

opt.generate_report()
