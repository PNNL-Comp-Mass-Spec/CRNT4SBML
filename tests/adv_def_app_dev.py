import sys
sys.path.insert(0, "..")

import crnt4sbml


network = crnt4sbml.CRNT("../sbml_files/conradi_paper.xml")
#network = crnt4sbml.CRNT("../sbml_files/Song2.xml")
#network = crnt4sbml.CRNT("../sbml_files/GuanyuWang.xml")

network.basic_report()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

opt = network.get_mass_conservation_approach()

#sys.exit()

bounds, concentration_bounds = opt.get_optimization_bounds()

print(bounds)
print(concentration_bounds)

print(opt.get_conservation_laws())
params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=1000,
                                                                     concentration_bounds=concentration_bounds)

# multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s1", parameters=params_for_global_min,
#                                                                                 auto_parameters={'PrincipalContinuationParameter': 'C2'})

# multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s3", parameters=params_for_global_min,
#                                                                                 auto_parameters={'PrincipalContinuationParameter': 'C2'})

# multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s1", parameters=params_for_global_min,
#                                                                                 auto_parameters={'PrincipalContinuationParameter': 'C1'})

# import sys
# sys.exit()

multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s2", parameters=params_for_global_min,
                                                                                auto_parameters={'PrincipalContinuationParameter': 'C1'})

# multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s7", parameters=params_for_global_min,
#                                                                                 auto_parameters={'PrincipalContinuationParameter': 'C4'})

opt.generate_report()

# network = crnt4sbml.CRNT("../sbml_files/HDA_3_2_1.xml")
#
# approach = network.get_advanced_deficiency_approach()
#
# approach.run_higher_deficiency_algorithm()

