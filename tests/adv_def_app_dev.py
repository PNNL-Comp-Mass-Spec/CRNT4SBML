import sys
sys.path.insert(0, "..")

import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")

opt = network.get_mass_conservation_approach()

bounds, concentration_bounds = opt.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds)

multistable_param_ind = opt.run_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 'C3'})

opt.generate_report()


# network = crnt4sbml.CRNT("../sbml_files/HDA_3_2_1.xml")
#
# approach = network.get_advanced_deficiency_approach()
#
# approach.run_higher_deficiency_algorithm()

