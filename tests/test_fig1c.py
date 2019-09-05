import sys
sys.path.insert(0, "..")

import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")

opt = network.get_mass_conservation_approach()

bounds, concentration_bounds = opt.get_optimization_bounds()

#bounds = [(1e-2, 1e2)]*12
#concentration_bounds = [(1e-2, 1e2)]*4

print("bounds")
print(bounds)
print("")
print("concentration_bounds")
print(concentration_bounds)
print("")

#params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=10,
#                                                                     concentration_bounds=concentration_bounds)

import numpy
#numpy.save('params.npy', params_for_global_min)
params_for_global_min = numpy.load('params.npy')

multistable_param_ind = opt.run_greedy_continuity_analysis(species="s15",
                                                           parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 'C3'})


#multistable_param_ind = opt.run_continuity_analysis(species='s15', parameters=params_for_global_min,
#                                                    auto_parameters={'PrincipalContinuationParameter': 'C3',
#                                                                     'RL0': 0.1, 'RL1': 30, 'DSMAX': 0.1})

opt.generate_report()

