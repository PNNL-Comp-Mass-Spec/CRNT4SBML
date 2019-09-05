import sys
sys.path.insert(0, "..")

import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Cii.xml")

opt = network.get_semi_diffusive_approach()

bounds = opt.get_optimization_bounds()
print("bounds")
print(bounds)

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=10, seed=0)

import numpy
numpy.save('params.npy',params_for_global_min)
#params_for_global_min = numpy.load('params.npy')

multistable_param_ind = opt.run_greedy_continuity_analysis(species="s7", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 're17'})

opt.generate_report()





