import sys
sys.path.insert(0, "..")

import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Cii.xml")

opt = network.get_semi_diffusive_approach()

opt.print_decision_vector()

print(opt.get_decision_vector())

bounds = opt.get_optimization_bounds()

print(opt.get_key_species())
print(opt.get_non_key_species())
print(opt.get_boundary_species())

import numpy
num_itr = 100
sys_min = numpy.finfo(float).eps
sd = 0
prnt_flg = False
num_dtype = numpy.float64

#params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=num_itr, seed=sd,
#                                                                     print_flag=prnt_flg, numpy_dtype=num_dtype,
#                                                                     sys_min_val=sys_min)

#numpy.save('params.npy', params_for_global_min)

params_for_global_min = numpy.load('params.npy')

opt.generate_report()

multistable_param_ind = opt.run_continuity_analysis(species='s7', parameters=params_for_global_min,
                                                    auto_parameters={'PrincipalContinuationParameter': 're17',
                                                                     'RL0': 0.1, 'RL1': 100, 'A0': 0.0,
                                                                     'A1': 10000})

# multistable_param_ind = opt.run_greedy_continuity_analysis(species='s7', parameters=params_for_global_min,
#                                                            auto_parameters={'PrincipalContinuationParameter': 're17'})

opt.generate_report()





