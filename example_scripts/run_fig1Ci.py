import crnt4sbml

c = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")

opt = c.get_mass_conservation_approach()

print(opt.get_decision_vector())

bnds = [(1e-2, 1e2)]*12

conc_bnds = [(1e-2, 1e2)]*7

num_itr = 100

import numpy

sys_min = numpy.finfo(float).eps
sd = 0
prnt_flg = False
num_dtype = numpy.float64

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bnds, concentration_bounds=conc_bnds,
                                                                     iterations=num_itr, seed=sd, print_flag=prnt_flg,
                                                                     numpy_dtype=num_dtype, sys_min_val=sys_min)

numpy.save('params.npy',params_for_global_min)
#params_for_global_min = numpy.load('params.npy')

print(opt.get_conservation_laws())

spcs = "s15"
PCP_x = "C3"

#multistable_param_ind = opt.run_continuity_analysis(species=spcs, parameters=params_for_global_min,
#                                                    auto_parameters={'PrincipalContinuationParameter': PCP_x,
#                                                                     'RL0': 0.1, 'RL1': 30, 'DSMAX': 0.1},
#                                                    print_lbls_flag=False, dir_path="./stability_graphs")

multistable_param_ind = opt.run_greedy_continuity_analysis(species=spcs, parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': PCP_x})

opt.generate_report()


