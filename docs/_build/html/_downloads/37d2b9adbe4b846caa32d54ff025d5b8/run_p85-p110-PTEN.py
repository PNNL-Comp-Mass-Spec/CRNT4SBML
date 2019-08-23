import crnt4sbml
import numpy

network = crnt4sbml.CRNT("../sbml_files/p85-p110-PTEN.xml")

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

bounds = [(1e-8, 1e-4),  # re1
          (1e-5, 1e-2),  # re1r
          (1e-8, 1e-4),  # re2
          (1e-5, 1e-2),  # re2r
          (1e-8, 1e-4),  # re3
          (1e-5, 1e-2),  # re3r
          (1e-8, 1e-4),  # re9
          (1e-5, 1e-2),  # re9r
          (1e-8, 1e-4),  # re10
          (1e-5, 1e-2),  # re10r
          (0.001, 1.0),  # re11
          (1e-8, 1e-4),  # re12
          (1e-5, 1e-2),  # re12r
          (0.001, 1.0),  # re13
          (1e-8, 1e-4),  # re14
          (1e-5, 1e-2),  # re14r
          (0.001, 1.0),  # re15
          (5e-1, 5e5),  # s8
          (5e-1, 5e5),  # s9
          (5e-1, 5e5),  # s23
          (5e-1, 5e5),  # s24
          (5e-1, 5e5)]  # s37
concentration_bounds = [(5e-1, 5e5)]*8

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds,
                                                                     iterations=100)

numpy.save('def_params.npy', params_for_global_min)
#params_for_global_min = numpy.load('def_params.npy')

#multistable_param_ind = opt.run_continuity_analysis(species="s37", parameters=params_for_global_min,
#                                                    auto_parameters={'DSMAX': 1e3,
#                                                                     'PrincipalContinuationParameter': "C4",
#                                                                     'RL0': 5e-2, 'RL1': 5e6, 'A0': 1e-6, 'A1': 1e6},
#                                                    print_lbls_flag=True)

multistable_param_ind = opt.run_greedy_continuity_analysis(species="s37", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': "C4"})

opt.generate_report()
