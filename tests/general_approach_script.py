import sys
sys.path.insert(0, "..")
import crnt4sbml
import numpy
import sympy


# 1.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/a_b.xml") # yes 10
# signal = "C1"
# #response = "s6"
# response = "s5"
# iters = 10

# 2.
# network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml") # yes 10
# signal = "C3"
# response = "s15"
# iters = 10

# 3.
# network = crnt4sbml.CRNT("../sbml_files/closed_fig5A.xml") # yes 10
# signal = "C2"
# response = "s9"
# iters = 10

# 4.
# network = crnt4sbml.CRNT("../sbml_files/irene2014.xml") # yes 10
# signal = "C1"
# response = "s1"
# iters = 10

# 5.
# network = crnt4sbml.CRNT("../sbml_files/irene2009.xml") # yes 10
# signal = "C1"
# response = "s3"
# iters = 10

# 6.
# network = crnt4sbml.CRNT("../sbml_files/hervagault_canu.xml") # no with 50
# signal = "C1"
# response = "s1"
# iters = 10

# 7.
# network = crnt4sbml.CRNT("../sbml_files/conradi2007.xml") # yes with 20
# signal = "C2"
# response = "s1"
# iters = 10

# 8.
# network = crnt4sbml.CRNT("../sbml_files/double_insulin_binding.xml") # yes with 10 but not great
# signal = "C2"
# response = "s5"
# iters = 10

# 9.
# network = crnt4sbml.CRNT("../sbml_files/p85-p110-PTEN.xml") # no with 10
# signal = "C4"
# response = "s37"
# iters = 10

# 10.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/positive_mTORC2_v2.xml") # no with 10
# signal = "C1"
# response = "s6"
# iters = 10

# 11.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/unfolded_classic_positive_negative_3.xml") # no with 10
# # yes with [1e-2,1e2] with 10
# signal = "C2"
# response = "s6"
# iters = 10

# 12.
# network = crnt4sbml.CRNT("../sbml_files/Song.xml") # yes for 100
# signal = "C1"
# response = "s2"
# iters = 100

# 13.
# network = crnt4sbml.CRNT("../sbml_files/DoublePhos.xml") # yes 10
# signal = "C2"
# response = "s4"
# iters = 10

# import math
# a = 1
# b = 5
# r = 1
# n_til = 10
#
# a_bar = a + b - 1
# b_bar = b - r - 1
#
# prob = 1 - (math.factorial(n_til + a_bar)*math.factorial(2*n_til + b_bar))/(math.factorial(2*n_til + a_bar)*math.factorial(n_til + b_bar))
#
# print(prob)
#
# sys.exit()


# 14.
# network = crnt4sbml.CRNT("../sbml_files/attempt_at_irreversible.xml")              ###################################################
# signal = "C1"
# response = "s26"
# iters = 10

# network = crnt4sbml.CRNT("../sbml_files/model_e_v1.xml") #
# signal = "C1"
# response = "s1"
# iters = 10

# 16.
# network = crnt4sbml.CRNT("../sbml_files/Fig4B_closed.xml")
# network = crnt4sbml.CRNT("../sbml_files/Fig4B_open.xml")
# network = crnt4sbml.CRNT("../sbml_files/small_non_bistable.xml")
# signal = "C1"
# response = "s1"
# iters = 10

# 17.
network = crnt4sbml.CRNT("../sbml_files/zeroed_reactions.xml")
signal = "C1"
response = "s2"

# print(len(network.get_c_graph().get_species()))
# print(len(network.get_c_graph().get_complexes()))
# print(len(network.get_c_graph().get_linkage_classes()))
#
# sympy.pprint(network.get_c_graph().get_s().rref())
#
# print(network.get_c_graph().get_s().rank())
#
# sympy.pprint(network.get_c_graph().get_ode_system())
#
# #sys.exit()
# print(network.get_c_graph().get_deficiency())

GA = network.get_general_approach(signal=signal, response=response)

print(GA.get_conservation_laws())
print(GA.get_fixed_reactions())
print(GA.get_solutions_to_fixed_reactions())




sys.exit()

# opt = network.get_mass_conservation_approach()
#
# bounds, concentration_bounds = opt.get_optimization_bounds()
#
# bounds = [(1e-2, 100.0)]*len(bounds)
#
# concentration_bounds = [(1e-2, 100.0)]*len(concentration_bounds)
#
# params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=100,
#                                                                      concentration_bounds=concentration_bounds,
#                                                                      print_flag=True, confidence_level_flag=True,
#                                                                      change_in_rel_error=1e-2)

# opt = network.get_semi_diffusive_approach()
#
# bounds = opt.get_optimization_bounds()
#
# bounds = [(1e-2, 100.0)]*len(bounds)
#
# params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=100,
#                                                                      print_flag=True, confidence_level_flag=True,
#                                                                      change_in_rel_error=1e-2)


# GA = network.get_general_approach(signal=signal, response=response)
#
# print(GA.get_input_vector())
# # print(GA.get_decision_vector())
# # print(GA.get_fixed_reactions())
# # print(GA.get_solutions_to_fixed_reactions())
#
# print(GA.get_determinant_of_jacobian())
#
# bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())
#
#
# #sys.exit()
#
# params, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, dual_annealing_iters=1000,
#                                            confidence_level_flag=True)
# GA.generate_report()
#
# sys.exit()

network.basic_report()
network.print_c_graph()

# print(network.get_c_graph().get_species())
# sympy.pprint(network.get_c_graph().get_b())

# var1 = numpy.float64(1.120003)
# var2 = numpy.float64(1.12)
#
# var1 = numpy.float64(1.120003e-30)
# var2 = numpy.float64(1.12e-30)
#
#
# # print(var1)
# # print(numpy.log10(var1))
# # print(var2)
# # print(numpy.log10(var2))
#
# abs_diff = abs(var1 - var2)
# relative_abs_diff = abs(var1 - var2)/var1
# print(f"Absolute difference: {abs_diff}")
# print(f"Relative absolute difference: {relative_abs_diff}")
#
# relative_abs_diff = abs(var2 - var1)/var2
# print(f"Relative absolute difference: {relative_abs_diff}")
#
# sys.exit()

# opt = network.get_mass_conservation_approach()
#
# bounds, concentration_bounds = opt.get_optimization_bounds()
#
# bounds = [(1e-3, 100.0)]*len(bounds)
#
# concentration_bounds = [(1e-3, 100.0)]*len(concentration_bounds)
#
# params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=100,
#                                                                      concentration_bounds=concentration_bounds, print_flag=True)
#
# multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min,
#                                                                                 auto_parameters={'PrincipalContinuationParameter': signal})

#sys.exit()

GA = network.get_general_approach(signal=signal, response=response)


# 1.
# bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())
# print(GA.get_input_vector())
# print(GA.get_decision_vector())
# print(GA.get_fixed_reactions())
# print(GA.get_solutions_to_fixed_reactions())

#sys.exit()

# 1.
# bnds = [(1e-3, 100.0)]*len(network.get_c_graph().get_reactions()) + [(10.0, 100.0)] + [(-100.0, -1.0)] + [(-100.0, -1.0)]

# 2.
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(-10.0, 10.0)]*6 + [(1e-2, 100.0)]

# 4.
#bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 1000.0)] + [(-100.0, 10.0)]*6

# 5.
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(-100.0, 100.0)]*3 + [(1e-2, 1000.0)] + [(-100.0, 100.0)]

# 6.
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)] + [(-100.0, 100.0)]*6

# 7.
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)] + [(-100.0, 100.0)]*8

# 13.
#bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(-100.0, 100.0)]*2 + [(-100.0, 100.0)] + [(-10.0, 100.0)]*5

# 14.
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)] + [(-10.0, 10.0)] + [(-10.0, 10.0)]

# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(-10.0, 10.0)]*4 + [(1e-2, 10.0)] + [(-100.0, 10.0)]*4

#bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*1 + [(1e-2, 100.0)] + [(1e-2, 100.0)]*4

# 2.
bnds = GA.get_optimization_bounds()

# 3.
#bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())

# 4.
#bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())

# 5.
#bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())

# 6.
#bnds = [(1e-4, 1e2)]*len(GA.get_input_vector())

# 7.
#bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())

# 8.
#bnds = GA.get_optimization_bounds()

# 9.
#bnds = GA.get_optimization_bounds()

# 10.
#bnds = GA.get_optimization_bounds()

# 11.
#bnds = GA.get_optimization_bounds()

# 12.
#bnds = GA.get_optimization_bounds()
#bnds = [(1e-3, 6.0)]*len(network.get_c_graph().get_reactions()) + [(1e-3, 1000.0)]*len(network.get_c_graph().get_species())
#bnds = [(1e-3, 10.0)]*len(network.get_c_graph().get_reactions()) + [(1e-3, 1000.0)]*len(network.get_c_graph().get_species())

# # 13.
# bnds = GA.get_optimization_bounds()

# 14.
#bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(-100.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)]
#bnds = [(1e-2, 100.0)]*len(GA.get_input_vector())
# bnds = GA.get_optimization_bounds()

# 16.
#bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())

# print(GA.get_input_vector())
#
# sympy.pprint(GA.get_indpendent_odes_subs())
#
# sympy.pprint(sympy.expand(GA.get_indpendent_odes_subs()))
#
# #sys.exit()
#
# print(bnds)
# print(GA.get_decision_vector())


params, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=False,
                                           dual_annealing_iters=1000, confidence_level_flag=True)


#numpy.save('params.npy', params)
#params = numpy.load('params.npy')

# print(det_point_sets)
# numpy.savetxt("./data_for_PCA/end_dec_vec_1000_samples.txt", det_point_sets)
# numpy.savetxt("./data_for_PCA/beginning_dec_vec_1000_samples.txt", feasible_point_sets)
# numpy.savetxt("./data_for_PCA/obj_function_1000_samples.txt", obj_fun_vals)


#numpy.save('params_intuitive.npy', params)

#params = numpy.load('params_intuitive.npy')

# print("params")
# print(params)
multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal})

GA.generate_report()

