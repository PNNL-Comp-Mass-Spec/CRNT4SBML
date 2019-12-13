import sys
sys.path.insert(0, "..")
import crnt4sbml
import numpy


# 1.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/a_b.xml") # yes 10
# signal = "C1"
# response = "s6"
# iters = 10

# 2.
network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml") # yes 10
signal = "C3"
response = "s15"
iters = 10

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
# iters = 50

# 7.
# network = crnt4sbml.CRNT("../sbml_files/conradi2007.xml") # yes with 20
# signal = "C2"
# response = "s1"
# iters = 20

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

network.basic_report()
network.print_c_graph()

GA = network.get_general_approach(signal=signal, response=response)

# 1.
#bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())

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


params, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=False, dual_annealing_iters=100)


#numpy.save('params_intuitive.npy', params)

#params = numpy.load('params_intuitive.npy')

# print("params")
# print(params)

multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal}, print_lbls_flag=False)


