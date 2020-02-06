import sys
sys.path.insert(0, "..")
import crnt4sbml
import numpy
import sympy


# 1.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/a_b.xml") # yes 10
# signal = "C1"
# response = "s5"
# iters = 100

# 2.
# network = crnt4sbml.CRNT("../sbml_files/hervagault_canu.xml") # yes with 100
# signal = "C1"
# response = "s1"
# iters = 10

# 3.
# network = crnt4sbml.CRNT("../sbml_files/p85-p110-PTEN.xml") # no with 10
# signal = "C4"
# response = "s37"
# iters = 100

# 4.
# network = crnt4sbml.CRNT("../sbml_files/Song.xml") # yes for 100
# signal = "C1"
# response = "s2"
# iters = 12

# #network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/open_fig5B_version2_clean.xml")  # graph is uniterminal
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/open_fig5B_version2.xml")
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/open_fig5B_version2_closed.xml")
# network = crnt4sbml.CRNT("../sbml_files/open_fig5B_no_sink.xml")

# 5.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/DoublePhos_DeadEnd.xml")
# signal = "C2"
# response = "s4"
# iters = 10

# 5 diff
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/DoublePhos_OGlcNAc.xml")
# signal = "C2"
# response = "s4"
# iters = 10

# 5. diff 2
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/DoublePhos_OGlcNAc_2.xml")
# signal = "C2"
# response = "s4"
# iters = 10


# 5. dou
network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/DoublePhos_DoublePath.xml")
signal = "C2"
response = "s4"
iters = 10

# 6.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/DoublePhos_DeadEnd_dbl_removed.xml")
# signal = "C2"
# response = "s5"
# iters = 1000

# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/irene_simple_bistable.xml")
# signal = "C1"
# response = "s1"
# iters = 10

network.basic_report()
network.print_c_graph()
#
# # network.get_network_graphml()
opt = network.get_mass_conservation_approach()
#
network.get_low_deficiency_approach().report_deficiency_one_theorem()
network.get_low_deficiency_approach().report_deficiency_zero_theorem()
sys.exit()
#

# print(network.get_c_graph().get_species())
# print(network.get_c_graph().get_reactions())

GA = network.get_general_approach(signal=signal, response=response, fix_reactions=True)

print(GA.get_conservation_laws())

sys.exit()

# 1.
# bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())

# 2.
# bnds = [(1e-4, 1e2)]*len(GA.get_input_vector())

# 3.
# bnds = GA.get_optimization_bounds()

# 4.
# bnds = [(1e-3, 6.0)]*len(network.get_c_graph().get_reactions()) + [(1e-3, 1000.0)]*len(network.get_c_graph().get_species())

# 5.
# bnds = GA.get_optimization_bounds()

# 6.
# bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())
bnds = GA.get_optimization_bounds()

# params, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=True,
#                                            dual_annealing_iters=5000, confidence_level_flag=True)

# params, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=False,
#                                            dual_annealing_iters=1000, confidence_level_flag=True)

# numpy.save('params.npy', params)

# numpy.save('params_DoublePhos_DeadEnd_removed.npy', params)
# params = numpy.load('params_DoublePhos_DeadEnd.npy')

# sys.exit()
params, obj_fun_vals, my_rank = GA.run_mpi_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=False,
                                                        dual_annealing_iters=1000, confidence_level_flag=True)

# params, obj_fun_vals, my_rank = GA.run_mpi_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=False,
#                                                         dual_annealing_iters=1000, confidence_level_flag=True)
GA.generate_report()
if my_rank == 0:
    #numpy.save('params_DoublePhos_DeadEnd_removed.npy', params)
    #print(params)
    print(obj_fun_vals)

# print(params)
# print(obj_fun_vals)

#sys.exit()

# multistable_param_ind, sample_portion, plot_specifications = GA.run_mpi_continuity_analysis(species=response, parameters=params, print_lbls_flag=False,
#                                                                                             auto_parameters={'PrincipalContinuationParameter': signal}, dir_path="./num_cont_graphs_parallel")

multistable_param_ind, sample_portion, plot_specifications = GA.run_mpi_greedy_continuity_analysis(species=response, parameters=params, print_lbls_flag=False,
                                                                                   auto_parameters={'PrincipalContinuationParameter': signal}, dir_path="./num_cont_graphs_parallel")

# multistable_param_ind, sample_portion, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params, print_lbls_flag=False,
#                                                                                    auto_parameters={'PrincipalContinuationParameter': signal})

# multistable_param_ind, plot_specifications = GA.run_continuity_analysis(species=response, parameters=params, print_lbls_flag=False,
#                                                                                         auto_parameters={'PrincipalContinuationParameter': signal})

# multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params, print_lbls_flag=True,
#                                                                                auto_parameters={'PrincipalContinuationParameter': signal}, dir_path="./num_cont_DoublePhos_DeadEnd_removed")#,
#                                                                                #plot_labels=['Rtot', 'S1*', None])


GA.generate_report()

