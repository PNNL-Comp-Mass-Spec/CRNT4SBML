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
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/DoublePhos_DoublePath.xml")
# signal = "C2"
# response = "s4"
# iters = 10

# 6.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/DoublePhos_DeadEnd_dbl_removed.xml")
# signal = "C2"
# response = "s5"
# iters = 1000

# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/irene_simple_bistable.xml")
# signal = "C1"
# response = "s1"
# iters = 10

# network.basic_report()
# network.print_c_graph()
#
# # network.get_network_graphml()
# opt = network.get_mass_conservation_approach()
# #
# network.get_low_deficiency_approach().report_deficiency_one_theorem()
# network.get_low_deficiency_approach().report_deficiency_zero_theorem()
# sys.exit()
#

# 7.
network = crnt4sbml.CRNT("../sbml_files/DoublePhos.xml")
signal = 'C2'
response = 's4'

# 8.
# network = crnt4sbml.CRNT("../sbml_files/Fig4C_closed.xml")
# signal = 'C1'
# response = 's1'

# 9.
# network = crnt4sbml.CRNT("../sbml_files/closed_fig5A.xml")
# signal = 'C2'
# response = 's9'

# 10.
# network = crnt4sbml.CRNT("../sbml_files/irene2014.xml")
# signal = 'C1'
# response = 's1'

# 11.
# network = crnt4sbml.CRNT("../sbml_files/irene2009.xml")
# signal = 'C1'
# response = 's3'

# 12.
# network = crnt4sbml.CRNT("../sbml_files/hervagault_canu.xml")
# signal = 'C1'
# response = 's1'

# 13.
# network = crnt4sbml.CRNT("../sbml_files/conradi2007.xml")
# signal = 'C2'
# response = 's1'

# 14.
# network = crnt4sbml.CRNT("../sbml_files/double_insulin_binding.xml")
# signal = 'C2'
# response = 's5'

# print(network.get_c_graph().get_species())
# print(network.get_c_graph().get_reactions())

GA = network.get_general_approach(signal=signal, response=response, fix_reactions=True)

#print(GA.get_conservation_laws())

# sys.exit()

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
# bnds = GA.get_optimization_bounds()

# 7.
bnds = GA.get_optimization_bounds()
iters = 10

# 8.
# bnds = GA.get_optimization_bounds()
# iters = 10

# 9.
# bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())
# iters = 100

# 10.
# bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())
# iters = 100

# 11.
# bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())
# iters = 100

# 12.
# bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())
# iters = 100

# 13.
# bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())
# iters = 100

# 14.
# bnds = GA.get_optimization_bounds()
# iters = 100

params, obj_fun_vals, my_rank = GA.run_mpi_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=True,
                                                        dual_annealing_iters=100, confidence_level_flag=True)

if my_rank == 0:
    # numpy.save('params.npy', params)
    #print(params)
    print(obj_fun_vals)
    print(len(obj_fun_vals))

# params = numpy.load('params.npy')

# sys.exit()

multistable_param_ind, sample_portion, plot_specifications = GA.run_mpi_greedy_continuity_analysis(species=response, parameters=params, print_lbls_flag=False,
                                                                                   auto_parameters={'PrincipalContinuationParameter': signal}, dir_path="./num_cont_graphs_parallel")

GA.generate_report()

