import sys
sys.path.insert(0, "..")
import crnt4sbml
import numpy
import sympy


# 1.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_1.xml")  # yes
# signal = "C1"
# response = "s11"
# iters = 50
# d_iters = 100
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species()))

# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_4c.xml")
# signal = "C1"
# response = "s11"
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species()))

# 2.
# network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml") # yes 500
# signal = "C3"
# response = "s15"
# iters = 500
# d_iters = 100
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species()))

# 3.
# network = crnt4sbml.CRNT("../sbml_files/closed_fig5A.xml") # yea
# signal = "C2"
# response = "s9"
# iters = 1000
# d_iters = 1000
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species()))

# 4.
network = crnt4sbml.CRNT("../sbml_files/irene2014.xml") # yes
signal = "C1"
response = "s1"
iters = 100
d_iters = 1000
bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species()))

# 5.
# network = crnt4sbml.CRNT("../sbml_files/irene2009.xml") # yes
# signal = "C1"
# response = "s3"
# iters = 100
# d_iters = 1000
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species()))

# 6.
# network = crnt4sbml.CRNT("../sbml_files/conradi2007.xml")
# signal = "C2"
# response = "s1"
# iters = 50
# d_iters = 1000
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species()))

# 7.
# network = crnt4sbml.CRNT("../sbml_files/Song.xml") # yes
# signal = "C1"
# response = "s2"
# iters = 20
# d_iters = 1000
# bnds = [(1e-3, 6.0)]*len(network.get_c_graph().get_reactions()) + [(1e-3, 1000.0)]*len(network.get_c_graph().get_species())

# 8.
# network = crnt4sbml.CRNT("../sbml_files/hervagault_canu.xml") # yes
# signal = "C1"
# response = "s1"
# iters = 10
# d_iters = 1000
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species()))

# 9.
# network = crnt4sbml.CRNT("../sbml_files/DoublePhos.xml") # yes
# signal = "C2"
# response = "s4"
# iters = 10
# d_iters = 1000
# bnds = GA.get_optimization_bounds()

# 10.
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts.xml")
# signal = "C2"
# response = "s11"
# iters = 75
# d_iters = 1000
# bnds = [(2.4, 2.42), (27.5, 28.1), (2.0, 2.15), (48.25, 48.4), (0.5, 1.1), (1.8, 2.1), (17.0, 17.5), (92.4, 92.6), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + \
#        [(0.0, 100.0), (18.0, 18.5), (0.0, 100.0), (0.0, 100.0), (27.0, 27.1), (8.2, 8.3), (90.0, 90.1), (97.5, 97.9), (30.0, 30.1)]

# network.basic_report()

GA = network.get_general_approach()

print(GA.get_conservation_laws())


# 9.
# bnds = GA.get_optimization_bounds()

GA.initialize_general_approach(signal=signal, response=response)

cons = [] #[{'type': 'ineq', 'fun': lambda x:  x[9] - 2.0*x[8]}, {'type': 'eq', 'fun': lambda x:  x[16]}]

# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=True,
#                                                           dual_annealing_iters=d_iters, confidence_level_flag=True,
#                                                           constraints=cons)
#
# GA.generate_report()

# numpy.save('./num_cont_direct_2/params.npy', params_for_global_min)

# numpy.save('./num_cont_direct_2/params_Fig1Ci.npy', params_for_global_min)

# numpy.save('./num_cont_direct_2/params_closed_fig5A.npy', params_for_global_min)

# numpy.save('./num_cont_direct_2/params_irene2014.npy', params_for_global_min)

# numpy.save('./num_cont_direct_2/params_irene2009.npy', params_for_global_min)

# numpy.save('./num_cont_direct_2/params_conradi2007.npy', params_for_global_min)

# numpy.save('./num_cont_direct_2/params_song.npy', params_for_global_min)

# numpy.save('./num_cont_direct_2/params_hervagault_canu.npy', params_for_global_min)

# numpy.save('./num_cont_direct_2/params_DoublePhos.npy', params_for_global_min)

# numpy.save('./num_cont_direct_2/params_Nuts.npy', params_for_global_min)

# numpy.save('./num_cont_lagrangian/params.npy', params_for_global_min)

# params_for_global_min = numpy.load('./num_cont_lagrangian/params.npy')


# sys.exit()

# print(network.get_c_graph().get_species())
# print(params_for_global_min)

# params_for_global_min = numpy.load('./num_cont_direct/params.npy')
# path = './num_cont_direct'

# params_for_global_min = numpy.load('./num_cont_direct_2/params_sub_mod_1.npy')

# params_for_global_min = numpy.load('./num_cont_direct_2/params_Fig1Ci.npy')

# params_for_global_min = numpy.load('./num_cont_direct_2/params_closed_fig5A.npy')

params_for_global_min = numpy.load('./num_cont_direct_2/params_irene2014.npy')

# params_for_global_min = numpy.load('./num_cont_direct_2/params_irene2009.npy')

# params_for_global_min = numpy.load('./num_cont_direct_2/params_conradi2007.npy')

# params_for_global_min = numpy.load('./num_cont_direct_2/params_song.npy')

# params_for_global_min = numpy.load('./num_cont_direct_2/params_hervagault_canu.npy')

# params_for_global_min = numpy.load('./num_cont_direct_2/params_DoublePhos.npy')

# params_for_global_min = numpy.load('./num_cont_direct_2/params_Nuts.npy')


path = './num_cont_direct_2'

# params_for_global_min = [params_for_global_min[0], params_for_global_min[1], params_for_global_min[2], params_for_global_min[3]]

# print(params_for_global_min)

# sys.exit()

GA.run_direct_simulation(params_for_global_min, path, change_in_relative_error=1e-4)


# multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
#                                                                                auto_parameters={'PrincipalContinuationParameter': signal},
#                                                                                dir_path=path)

# multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
#                                                                         auto_parameters={'PrincipalContinuationParameter': signal, 'ISW': -1, 'ISP': 0},
#                                                                         dir_path='./num_cont_lagrangian')