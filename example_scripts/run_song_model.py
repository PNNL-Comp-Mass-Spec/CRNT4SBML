import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Song.xml")

network.basic_report()

network.print_c_graph()

network.get_network_graphml()

approach = network.get_general_approach()
approach.initialize_general_approach(signal="C1", response="s2", fix_reactions=True)

print(approach.get_input_vector())

bnds = [(1e-3, 6.0)]*len(network.get_c_graph().get_reactions()) + [(1e-3, 1000.0)]*len(network.get_c_graph().get_species())

params, obj_fun_vals = approach.run_optimization(bounds=bnds, iterations=100)

multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s2", parameters=params,
                                                                                     auto_parameters={'PrincipalContinuationParameter': "C1"})

approach.generate_report()

