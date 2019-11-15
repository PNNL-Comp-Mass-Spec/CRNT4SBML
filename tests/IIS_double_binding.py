import sys
sys.path.insert(0, "..")
import crnt4sbml


network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/IIS_double_binding.xml")

#network.get_network_graphml()


# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/a_b.xml")

#network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/p85-p110-PTEN_v3.xml")

#network = crnt4sbml.CRNT("../sbml_files/closed_fig5A.xml")

print(network.get_c_graph().get_network_dimensionality_classification())

#sys.exit()

network.basic_report()

network.print_biological_reaction_types()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

# optimization approach
opt = network.get_mass_conservation_approach()

# the decision vector
opt.get_decision_vector()

# this function suggests physiological bounds
bounds, concentration_bounds = opt.get_optimization_bounds()

# overwriting specie concentration bounds for s4. Concentrations are in pM.
opt.get_concentration_bounds_species()

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds)

opt.generate_report()


# The reponse-related specie should be picked based on CellDesigner IDs. In our case phoshorylated A is s2.
# How to pick continuation parameter? In our case it is the amount of A protein, thus the conservation law 3.
print(opt.get_conservation_laws())


multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s6", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 'C2'},
                                                           print_lbls_flag=False, dir_path='./abstracted_cont_graphs')

# multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
#                                                            auto_parameters={'PrincipalContinuationParameter': 'C5'},
#                                                            print_lbls_flag=False, dir_path='./abstracted_cont_graphs')

opt.generate_report()
