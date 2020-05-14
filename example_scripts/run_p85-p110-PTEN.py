import sys
sys.path.insert(0, "..")
import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/p85-p110-PTEN.xml")

approach = network.get_mass_conservation_approach()

bounds, concentration_bounds = approach.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds, iterations=5000,
                                                                          concentration_bounds=concentration_bounds,
                                                                          parallel_flag=True)

if approach.get_my_rank() == 0:

    network.basic_report()

    network.print_c_graph()

    ldt = network.get_low_deficiency_approach()
    ldt.report_deficiency_zero_theorem()
    ldt.report_deficiency_one_theorem()

    print("")

    print("Decision Vector:")
    print(approach.get_decision_vector())
    print("")

    print("Species for concentration bounds:")
    print(approach.get_concentration_bounds_species())

    multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s37", parameters=params_for_global_min,
                                                                                         auto_parameters={'PrincipalContinuationParameter': "C4"})

    approach.generate_report()
