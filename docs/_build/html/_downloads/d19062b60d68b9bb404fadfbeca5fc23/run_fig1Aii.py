import crnt4sbml_test

network = crnt4sbml_test.CRNT("../sbml_files/Fig_1Aii.xml")

network.basic_report()

network.print_c_graph()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

opt = network.get_mass_conservation_approach()

