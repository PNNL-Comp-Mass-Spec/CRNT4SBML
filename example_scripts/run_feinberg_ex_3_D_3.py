import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/feinberg_ex_3_D_3.xml")

network.basic_report()

network.print_c_graph()

ldt = network.get_low_deficiency_approach()
ldt.report_deficiency_zero_theorem()
ldt.report_deficiency_one_theorem()

approach = network.get_mass_conservation_approach()
