import sys
sys.path.insert(0, "..")

import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/HDA_3_2_1.xml")

approach = network.get_advanced_deficiency_approach()

approach.run_higher_deficiency_algorithm()

