import sys
sys.path.insert(0, "..")

import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/open_fig5B.xml")

network.get_network_graphml()

