
import re
import os
import warnings
import libsbml
import networkx

from .c_graph import Cgraph
from .low_deficiency_approach import LowDeficiencyApproach
from .mass_conservation_approach import MassConservationApproach
from .semi_diffusive_approach import SemiDiffusiveApproach
from .advanced_deficiency_approach import AdvancedDeficiencyApproach


class CRNT:
    """
    Class for managing CRNT methods.
    """
    def __init__(self, path):
        """
        Initialization of CRNT class.

        Parameters
        -----------
        path: string
            String representation of the path to the XML file.
        Example
        ---------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        """
        self.__path = path
        self.__sbml_model = None

        self.__validate_path()
        # self.__read_sbml_model()
        reader = libsbml.SBMLReader()
        sbml_doc = reader.readSBMLFromFile(self.__path)
        self.__sbml_model = sbml_doc.getModel()
        self.__cgraph = Cgraph(self.__sbml_model)  # init Cgraph object

    def __validate_path(self):
        if not os.path.isfile(self.__path):
            raise Exception("file is missing")

    def __read_sbml_model(self):
        reader = libsbml.SBMLReader()
        sbml_doc = reader.readSBMLFromFile(self.__path)
        # This is just to warn user. Not considered as a failure at this point.
        if sbml_doc.validateSBML() > 0:
            warnings.warn("Warning, SBML did not pass validation!\nProceeding anyway.")
        if sbml_doc.checkConsistency() > 0:
            warnings.warn("Warning, SBML did not pass checking of consistency!\nProceeding anyway.")
        self.__sbml_model = sbml_doc.getModel()
        pass

    # public methods

    def plot_c_graph(self):
        """
        Generates a matplotlib plot for the C-graph of the network using the networkx.draw function with circular and
        Kamada Kawai layout.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.plot_c_graph()
        """
        self.__cgraph.plot()

    def plot_save_c_graph(self):
        """
        Saves the matplotlib plot for the C-graph of the network using the networkx.draw function with circular and
        Kamada Kawai layout to the file network_cgraph.png

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.plot_save_c_graph()
        """
        self.__cgraph.plot_save()

    def print_c_graph(self):
        """
        Prints the reactions and reaction labels for the network.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> network.print_c_graph()
            Reaction graph of the form
            reaction -- reaction label:
            s1+s2 -> s3  --  re1
            s3 -> s1+s2  --  re1r
            s3 -> s6+s2  --  re2
            s6+s7 -> s16  --  re3
            s16 -> s6+s7  --  re3r
            s16 -> s7+s1  --  re4
            s1+s6 -> s15  --  re6
            s15 -> s1+s6  --  re6r
            s15 -> 2*s6  --  re8
        """
        self.__cgraph.print()

    def get_network_graphml(self):
        """
        Writes the NetworkX Digraph to the file network.graphml. Note that this generation only includes the names of
        the nodes, edges, and edge reaction names, it does not include other list attributes of the nodes and edges.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_network_graphml()
        """
        c_graph = self.__cgraph.get_graph()
        g = networkx.DiGraph()
        g.add_nodes_from(c_graph.nodes)
        g.add_edges_from(c_graph.edges)

        # add reaction names to graph
        edge_names = []
        for i in c_graph.edges:
            if (i[1], i[0]) in c_graph.edges:
                edge_names.append({"label": c_graph.edges[i]["label"], 'Bend': True})
            else:
                edge_names.append({"label": c_graph.edges[i]["label"], 'Bend': False})

        zip_obj = zip([i for i in g.edges], edge_names)
        edge_dict = dict(zip_obj)
        networkx.set_edge_attributes(g, edge_dict)

        #add if a complex is a boundary condition
        node_bc = [{'BC': c_graph.nodes[i]["species_bc"]} for i in c_graph.nodes]
        zip_obj = zip([i for i in g.nodes], node_bc)
        node_dict = dict(zip_obj)
        networkx.set_node_attributes(g, node_dict)

        networkx.write_graphml(g, "network.graphml")

    def basic_report(self):
        """
        Prints out basic CRNT properties of the network.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> network.basic_report()
            Number of species: 7
            Number of complexes: 9
            Number of reactions: 9
            Network deficiency: 2
        """

        report = """
        Number of species: %i
        Number of complexes: %i
        Number of reactions: %i 
        Network deficiency: %i
        """ % (
            len(self.__cgraph.get_species()),
            len(self.__cgraph.get_complexes()),
            len(self.__cgraph.get_reactions()),
            self.__cgraph.get_deficiency()
        )
        print("")
        print(re.sub(r"^\s+", "", report, flags=re.MULTILINE))

    def get_low_deficiency_approach(self):
        """
        Initializes and creates an object for the class LowDeficiencyApproach for the CRNT object constructed.

        See also
        ---------
        crnt4sbml.LowDeficiencyApproach

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> approach = network.get_low_deficiency_approach()
        """
        return LowDeficiencyApproach(self.__cgraph)

    def get_mass_conservation_approach(self):
        """
        Initializes and creates an object for the class MassConservationApproach for the CRNT object constructed.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        See also
        ---------
        crnt4sbml.MassConservationApproach

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 0.8010799999999998
            Solving for species' concentrations ...
            Elapsed time for finding species' concentrations: 1.1221319999999997
        """
        ldh = LowDeficiencyApproach(self.__cgraph)
        if ldh.does_satisfy_any_low_deficiency_theorem():
            message = """
            Network satisfies one of the low deficiency theorems.
            One should not run the optimization-based methods.
            """
            print(re.sub(r"^\s+", "", message, flags=re.MULTILINE))
            return None
        else:
            return MassConservationApproach(self.__cgraph)

    def get_semi_diffusive_approach(self):
        """
        Initializes and creates an object for the class SemiDiffusiveApproach for the CRNT object constructed.

        See also
        ---------
        crnt4sbml.SemiDiffusiveApproach

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/semi_diffusive_sbml_file.xml")
        >>> approach = network.get_semi_diffusive_approach()
        """
        ldh = LowDeficiencyApproach(self.__cgraph)
        if ldh.does_satisfy_any_low_deficiency_theorem():
            message = """
            Network satisfies one of the low deficiency theorems.
            One should not run the optimization-based methods.
            """
            print(re.sub(r"^\s+", "", message, flags=re.MULTILINE))
            return None
        else:
            return SemiDiffusiveApproach(self.__cgraph)

    def get_advanced_deficiency_approach(self):
        """
        Placeholder for Advanced Deficiency Approach. Future version of crnt4sbml will include the implementation of
        the Higher Deficiency Algorithm.
        """
        ldh = LowDeficiencyApproach(self.__cgraph)
        if ldh.does_satisfy_any_low_deficiency_theorem():
            message = """
            Network satisfies one of the low deficiency theorems.
            One should not  run the advanced deficiency method.
            """
            print(re.sub(r"^\s+", "", message, flags=re.MULTILINE))
            return None
        else:
            return AdvancedDeficiencyApproach(self.__cgraph)

    def get_c_graph(self):
        """
        Allows access to the class C-graph for the constructed CRNT object. Returns C-graph object for the provided
        CRNT object.

        See also
        ---------
        crnt4sbml.Cgraph

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> c_graph = network.get_c_graph()
        """
        return self.__cgraph

    @staticmethod
    def get_physiological_range(for_what=None):
        """Obtains physiological ranges.

        Parameters
        -----------
        for_what: string
            Accepted values: "concentration", "complex formation", "complex dissociation", or "catalysis"

        Returns
        --------
        concentration:
            5e-1,5e5 pM
        complex formation:
            1e-8,1e-4  pM^-1s^-1
        complex dissociation:
            1e-5,1e-3 s^-1
        catalysis:
            1e-3,1 s^-1

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_physiological_range("concentration")
        """
        valid = {"concentration", "complex_formation", "complex_dissociation", "monomolecular"}
        if for_what is None:
            warnings.warn("Please provide the argument what the range should be provided for.")
            return None
        if for_what not in valid:
            raise ValueError("for_what argument must be one of %r." % valid)
        if for_what == "concentration":
            return 5e-1, 5e5
        if for_what == "complex_formation":
            return 1e-8, 1e-4
        if for_what == "complex_dissociation":
            return 1e-5, 1e-3
        if for_what == "catalysis":
            return 1e-3, 1e0
