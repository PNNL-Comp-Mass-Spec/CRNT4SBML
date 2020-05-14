# CRNT class
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
from .general_approach import GeneralApproach


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

    def get_celldesigner_reaction_names(self):

        reaction_type_and_modification = []

        for i in self.__sbml_model.getListOfReactions():

            i_initial = i.getAnnotation().getChild(0)

            child_number = self.__get_child_number(i_initial, 'reactionType')

            if child_number is not None:

                reaction_type = i_initial.getChild(child_number).getChild(0).toXMLString()

                if i.getNumModifiers() != 0:

                    child_number1 = self.__get_child_number(i_initial, 'listOfModification')

                    if child_number1 is not None:

                        child_number2 = self.__get_child_number(i_initial.getChild(child_number1), 'modification')

                        ii = i_initial.getChild(child_number1).getChild(child_number2)
                        modification = self.__get_modification_type(ii)

                        reaction_type_and_modification.append((reaction_type, modification))
                    else:
                        reaction_type_and_modification.append((reaction_type, None))
                else:
                    reaction_type_and_modification.append((reaction_type, None))

            else:
                reaction_type_and_modification.append((None, None))

            reactants = [ii.getSpecies() for ii in i.getListOfReactants()]
            products = [ii.getSpecies() for ii in i.getListOfProducts()]
            modifiers = [ii.getSpecies() for ii in i.getListOfModifiers()]
            print([i.getModifier(ii) for ii in modifiers])
            print("reactants")
            print(reactants)
            print("products")
            print(products)
            print("modifiers")
            print(modifiers)


        is_any_none = any([i[0] is None for i in reaction_type_and_modification])

        if is_any_none:

            print("Some reaction types cannot be found, CellDesigner may not have been used to generate SBML file.")

        return reaction_type_and_modification

    def __get_child_number(self, i, string_value):

        num_children = i.getNumChildren()
        children_strings = [i.getChild(ii).getName() for ii in range(num_children)]
        if string_value in children_strings:
            child_number = children_strings.index(string_value)
        else:
            child_number = None

        return child_number

    def __get_modification_type(self, ii):

        attribute_names = [ii.getAttrName(iii) for iii in range(ii.getAttributesLength())]

        if 'type' in attribute_names:

            indx_type = attribute_names.index('type')
            attribute_values = [ii.getAttrValue(iii) for iii in attribute_names]

            return attribute_values[indx_type]
        else:
            return None

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
            s1+s6 -> s15  --  re5
            s15 -> s1+s6  --  re5r
            s15 -> 2*s6  --  re6
        """
        self.__cgraph.print()

    def print_biological_reaction_types(self):
        """
        Prints the reactions, reaction labels, and biological reaction type for the network.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> network.print_biological_reaction_types()
            Reaction graph of the form
            reaction -- reaction label -- biological reaction type:
            s1+s2 -> s3  --  re1 -- complex formation
            s3 -> s1+s2  --  re1r -- complex dissociation
            s3 -> s6+s2  --  re2 -- catalysis
            s6+s7 -> s16  --  re3 -- complex formation
            s16 -> s6+s7  --  re3r -- complex dissociation
            s16 -> s7+s1  --  re4 -- catalysis
            s1+s6 -> s15  --  re5 -- complex formation
            s15 -> s1+s6  --  re5r -- complex dissociation
            s15 -> 2*s6  --  re6 -- catalysis
        """

        print("")
        print("Reaction graph of the form")
        print("reaction -- reaction label -- biological reaction type:")
        for e in self.__cgraph.get_g_edges():
            print(e[0] + ' -> ' + e[1] + '  --  ' + self.__cgraph.get_graph().edges[e]['label'] + ' -- ' + \
                  self.__cgraph.get_graph().edges[e]['type'])
        print("")


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
            Elapsed time for creating Equilibrium Manifold: 2.060944
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
            return MassConservationApproach(self.__cgraph, self.get_physiological_range)

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
            return SemiDiffusiveApproach(self.__cgraph, self.get_physiological_range)

    def get_general_approach(self):
        """
        Initializes and creates an object for the class GeneralApproach for the CRNT object constructed.

        See also
        ---------
        crnt4sbml.GeneralApproach

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        """
        ldh = LowDeficiencyApproach(self.__cgraph)
        if ldh.does_satisfy_any_low_deficiency_theorem():
            message = """
            Network satisfies one of the low deficiency theorems.
            One should not run the optimization-based methods.
            """
            print(re.sub(r"^\s+", "", message, flags=re.MULTILINE))
            return GeneralApproach(self.__cgraph, self.get_physiological_range)
        else:
            return GeneralApproach(self.__cgraph, self.get_physiological_range)

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
            Accepted values: "concentration", "complex formation", "complex dissociation", "catalysis", or "flux"

        Returns
        --------
        concentration: tuple
            (5e-1,5e5) pM
        complex formation: tuple
            (1e-8,1e-4)  pM^-1s^-1
        complex dissociation: tuple
            (1e-5,1e-3) s^-1
        catalysis: tuple
            (1e-3,1) s^-1
        flux: tuple
            (0, 55) M s^-1


        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_physiological_range("concentration")
        """
        valid = {"concentration", "complex formation", "complex dissociation", "catalysis", "flux"}
        if for_what is None:
            warnings.warn("Please provide the argument what the range should be provided for.")
            return None
        if for_what not in valid:
            raise ValueError("for_what argument must be one of %r." % valid)
        if for_what == "concentration":
            return (5e-1, 5e5)
        if for_what == "complex formation":
            return (1e-8, 1e-4)
        if for_what == "complex dissociation":
            return (1e-5, 1e-3)
        if for_what == "catalysis":
            return (1e-3, 1e0)
        if for_what == "flux":
            return (0, 55)
