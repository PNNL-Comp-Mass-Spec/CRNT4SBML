import networkx
import numpy
import sympy
import scipy.optimize
import sys
import matplotlib.cbook
import warnings
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


class Cgraph:
    """
    Class for constructing core CRNT values and C-graph of the network.
    """
    def __init__(self, model):
        """
        Initialization of Cgraph class.

        See also
        ---------
        crnt4sbml.CRNT.get_c_graph()
        """
        self.__g = networkx.DiGraph()

        # list of nodes and edges in the order they were added
        self.__g_nodes = []
        self.__g_edges = []

        self.__dic_species_id_bcs = {}
        for i in model.getListOfSpecies():
            self.__dic_species_id_bcs.update({i.getId(): i.getBoundaryCondition()})

        self.__species = []

        self.__parse_reactions(model)

        self.__assign_reaction_type()

        # self.__species = [i.getId() for i in model.getListOfSpecies() if not i.getBoundaryCondition()]

        # temporary makes testing of higher deficiency theory easier
        ######################################################################################
        #self.__species_name = [i.getName() for i in model.getListOfSpecies()]
        #self.__dic_species_id_name = {}
        #for i in range(len(self.__species)):
        #    self.__dic_species_id_name.update({self.__species[i]: self.__species_name[i]})
        ######################################################################################

        self.__linkage_classes = [self.__g.subgraph(c) for c in networkx.weakly_connected_components(self.__g)]

        self.__complexes = [self.__g.nodes[n]['label'] for n in self.__g_nodes]
        self.__reactions = [self.__g.edges[e]['label'] for e in self.__g_edges]

        # self.__check_weak_reversibility_of_linkage_classes() # todo: candidate for deletion

        # core calculations
        self.__create_y_matrix()
        self.__create_s_matrix()
        self.__create_a_matrix()
        self.__create_mass_action_monomials()

        if all([i <= 1 for i in self.get_number_of_terminal_strong_lc_per_lc()]) and \
            self.get_dim_equilibrium_manifold() > 0:
            self.__create_BT_matrix()
            self.__create_lambda_matrix()

        # self.__calculate_deficiency() # todo: candidate for deletion
        # self.__calculate_number_of_conservation_relations() # todo: candidate for deletion
        # self.__classify_dimensionality() # todo: candidate for deletion

    def __parse_reactants_and_products(self, reaction):

        # extracting reactants
        reactants = [i.getSpecies() for i in reaction.getListOfReactants()]

        bound_cond = any([self.__dic_species_id_bcs.get(n, n) for n in reactants])
        if bound_cond:
            color = 'red'
        else:
            color = 'green'

        stoichiometries = [int(i.getStoichiometry()) for i in reaction.getListOfReactants()]
        temp = zip([repr(i) if i != 1 else '' for i in stoichiometries], reactants)
        reactant_complex_name = '+'.join(['*'.join(filter(lambda x: x != '', i)) for i in temp])
        if self.__g.number_of_nodes() == 0 or reactant_complex_name not in [self.__g.nodes[v]['label'] for
                                                                            v in self.__g.nodes]:
            self.__g.add_node(reactant_complex_name,
                              label=reactant_complex_name,
                              species=reactants,
                              stoichiometries=stoichiometries,
                              species_bc=bound_cond,
                              color=color)
            self.__g_nodes.append(reactant_complex_name)

        # extracting products
        products = [i.getSpecies() for i in reaction.getListOfProducts()]

        bound_cond = any([self.__dic_species_id_bcs.get(n, n) for n in products])
        if bound_cond:
            color = 'red'
        else:
            color = 'green'

        stoichiometries = [int(i.getStoichiometry()) for i in reaction.getListOfProducts()]
        temp = zip([repr(i) if i != 1 else '' for i in stoichiometries], products)
        product_complex_name = '+'.join(['*'.join(filter(lambda x: x != '', i)) for i in temp])
        if product_complex_name not in [self.__g.nodes[v]['label'] for v in self.__g.nodes]: 
            self.__g.add_node(product_complex_name,
                              label=product_complex_name,
                              species=products,
                              stoichiometries=stoichiometries,
                              species_bc=bound_cond,
                              color=color)
            self.__g_nodes.append(product_complex_name)

        [self.__species.append(i) for i in reactants if (i not in self.__species) and (not self.__dic_species_id_bcs[i])]
        [self.__species.append(i) for i in products if (i not in self.__species) and (not self.__dic_species_id_bcs[i])]

        return [reactant_complex_name, product_complex_name]

    def __parse_reactants_products_and_modifiers(self, reaction):

        modifiers = [i.getSpecies() for i in reaction.getListOfModifiers()]

        # extracting reactants
        reactants = [i.getSpecies() for i in reaction.getListOfReactants()]

        reactant_complex_name = self.__create_reactant_or_product_complex(reaction, reactants, modifiers, "reactant")

        # extracting products
        products = [i.getSpecies() for i in reaction.getListOfProducts()]

        product_complex_name = self.__create_reactant_or_product_complex(reaction, products, modifiers, "product")

        # extracting modifiers
        bound_cond = any([self.__dic_species_id_bcs.get(n, n) for n in modifiers])
        if bound_cond:
            color = 'red'
        else:
            color = 'green'

        stoichiometries = [int(1.0)]
        modifiers_complex_name = "".join(modifiers) + "".join(reactants)

        if modifiers_complex_name not in [self.__g.nodes[v]['label'] for v in self.__g.nodes]:
            self.__g.add_node(modifiers_complex_name,
                              label=modifiers_complex_name,
                              species=[modifiers_complex_name],
                              stoichiometries=stoichiometries,
                              species_bc=bound_cond,
                              color=color)
            self.__g_nodes.append(modifiers_complex_name)

        [self.__species.append(i) for i in reactants if (i not in self.__species) and (not self.__dic_species_id_bcs[i])]
        [self.__species.append(i) for i in products if (i not in self.__species) and (not self.__dic_species_id_bcs[i])]
        [self.__species.append(i) for i in modifiers if (i not in self.__species) and (not self.__dic_species_id_bcs[i])]
        [self.__species.append(i) for i in [modifiers_complex_name] if (i not in self.__species)]

        return [reactant_complex_name, product_complex_name, modifiers_complex_name]

    def __create_reactant_or_product_complex(self, reaction, non_modifier, modifiers, comp):

        bound_cond = any([self.__dic_species_id_bcs.get(n, n) for n in non_modifier + modifiers])
        if bound_cond:
            color = 'red'
        else:
            color = 'green'

        # this statement assumes one species in modifier
        indices = [i for i in range(len(non_modifier)) if modifiers[0] == non_modifier[i]]
        if comp == 'reactant':

            reactants_stoich = [int(i.getStoichiometry()) for i in reaction.getListOfReactants()]

            # assumes only one modifier
            modifiers_sotich = int(1.0)

            if indices:
                stoich_sum = sum([reactants_stoich[i] for i in indices]) + modifiers_sotich
                stoichiometries = [reactants_stoich[i] for i in range(len(reactants_stoich)) if i not in indices] + \
                                  [stoich_sum]
                species_list = [non_modifier[i] for i in range(len(non_modifier)) if i not in indices] + modifiers

            else:
                stoichiometries = reactants_stoich + [modifiers_sotich]
                species_list = non_modifier + modifiers
        else:
            products_stoich = [int(i.getStoichiometry()) for i in reaction.getListOfProducts()]
            # assumes only one modifier
            modifiers_sotich = int(1.0)

            if indices:
                stoich_sum = sum([products_stoich[i] for i in indices]) + modifiers_sotich
                stoichiometries = [products_stoich[i] for i in range(len(products_stoich)) if i not in indices] + \
                                  [stoich_sum]
                species_list = [non_modifier[i] for i in range(len(non_modifier)) if i not in indices] + modifiers
            else:
                stoichiometries = products_stoich + [modifiers_sotich]
                species_list = non_modifier + modifiers

        temp = zip([repr(i) if i != 1 else '' for i in stoichiometries], species_list)
        complex_name = '+'.join(['*'.join(filter(lambda x: x != '', i)) for i in temp])

        if self.__g.number_of_nodes() == 0 or complex_name not in [self.__g.nodes[v]['label'] for
                                                                   v in self.__g.nodes]:
            self.__g.add_node(complex_name, label=complex_name, species=species_list, stoichiometries=stoichiometries,
                              species_bc=bound_cond, color=color)
            self.__g_nodes.append(complex_name)

        return complex_name

    def __extract_direct_reaction(self, reaction):
        r, p = self.__parse_reactants_and_products(reaction)
        if (r, p) not in [i for i in self.__g.edges]:
            self.__g.add_edge(r, p,
                              label=reaction.getId(),
                              k=None,
                              sbml_label=reaction.getId(),
                              type=None)
            self.__g_edges.append((r, p))

    def __extract_reverse_reaction(self, reaction):
        r, p = self.__parse_reactants_and_products(reaction)
        if (p, r) not in [i for i in self.__g.edges]:
            self.__g.add_edge(p, r,
                              label=reaction.getId() + 'r',
                              k=None,
                              sbml_label=reaction.getId(),
                              type=None)
            self.__g_edges.append((p, r))

    def __extract_catalysis_reaction(self, reaction):

        modifier = self.__get_celldesigner_reaction_type(reaction)

        if modifier != 'CATALYSIS':

            exception_statement = "For reaction " + reaction.getId() + " there is a modifier that cannot be " \
                                                                       "identified, there are multiple modifiers," \
                                                                       " or is not catalysis. \nThis may be" \
                                                                       " a result of not constructing the SBML file " \
                                                                       "using CellDesigner. \nPlease make the appropriate" \
                                                                       " changes to this reaction to continue."
            raise Exception(exception_statement)

        r, p, m = self.__parse_reactants_products_and_modifiers(reaction)

        # add complex formation
        if (r, m) not in [i for i in self.__g.edges]:
            self.__g.add_edge(r, m,
                              label=reaction.getId() + 'f',
                              k=None,
                              sbml_label=reaction.getId(),
                              type="complex formation")
            self.__g_edges.append((r, m))

        # add complex dissociation
        if (m, r) not in [i for i in self.__g.edges]:
            self.__g.add_edge(m, r,
                              label=reaction.getId() + 'd',
                              k=None,
                              sbml_label=reaction.getId(),
                              type="complex dissociation")
            self.__g_edges.append((m, r))

        # add catalysis
        if (m, p) not in [i for i in self.__g.edges]:
            self.__g.add_edge(m, p,
                              label=reaction.getId() + 'c',
                              k=None,
                              sbml_label=reaction.getId(),
                              type="catalysis")
            self.__g_edges.append((m, p))

    def __parse_reactions(self, model):
        for i in model.getListOfReactions():
            if i.getNumModifiers() != 0:
                self.__extract_catalysis_reaction(i)
            else:
                self.__extract_direct_reaction(i)
                if i.getReversible():
                    self.__extract_reverse_reaction(i)

    def __get_celldesigner_reaction_type(self, i):

        i_initial = i.getAnnotation().getChild(0)

        child_number = self.__get_child_number(i_initial, 'reactionType')

        if child_number is not None:

            if i.getNumModifiers() == 1:

                child_number1 = self.__get_child_number(i_initial, 'listOfModification')

                if child_number1 is not None:

                    child_number2 = self.__get_child_number(i_initial.getChild(child_number1), 'modification')

                    ii = i_initial.getChild(child_number1).getChild(child_number2)
                    modification = self.__get_modification_type(ii)

                    return modification
                else:
                    return None
            else:
                return None

        else:
            return None

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

    def __assign_reaction_type(self):

        for i in self.get_g_edges():

            current_edge = self.__g.edges[i]

            if current_edge['type'] is None:

                reactant_node = self.__g.nodes[i[0]]
                product_node = self.__g.nodes[i[1]]

                if sum(reactant_node['stoichiometries']) == 2 and sum(product_node['stoichiometries']) == 1:
                    current_edge['type'] = "complex formation"

                elif sum(reactant_node['stoichiometries']) == 1 and sum(product_node['stoichiometries']) == 2:
                    if (i[1], i[0]) in self.get_g_edges():
                        current_edge['type'] = "complex dissociation"
                    else:
                        current_edge['type'] = "catalysis"

                elif sum(reactant_node['stoichiometries']) == 1 and sum(product_node['stoichiometries']) == 1:
                    current_edge['type'] = "catalysis"

                elif sum(reactant_node['stoichiometries']) == 1 and sum(product_node['stoichiometries']) == 3:
                    current_edge['type'] = "catalysis"

                else:
                    output_statement = "The reaction type of reaction " + current_edge['label'] + " could not be " +\
                                       "identified as it does not fit any biological criteria established. \n" +\
                                       "Automatic creation of bounds for optimization routines cannot be completed. \n"

                    print(output_statement)

    # core CRNT calcs
    def __create_y_matrix(self):

        self.__Y = sympy.zeros(len(self.__species), len(self.__complexes))
        count = 0
        for i in self.__g_nodes:
            v = self.__g.nodes[i]
            stoich = dict(zip(v['species'], v['stoichiometries']))
            self.__Y[:, count] = [stoich[i] if i in stoich.keys() else 0 for i in self.__species]
            count += 1

    def __create_s_matrix(self):

        self.__S = sympy.zeros(len(self.__species), len(self.__reactions))
        count = 0 
        for i in self.__g_edges:
            target = self.__complexes.index(i[1])
            source = self.__complexes.index(i[0])
            self.__S[:, count] = self.__Y[:, target] - self.__Y[:, source]
            count += 1

    def __create_lambda_matrix(self):
        self.__Lambda = sympy.zeros(len(self.__complexes), len(self.__linkage_classes))
        for li in range(len(self.__linkage_classes)):
            l_complexes = [i for i in self.__linkage_classes[li].nodes]
            self.__Lambda[:, li] = [1 if i in l_complexes else 0 for i in self.__complexes]

    def __create_a_matrix(self):
        self.__A = sympy.zeros(len(self.__complexes), len(self.__complexes))
        inc_list = [list(self.__g.out_edges(n, data='label')) for n in self.__g_nodes]

        for i in range(len(self.__complexes)):
            for ii in inc_list[i]:
                self.__A[i, i] = self.__A[i, i] - sympy.Symbol(ii[2], positive=True)

        for i in range(len(self.__complexes)):
            for j in range(len(self.__complexes)):
                if self.__g.has_edge(self.__complexes[j], self.__complexes[i]):
                    self.__A[i, j] = self.__A[i, j] + sympy.Symbol(self.__g.get_edge_data(self.__complexes[j],
                                                                                          self.__complexes[i])['label'],
                                                                   positive=True)

    def __create_BT_matrix(self):
        # (bt = NullSpace[Transpose[Y.A]]) // MatrixForm
        the_null_space = (self.__Y * self.__A).T.nullspace()

        # getting the number of columns given in nullspace computation
        sizes = len(the_null_space)

        # simplifying the entries given in the columns
        # produced by the nullspace
        for i in range(0, sizes):
            the_null_space[i] = sympy.simplify(the_null_space[i])

        bt_temp = self.__create_non_negative_b_matrix(the_null_space, sizes)
        
        # taking out zero rows if they are in unique vectors
        bt = sympy.Matrix(bt_temp[~numpy.all(bt_temp == 0, axis=1)]).evalf(10)

        # making the dimension of bt be lambda by N, just incase
        # we have more vectors than we should
        bt = bt[0:self.get_dim_equilibrium_manifold(), :]

        # putting in a check to make sure that the rows
        # of bt are linearly independent
        _, vals = bt.T.rref()
        if len(vals) != self.get_dim_equilibrium_manifold():
            print(" ")
            print("Error: Not all rows of B^T are linearly independent!")
            print(" ")
            sys.exit()

        self.__B = bt

    def __simplify_nullspace(self, the_null_space):

        # Flag that says a species or reaction was found in one of the column vectors of the nullspace
        denom_flag = False

        print("beginning null space")
        #print(the_null_space)
        print("")
        for i in range(0, len(the_null_space)):
            #print(the_null_space[i])
            the_null_space[i] = sympy.simplify(the_null_space[i])

        print(the_null_space)

        # find those columns that have a common denominator
        denom_list = []
        for ii in the_null_space:

            # getting common denominator
            denom_in_column = [sympy.denom(i) for i in ii]
            denom_list.append(denom_in_column)

        print(f"denom_list {denom_list}")

        like_columns = []
        for i in denom_list:
            same_denoms = [ii for ii in range(len(denom_list)) if denom_list[ii] == i]
            if same_denoms not in like_columns:
                like_columns.append(same_denoms)


        print(like_columns)
        temp_null = [sympy.zeros(len(self.__species), 1)] * (len(self.__species) - self.__S.rank())
        count = 0
        for ii in like_columns:

            for i in ii:
                temp_null[count] += the_null_space[i]
                temp_null[count] = sympy.simplify(temp_null[count])
            count += 1

        print(temp_null)
        sys.exit()

        # collecting nullspace vectors with common denominators
        all_denom = list(set(denom_list))
        temp_null = [sympy.zeros(len(self.__species), 1)] * (len(self.__species) - self.__S.rank())


        try:
            for ii in range(len(all_denom)):
                indices = [i for i, x in enumerate(denom_list) if x == all_denom[ii]]
                print(f"indices {indices}")
                for i in indices:
                    temp_null[ii] += the_null_space[i]
                    temp_null[ii] = sympy.simplify(temp_null[ii])
        except Exception as e:
            print("hi")

        print(f"final temp null {temp_null}")
        sys.exit()

        # checking to see if a reaction or species is in the reduced nullspace
        for i in temp_null:
            atoms_of_col_vec = [ii.atoms() for ii in i[:, 0]]
            for ii in atoms_of_col_vec:
                if any([iii in self.__reactions + self.__species for iii in ii]):
                    denom_flag = True

        if denom_flag:
            raise Exception("Nullspace calculation from S contains SymPy variables and this could not be resolved.")

        return temp_null

    def __create_non_negative_b_matrix(self, the_null_space, sizes):

        a_null = numpy.zeros((len(self.__species), sizes))
        for i in range(len(self.__species)):
            temp_vec = []
            for j in range(sizes):
                temp_vec.append(the_null_space[j][i])
            a_null[i, :] = temp_vec

        a_eq = numpy.array([numpy.sum(a_null, axis=0)])
        # must multiply by negative one because in optimization we have the inequality <= 0.0
        a_ub = -1.0*a_null
        b_ub = numpy.zeros(len(self.__species))
        b_eq = numpy.array([1.0])        
        
        # defining the number of solutions to simulate
        num_sols = ((a_ub.shape[0]+1)*(a_ub.shape[1]+1))*10

        # a matrix to hold the different solutions of the method
        sol = numpy.zeros((num_sols, len(self.__species)))

        # searching the nullspace for nonnegative vectors
        for i in range(num_sols):
            # generating a multivariate normal random variates
            crit = numpy.random.normal(0, 1.0, a_ub.shape[1])

            # Solving the linear programming problem
            out = scipy.optimize.linprog(crit, A_eq=a_eq, b_eq=b_eq, A_ub=a_ub, b_ub=b_ub, bounds=(None, None),
                                         method='simplex')
            # multiplying our nullspace vectors with the minimized value
            # to create a nonnegative vector
            sol[i, :] = numpy.dot(a_null, out.x)
            sol[i, [numpy.abs(sol[i, :]) < numpy.finfo(float).eps][0]] = numpy.float64(0.0)

        # getting the unique vectors
        unique_vecs = numpy.unique(sol, axis=0)
        num_rows = numpy.size(unique_vecs, 0)

        # taking the smallest nonzero entry of the unique vectors
        # and dividing by it this will hopfully create nice
        # looking vectors that are whole numbers
        for i in range(num_rows):
            minval = numpy.min(unique_vecs[i, numpy.nonzero(unique_vecs[i, :])])
            unique_vecs[i, :] = unique_vecs[i, :]/numpy.float64(minval)
            
        unique_vecs = numpy.unique(unique_vecs.round(decimals=10), axis=0)

        return unique_vecs

    def __create_mass_action_monomials(self):
        self.__psi = sympy.ones(len(self.__complexes), 1)
        for j in range(len(self.__complexes)):
            for i in range(len(self.__species)):
                self.__psi[j] *= sympy.Symbol(self.__species[i], positive=True)**self.__Y[i, j]
                
                # TODO: possibly make concetrations positive

    # public getters
    def get_dict_id_name(self):
        return self.__dic_species_id_name

    def get_ode_system(self):
        """
        Returns SymPy matrix representing the ODE system.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the
        provided example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> sympy.pprint(network.get_c_graph().get_ode_system())
            ⎡    -re₁⋅s₁⋅s₂ + re1r⋅s₃ + re₄⋅s₁₆ - re₅⋅s₁⋅s₆ + re5r⋅s₁₅     ⎤
            ⎢                                                              ⎥
            ⎢                 -re₁⋅s₁⋅s₂ + s₃⋅(re1r + re₂)                 ⎥
            ⎢                                                              ⎥
            ⎢                 re₁⋅s₁⋅s₂ + s₃⋅(-re1r - re₂)                 ⎥
            ⎢                                                              ⎥
            ⎢re₂⋅s₃ - re₃⋅s₆⋅s₇ + re3r⋅s₁₆ - re₅⋅s₁⋅s₆ + s₁₅⋅(re5r + 2⋅re₆)⎥
            ⎢                                                              ⎥
            ⎢                -re₃⋅s₆⋅s₇ + s₁₆⋅(re3r + re₄)                 ⎥
            ⎢                                                              ⎥
            ⎢                re₃⋅s₆⋅s₇ + s₁₆⋅(-re3r - re₄)                 ⎥
            ⎢                                                              ⎥
            ⎣                re₅⋅s₁⋅s₆ + s₁₅⋅(-re5r - re₆)                 ⎦
        """
        return self.__Y*self.__A*self.__psi

    def get_graph(self):
        """
        Returns the NetworkX DiGraph representation of the network.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_c_graph().get_graph()
        """
        return self.__g

    def get_species(self):
        """
        Returns Python list of strings representing the species of the network.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> print(network.get_c_graph().get_species())
            ['s1', 's2', 's3', 's6', 's7', 's16', 's15']
        """
        return self.__species

    def get_complexes(self):
        """
        Returns Python list of strings representing the complexes of the network.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> print(network.get_c_graph().get_complexes())
            ['s1+s2', 's3', 's6+s2', 's6+s7', 's16', 's7+s1', 's1+s6', 's15', '2*s6']
        """
        return self.__complexes

    def get_reactions(self):
        """
        Returns Python list of strings representing the reactions of the network.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> print(network.get_c_graph().get_reactions())
            ['re1', 're1r', 're2', 're3', 're3r', 're4', 're5', 're5r', 're6']
        """
        return self.__reactions

    def get_a(self):
        """
        Returns SymPy matrix representing the kinetic constant matrix, :math:`A`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> sympy.pprint(network.get_c_graph().get_a())
            ⎡-re₁     re1r      0   0         0       0   0         0       0⎤
            ⎢                                                                ⎥
            ⎢re₁   -re1r - re₂  0   0         0       0   0         0       0⎥
            ⎢                                                                ⎥
            ⎢ 0        re₂      0   0         0       0   0         0       0⎥
            ⎢                                                                ⎥
            ⎢ 0         0       0  -re₃     re3r      0   0         0       0⎥
            ⎢                                                                ⎥
            ⎢ 0         0       0  re₃   -re3r - re₄  0   0         0       0⎥
            ⎢                                                                ⎥
            ⎢ 0         0       0   0        re₄      0   0         0       0⎥
            ⎢                                                                ⎥
            ⎢ 0         0       0   0         0       0  -re₅     re5r      0⎥
            ⎢                                                                ⎥
            ⎢ 0         0       0   0         0       0  re₅   -re5r - re₆  0⎥
            ⎢                                                                ⎥
            ⎣ 0         0       0   0         0       0   0        re₆      0⎦
        """
        return self.__A

    def get_y(self):
        """
        Returns SymPy matrix representing the molecularity matrix, :math:`Y`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> sympy.pprint(network.get_c_graph().get_y())
            ⎡1  0  0  0  0  1  1  0  0⎤
            ⎢                         ⎥
            ⎢1  0  1  0  0  0  0  0  0⎥
            ⎢                         ⎥
            ⎢0  1  0  0  0  0  0  0  0⎥
            ⎢                         ⎥
            ⎢0  0  1  1  0  0  1  0  2⎥
            ⎢                         ⎥
            ⎢0  0  0  1  0  1  0  0  0⎥
            ⎢                         ⎥
            ⎢0  0  0  0  1  0  0  0  0⎥
            ⎢                         ⎥
            ⎣0  0  0  0  0  0  0  1  0⎦
        """
        return self.__Y

    def get_s(self):
        """
        Returns SymPy matrix representing the stoichiometric matrix, :math:`S`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> sympy.pprint(network.get_c_graph().get_s())
            ⎡-1  1   0   0   0   1   -1  1   0 ⎤
            ⎢                                  ⎥
            ⎢-1  1   1   0   0   0   0   0   0 ⎥
            ⎢                                  ⎥
            ⎢1   -1  -1  0   0   0   0   0   0 ⎥
            ⎢                                  ⎥
            ⎢0   0   1   -1  1   0   -1  1   2 ⎥
            ⎢                                  ⎥
            ⎢0   0   0   -1  1   1   0   0   0 ⎥
            ⎢                                  ⎥
            ⎢0   0   0   1   -1  -1  0   0   0 ⎥
            ⎢                                  ⎥
            ⎣0   0   0   0   0   0   1   -1  -1⎦
        """
        return self.__S

    def get_b(self):
        """
        Returns SymPy matrix representing the mass conservation matrix, :math:`B`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> sympy.pprint(network.get_c_graph().get_b())
            ⎡ 0    0    0    0   1.0  1.0   0 ⎤
            ⎢                                 ⎥
            ⎢ 0   1.0  1.0   0    0    0    0 ⎥
            ⎢                                 ⎥
            ⎣1.0   0   1.0  1.0   0   1.0  2.0⎦
        """
        return self.__B

    def get_lambda(self):
        """
        Returns SymPy matrix representing the linkage class matrix, :math:`\Lambda`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> sympy.pprint(network.get_c_graph().get_lambda())
            ⎡1  0  0⎤
            ⎢       ⎥
            ⎢1  0  0⎥
            ⎢       ⎥
            ⎢1  0  0⎥
            ⎢       ⎥
            ⎢0  1  0⎥
            ⎢       ⎥
            ⎢0  1  0⎥
            ⎢       ⎥
            ⎢0  1  0⎥
            ⎢       ⎥
            ⎢0  0  1⎥
            ⎢       ⎥
            ⎢0  0  1⎥
            ⎢       ⎥
            ⎣0  0  1⎦
        """
        return self.__Lambda

    def get_psi(self):
        """
        Returns SymPy matrix representing the mass action monomials, :math:`\psi`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> sympy.pprint(network.get_c_graph().get_psi())
            ⎡s₁⋅s₂⎤
            ⎢     ⎥
            ⎢ s₃  ⎥
            ⎢     ⎥
            ⎢s₂⋅s₆⎥
            ⎢     ⎥
            ⎢s₆⋅s₇⎥
            ⎢     ⎥
            ⎢ s₁₆ ⎥
            ⎢     ⎥
            ⎢s₁⋅s₇⎥
            ⎢     ⎥
            ⎢s₁⋅s₆⎥
            ⎢     ⎥
            ⎢ s₁₅ ⎥
            ⎢     ⎥
            ⎢   2 ⎥
            ⎣ s₆  ⎦
        """
        return self.__psi

    def get_g_nodes(self):
        """
        Returns a list of strings that represent the order of the nodes of the NetworkX DiGraph.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_c_graph().get_g_nodes()
        """
        return self.__g_nodes

    def get_g_edges(self):
        """
        Returns a list of tuples of strings that represent the order of the edges of the NetworkX DiGraph.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_c_graph().get_g_edges()
        """
        return self.__g_edges

    def get_weak_reversibility_of_linkage_classes(self):
        """
        Returns list of Python boolean types for the weak reversibility of each linkage class. If the linkage class is
        weakly reversible then the entry in the list is True, False otherwise with order as defined by
        :func:`crnt4sbml.Cgraph.get_linkage_classes`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_c_graph().get_weak_reversibility_of_linkage_classes()
        """
        return [networkx.is_strongly_connected(i) for i in self.__linkage_classes]

    def get_if_cgraph_weakly_reversible(self):
        """
        Returns weak reversibility of the network. If the network is weakly reversible True is returned, False otherwise.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_c_graph().get_if_cgraph_weakly_reversible()
        """
        return all(self.get_weak_reversibility_of_linkage_classes())

    def get_dim_equilibrium_manifold(self):
        """
        Returns integer value representing the dimension of the equilibrium manifold, :math:`\lambda`. This value
        is the number of mass conservation relationships.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> print(network.get_c_graph().get_dim_equilibrium_manifold())
            3
        """
        return len(self.__species) - self.__S.rank()

    def get_deficiency(self):
        """
        Returns integer value representing the deficiency of the network, :math:`\delta`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> print(network.get_c_graph().get_deficiency())
            2
        """
        return len(self.__complexes) - len(self.__linkage_classes) - self.__S.rank()

    def get_linkage_classes(self):
        """
        Returns list of NetworkX subgraphs representing the linkage classes.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_c_graph().get_linkage_classes()
        """
        return self.__linkage_classes

    def get_network_dimensionality_classification(self):
        """
        Returns a two element list specifying the dimensionality of the network.
        Possible output:
        ["over-dimensioned",0]

        or

        ["proper",1]

        or

        ["under-dimensioned",2]

        or

        ["NOT DEFINED!",3]

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_c_graph().get_network_dimensionality_classification()
        """
        temp = self.get_dim_equilibrium_manifold() - self.get_deficiency()
        if temp < 0:
            network_class = "over-dimensioned"
            classification = 0
        elif temp == 0:
            network_class = "proper"
            classification = 1
        elif temp > 0:
            network_class = "under-dimensioned"
            classification = 2
        else:
            network_class = "NOT DEFINED!"
            classification = 3
        return [network_class, classification]

    def get_linkage_classes_deficiencies(self):
        """
        Returns an interger list of each linkage class deficiency. Here, the first element corresponds to the first
        linkage class with order as defined by :func:`crnt4sbml.Cgraph.get_linkage_classes`.

        Example
        ---------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_c_graph().get_linkage_classes_deficiencies()
        """
        deficiencies = []
        for l in self.__linkage_classes:
            temp = [self.__g.edges[e]['label'] for e in l.edges]
            reactions_idx = [self.__reactions.index(i) for i in temp]
            deficiencies.append(networkx.number_of_nodes(l) - 1 - self.__S[:, reactions_idx].rank())
        return deficiencies

    def get_number_of_terminal_strong_lc_per_lc(self):
        """
        Returns an integer list stating the number of terminally strong linkage classes per linkage class. Here,
        the first element corresponds to the first linkage class with order as defined by
        :func:`crnt4sbml.Cgraph.get_linkage_classes`.

        Example
        ---------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> network.get_c_graph().get_number_of_terminal_strong_lc_per_lc()
        """
        number_of_terminal_strong_lc_per_lc = []
        for lc_i in self.__linkage_classes:
            num_term_strong = 0
            strong_lcs = [lc_i.subgraph(c).copy() for c in networkx.strongly_connected_components(lc_i)]
            for slc_j in strong_lcs:
                slc_j_lbls = [i for i in slc_j.nodes]
                x = [lc_i.successors(v_i) for v_i in slc_j_lbls]
                connected_out_vertices = [item for sublist in x for item in sublist]
                d = set(connected_out_vertices) - set(slc_j_lbls)
                if d == set():
                    num_term_strong += 1
            number_of_terminal_strong_lc_per_lc.append(num_term_strong)
        return number_of_terminal_strong_lc_per_lc

    # public methods
    def print(self):
        """
        Prints edges and nodes of NetworkX DiGraph.

        See also
        ---------
        crnt4sbml.CRNT.print_c_graph
        """
        print("")
        print("Reaction graph of the form") 
        print("reaction -- reaction label:")
        for e in self.__g_edges:
            print(e[0] + ' -> ' + e[1] + '  --  ' + self.__g.edges[e]['label'])
        print("")

    def plot(self):
        """
        Plots NetworkX DiGraph.

        See also
        ---------
        crnt4sbml.CRNT.plot_c_graph
        """
        pos = networkx.circular_layout(self.__g, scale=1.5)
        pos = networkx.kamada_kawai_layout(self.__g, pos=pos, scale=1.5)
        node_colors = [self.__g.nodes[n]['color'] for n in self.__g_nodes]
        networkx.draw(self.__g, pos, node_color=node_colors)
        node_labels = networkx.get_node_attributes(self.__g, 'label')
        networkx.draw_networkx_labels(self.__g, pos, labels=node_labels)
        edge_labels = networkx.get_edge_attributes(self.__g, 'label')
        networkx.draw_networkx_edge_labels(self.__g, pos, edge_labels, label_pos=0.3)
        plt.show()

    def plot_save(self):
        """
        Saves the plot of the NetworkX DiGraph.

        See also
        ---------
        crnt4sbml.CRNT.plot_save_c_graph
        """
        pos = networkx.circular_layout(self.__g, scale=1.5)
        pos = networkx.kamada_kawai_layout(self.__g, pos=pos, scale=1.5)
        node_colors = [self.__g.nodes[n]['color'] for n in self.__g_nodes]
        networkx.draw(self.__g, pos, node_color=node_colors)
        node_labels = networkx.get_node_attributes(self.__g, 'label')
        networkx.draw_networkx_labels(self.__g, pos, labels=node_labels)
        edge_labels = networkx.get_edge_attributes(self.__g, 'label')
        networkx.draw_networkx_edge_labels(self.__g, pos, edge_labels, label_pos=0.3)
        plt.savefig('network_cgraph.png')
