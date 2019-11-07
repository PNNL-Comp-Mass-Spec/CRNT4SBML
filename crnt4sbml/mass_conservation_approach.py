import os
import numpy
import sympy
import sympy.utilities.lambdify
import sys 
import time
import numpy.linalg
import itertools
import warnings
import math
from .bistability_finder import BistabilityFinder


# TODO: see why these imports produce
# segmentation fault: 11
##########################
# import pylab as plt
# import tellurium as te
# from rrplugins import
##########################

class MassConservationApproach:
    """
    Class for constructing variables and methods needed for the mass conservation approach.
    """
    def __init__(self, cgraph, get_physiological_range):
        """
        Initialization of the MassConservationApproach class.

        See also
        ---------
        crnt4sbml.CRNT.get_mass_conservation_approach()
        """
        self.__cgraph = cgraph

        self.get_physiological_range = get_physiological_range

        if not all([i <= 1 for i in self.__cgraph.get_number_of_terminal_strong_lc_per_lc()]):
            print("The network is not uniterminal!")
            sys.exit()

        if not self.__cgraph.get_dim_equilibrium_manifold() > 0:
            print("# of species - rank(S) is not greater than zero!")
            print("The mass conservation approach cannot be ran!")
            sys.exit()

        # declare key fields
        self.__deficiency_pars = None
        self.__concentration_pars = None
        self.__reaction_pars = None 
        self.__W = None  # matrix
        self.__W_nullspace = None
        self.__H = None  # vector
        self.__G = None  # matrix
        self.__symbolic_objective_fun = None
        self.__concentration_vals = None
        self.__decision_vector_x = None
        self.__concentration_funs = None 
        self.__objective_fun_params = None
        self.__lambda_objective_fun = None 
        self.__important_info = ""
        self.__numpy_dtype = None
        self.__independent_odes = None
        self.__independent_species = None

        # vars used frequently
        self.__N = len(self.__cgraph.get_species())
        self.__R = len(self.__cgraph.get_reactions())
        self.__species = self.__cgraph.get_species()
        self.__reactions = self.__cgraph.get_reactions()
        self.__delta = self.__cgraph.get_deficiency()
        self.__M = len(self.__cgraph.get_complexes())
        self.__ell = len(self.__cgraph.get_linkage_classes())
        self.__lambda = self.__cgraph.get_dim_equilibrium_manifold()
        self.__classification = self.__cgraph.get_network_dimensionality_classification()
        # compute necessary vectors and matrices
        self.__create_deficiency_pars()
        self.__create_concentration_pars()
        self.__create_reaction_pars()
        self.__create_w_matrix()
        self.__create_w_nullspace()
        self.__create_h_vector()
        self.__create_g_matrix()
        self.__create_symbolic_objective_fun()
        self.__create_decision_vector_x()
        self.__create_concentration_bounds_species()
        self.__create_concentration_lambda_fun()
        self.__create_objective_fun__lambda_fun()
        self.__create_g_matrix_lambda_fun()
        self.__create_dch_matrix_lambda_fun()

    def __create_deficiency_pars(self):
        # creating a vector of the deficiency parameters
        # \alpha_1, ..., \alpha_\delta
        self.__deficiency_pars = [sympy.symbols('a' + str(i + 1), real=True) for i
                                  in range(self.__delta)]

    def __create_concentration_pars(self):
        # putting the species in a list for more readable code
        self.__concentration_pars = [sympy.Symbol(self.__species[i], positive=True) for i in range(self.__N)]

    def __create_reaction_pars(self):
        self.__reaction_pars = [sympy.Symbol(self.__reactions[i], positive=True) for i in range(self.__R)]
        
    def __create_w_matrix(self):
        # concatenating Y and Lambda_T columnwise
        # to create [Y,Lambda_T]^T # TODO: is W transposed?
        self.__W = self.__cgraph.get_y().col_join(self.__cgraph.get_lambda().T)

    def __create_w_nullspace(self):
        # Finding the null space of [Y,Lambda_T]^T,
        # i.e. \omega_i i=1,...,\delta
        self.__W_nullspace = self.__W.nullspace()

    def __create_h_vector(self):

        print("Creating Equilibrium Manifold ...")
        start = time.process_time()
        # creates the H vector in equation (10) by finding the linearly
        # independent rows of (9)
        # symbolic for of psi vector
        psi_symbol_vec = sympy.zeros(self.__M, 1)
        for i in range(self.__M):
            psi_symbol_vec[i] = sympy.Symbol('psi' + str(i))

        # creating the right-hand side defined by
        # \sum_{i=1}^\delta \alpha_i \omega_i
        rhs = sympy.zeros(self.__M, 1)
        temp = sympy.zeros(self.__M, 1)
        for i in range(self.__delta):
            for j in range(self.__M):
                temp[j] = self.__deficiency_pars[i] * self.__W_nullspace[i][j]
            rhs += temp
        # creating H(c,\alpha,k) with rows that might
        # be linearly dependent
        temp_vec = self.__cgraph.get_a() * self.__cgraph.get_psi() - rhs
        temp_vec2 = self.__cgraph.get_a() * psi_symbol_vec - rhs

        # creating a matrix of coefficients for the variables
        # psi_1, ..., psi_M, alpha_1, ... alpha_delta
        variables = [sympy.Symbol('psi' + str(i)) for i in range(self.__M)] + self.__deficiency_pars
        a, _ = sympy.linear_eq_to_matrix(temp_vec2, variables)

        temp_mat_2 = sympy.zeros(self.__M, self.__M + self.__delta)

        # preallocating the H vector
        self.__H = sympy.zeros(self.__M - self.__ell, 1)

        # continue loop until acceptable concentration solutions are found
        flag = True
        indicies_to_skip = []
        while flag:

            reordered_indices, chosen_index = self.__create_fixed_free_pars_and_reordered_ind(temp_vec, indicies_to_skip)

            for i in range(self.__M):
                temp_mat_2[i, :] = a[reordered_indices[i], :]

            # gives vals which is the rows of temp_mat that are
            # independent from the others this allows us
            # to create H(c, \alpha, k)
            _, temp_vals = temp_mat_2.T.rref()
            vals = [reordered_indices[i] for i in temp_vals]

            # Filling the H vector with linearly independent rows
            for i in range(len(vals)):
                self.__H[i] = temp_vec[vals[i]]

            flag = self.__create_concentration_values()

            if flag:
                indicies_to_skip.append(chosen_index)

            nn = len(self.__concentration_pars)
            rr = self.__M - self.__ell - self.__delta
            if len(indicies_to_skip) == math.factorial(nn) / (math.factorial(rr) * math.factorial(nn-rr)):
                flag = False
                raise Exception("An analytic solution for the concentrations could not be found!")

        end = time.process_time()
        print("Elapsed time for creating Equilibrium Manifold: " + str(end - start))
        print("")

    def __create_fixed_free_pars_and_reordered_ind(self, temp_vec, indicies_to_skip):

        # determining the different combinations of concentrations
        # present in the independent variables once the deficiency
        # parameters are chosen to be independent
        leng = self.__M - self.__ell - self.__delta
        comb = list(itertools.combinations(self.__concentration_pars, leng))

        # building the different possible independent variable sets
        indp_vars = []
        for i in comb:
            indp_vars.append(list(i)+self.__deficiency_pars)

        # finding the number of linear equations produced by a
        # given independent variable set
        counts = []
        for jj in range(len(indp_vars)):

            if jj not in indicies_to_skip:

                num_lin_entries = [self.__is_linear(temp_vec[j], indp_vars[jj][0:leng]) for j in
                                   range(temp_vec.shape[0])].count(True)
                counts.append(num_lin_entries)

                # if all of the equations are linear stop,
                # prevents long run times
                if num_lin_entries == self.__M:
                    break

        # picking the independent variable set that has the most
        # amount of linear equations
        chosen_index = counts.index(max(counts))

        self.__fixed_pars = indp_vars[chosen_index]
        self.__free_pars = [i for i in self.__deficiency_pars + self.__concentration_pars if i not in self.__fixed_pars]

        # rearranging A s.t. the linear equations are in the top
        # rows, this is for convenience and easier solutions for
        # the independent variables
        out = [self.__is_linear(temp_vec[j], indp_vars[chosen_index]) for j in range(temp_vec.shape[0])]
        reordered_indices = [i for i, x in enumerate(out) if x] + [i for i, x in enumerate(out) if not x]

        return reordered_indices, chosen_index

    # routine that determines if a sympy expression is jointly
    # linear with respect to a given set of variables. This
    # test is conducted by seeing if the second order derivatives
    # are zero.
    def __is_linear(self, expr, variables):

        for x in variables:
            for y in variables:
                try: 
                    if not sympy.Eq(sympy.diff(expr, x, y), 0):
                        return False
                except TypeError:
                    return False
        return True

    def __create_g_matrix(self):

        # creating the matrix DCH which is the jacobian of
        # the vector H with respect to the concentration vector
        self.__DCH = self.__H.jacobian(sympy.Matrix(self.__concentration_pars))

        # creation of the matrix dah which is the jacobian of
        # the vector H with respect to the deficiency param vector
        dah = self.__H.jacobian(sympy.Matrix(self.__deficiency_pars))
        # creation of the matrix daw which is the jacobian of
        # the vector W with respect to the deficiency param vector
        # However, as given in page 11, D_\alpha W = 0, thus
        # just a zero matrix of size \lambda by \delta
        daw = sympy.zeros(self.__lambda, self.__delta)
        # creating the upper half of the matrix G i.e. [DCH dah]
        g_upper = self.__DCH.row_join(dah)
        # creating the lower half of the matrix G i.e. [DCW daw]
        # Note that D_c W = B^T
        g_lower = self.__cgraph.get_b().row_join(daw)
        # putting the upper and lower half of the matrix together
        # this forms the full G(c, \alpha, k) matrix
        self.__G = g_upper.col_join(g_lower)

    def __create_symbolic_objective_fun(self):
        # computing the simplified version of the objective
        # function defined as: det(G(c, \alpha, k))^2
        self.__symbolic_objective_fun = (self.__G.det(method='lu'))**2

    def __create_concentration_values(self):

        # Putting the concentrations in terms of the kinetic
        # constants and deficiency parameters using the H
        # vector of the equilibrium manifold
        temp_solution_tuple = sympy.solve(self.__H, self.__fixed_pars, dict=True)

        if not temp_solution_tuple:
            raise Exception("A solution for the fixed parameters could not be found!") 

        # constructing concentrations if only one solution is found
        if isinstance(temp_solution_tuple, dict):
            self.__concentration_vals = []
            for i in self.__concentration_pars:
                if i in temp_solution_tuple:
                    self.__concentration_vals.append(temp_solution_tuple[i])
                else:
                    self.__concentration_vals.append(i)
        # multiple solutions found
        else:
            solution_list = []
            for i in temp_solution_tuple:
                temp = []
                for j in self.__concentration_pars:
                    if j in i:
                        temp.append(i[j])
                    else:
                        temp.append(j)
                solution_list.append(temp)
            self.__concentration_vals = self.__pick_solution_set(solution_list)

        for i in self.__concentration_vals:

            deficiency_pars_found = [i.count(j) > 0 for j in self.__deficiency_pars]

            if True in deficiency_pars_found:
                flag = True
                break
            else:
                flag = False

        return flag

    # chose solution set that is most likely to produce
    # positive concentrations
    def __pick_solution_set(self, solution_list):
        positivity = []
        for i in solution_list:
            temp = []
            for j in i:
                temp.append(j.is_positive)
            positivity.append(temp)

        verdict = []
        for i in positivity:
            if False not in i:
                if all(i):
                    # positive concentrations achieved
                    verdict.append("P")
                elif True in i:
                    # positive concentrations and Nones
                    verdict.append("PU")
                else:
                    # All entries are None
                    verdict.append("U")
            else:
                # negative concentration given
                verdict.append("N")

        if "P" in verdict:
            indx = verdict.index("P")
            choice = solution_list[indx]
        elif "PU" in verdict:
            indx = verdict.index("PU")
            choice = solution_list[indx]
        elif "U" in verdict:
            indx = verdict.index("U")
            choice = solution_list[indx]
        else:
            print("Solution chosen produces all negative concentrations!")
            sys.exit()

        return choice 

    def __create_decision_vector_x(self):
        # if it is a proper/over-dimensioned network let
        # xvec = (k_1 , ... k_R, alpha_1, ... alpha_lambda)
        # else let
        # xvec = (k_1 , ... k_R, alpha_1, ... alpha_delta,
        # c_1, ..., c_(lambda - delta))
        self.__decision_vector_x = self.__reaction_pars + self.__free_pars
        
        self.__important_info = "The network is " + self.__classification[0] + "\n"

        self.__d_len = len(self.__decision_vector_x)

    def __create_concentration_bounds_species(self):
        self.__concentration_bounds_species = [i for i in self.__concentration_pars
                                               if i not in self.__decision_vector_x]

    def __create_concentration_lambda_fun(self):
        self.__concentration_funs = []
        for i in range(self.__N):
            self.__concentration_funs += [sympy.utilities.lambdify(self.__decision_vector_x,
                                                                   self.__concentration_vals[i])]

    def __create_objective_fun__lambda_fun(self):
        self.__objective_fun_params = self.__reaction_pars + self.__concentration_pars
        self.__lambda_objective_fun = sympy.utilities.lambdify(self.__objective_fun_params,
                                                               self.__symbolic_objective_fun)

    def __create_g_matrix_lambda_fun(self):
        self.__lambda_G_matrix = sympy.utilities.lambdify(self.__objective_fun_params, self.__G)
        
    def __create_dch_matrix_lambda_fun(self):
        self.__lambda_DCH_matrix = sympy.utilities.lambdify(self.__objective_fun_params, self.__DCH)

    # getters
    def get_w_matrix(self):
        """
        Returns SymPy matrix :math:`[Y, \Lambda^T]^T`, which we call the W matrix.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> sympy.pprint(approach.get_w_matrix())
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
            ⎢0  0  0  0  0  0  0  1  0⎥
            ⎢                         ⎥
            ⎢1  1  1  0  0  0  0  0  0⎥
            ⎢                         ⎥
            ⎢0  0  0  1  1  1  0  0  0⎥
            ⎢                         ⎥
            ⎣0  0  0  0  0  0  1  1  1⎦
        """
        return self.__W

    def get_w_nullspace(self):
        """
        Returns a list of SymPy column vectors representing :math:`Null([Y, \Lambda^T]^T)`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> sympy.pprint(approach.get_w_nullspace())
            ⎡⎡-1⎤  ⎡1 ⎤⎤
            ⎢⎢  ⎥  ⎢  ⎥⎥
            ⎢⎢0 ⎥  ⎢0 ⎥⎥
            ⎢⎢  ⎥  ⎢  ⎥⎥
            ⎢⎢1 ⎥  ⎢-1⎥⎥
            ⎢⎢  ⎥  ⎢  ⎥⎥
            ⎢⎢-1⎥  ⎢0 ⎥⎥
            ⎢⎢  ⎥  ⎢  ⎥⎥
            ⎢⎢0 ⎥, ⎢0 ⎥⎥
            ⎢⎢  ⎥  ⎢  ⎥⎥
            ⎢⎢1 ⎥  ⎢0 ⎥⎥
            ⎢⎢  ⎥  ⎢  ⎥⎥
            ⎢⎢0 ⎥  ⎢-1⎥⎥
            ⎢⎢  ⎥  ⎢  ⎥⎥
            ⎢⎢0 ⎥  ⎢0 ⎥⎥
            ⎢⎢  ⎥  ⎢  ⎥⎥
            ⎣⎣0 ⎦  ⎣1 ⎦⎦
        """
        return self.__W_nullspace

    def get_h_vector(self):
        """
        Returns a SymPy matrix representing the equilibrium manifold.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> sympy.pprint(approach.get_h_vector())
            ⎡a₁ - a₂ - re₁⋅s₁⋅s₂ + re1r⋅s₃⎤
            ⎢                             ⎥
            ⎢re₁⋅s₁⋅s₂ + s₃⋅(-re1r - re₂) ⎥
            ⎢                             ⎥
            ⎢  a₁ - re₃⋅s₆⋅s₇ + re3r⋅s₁₆  ⎥
            ⎢                             ⎥
            ⎢re₃⋅s₆⋅s₇ + s₁₆⋅(-re3r - re₄)⎥
            ⎢                             ⎥
            ⎢  a₂ - re₅⋅s₁⋅s₆ + re5r⋅s₁₅  ⎥
            ⎢                             ⎥
            ⎣re₅⋅s₁⋅s₆ + s₁₅⋅(-re5r - re₆)⎦
        """
        return self.__H

    def get_g_matrix(self):
        """
        Returns a SymPy matrix representing the G matrix of the defined optimization problem.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> sympy.pprint(approach.get_g_matrix())
            ⎡-re₁⋅s₂  -re₁⋅s₁     re1r         0        0          0            0       1  -1⎤
            ⎢                                                                                ⎥
            ⎢re₁⋅s₂   re₁⋅s₁   -re1r - re₂     0        0          0            0       0  0 ⎥
            ⎢                                                                                ⎥
            ⎢   0        0          0       -re₃⋅s₇  -re₃⋅s₆     re3r           0       1  0 ⎥
            ⎢                                                                                ⎥
            ⎢   0        0          0       re₃⋅s₇   re₃⋅s₆   -re3r - re₄       0       0  0 ⎥
            ⎢                                                                                ⎥
            ⎢-re₅⋅s₆     0          0       -re₅⋅s₁     0          0          re5r      0  1 ⎥
            ⎢                                                                                ⎥
            ⎢re₅⋅s₆      0          0       re₅⋅s₁      0          0       -re5r - re₆  0  0 ⎥
            ⎢                                                                                ⎥
            ⎢   0        0          0          0       1.0        1.0           0       0  0 ⎥
            ⎢                                                                                ⎥
            ⎢   0       1.0        1.0         0        0          0            0       0  0 ⎥
            ⎢                                                                                ⎥
            ⎣  1.0       0         1.0        1.0       0         1.0          2.0      0  0 ⎦
        """
        return self.__G

    def get_dch_matrix(self):
        """
        Returns a SymPy matrix representing the Jacobian of the equilibrium manifold with respect to the species.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> import sympy
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> sympy.pprint(approach.get_dch_matrix())
            ⎡-re₁⋅s₂  -re₁⋅s₁     re1r         0        0          0            0     ⎤
            ⎢                                                                         ⎥
            ⎢re₁⋅s₂   re₁⋅s₁   -re1r - re₂     0        0          0            0     ⎥
            ⎢                                                                         ⎥
            ⎢   0        0          0       -re₃⋅s₇  -re₃⋅s₆     re3r           0     ⎥
            ⎢                                                                         ⎥
            ⎢   0        0          0       re₃⋅s₇   re₃⋅s₆   -re3r - re₄       0     ⎥
            ⎢                                                                         ⎥
            ⎢-re₅⋅s₆     0          0       -re₅⋅s₁     0          0          re5r    ⎥
            ⎢                                                                         ⎥
            ⎣re₅⋅s₆      0          0       re₅⋅s₁      0          0       -re5r - re₆⎦
        """
        return self.__DCH

    def get_lambda_g_matrix(self):
        """
        Returns a lambda function representation of the G matrix. Here the arguments of the lambda function are given
        by the values provided by :func:`crnt4sbml.MassConservationApproach.get_objective_fun_params`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_lambda_g_matrix())
            <function _lambdifygenerated at 0x13248ac80>
        """
        return self.__lambda_G_matrix

    def get_lambda_dch_matrix(self):
        """
        Returns a lambda function representation of the Jacobian of the equilibrium manifold matrix. Here the
        arguments of the lambda function are given by the values provided by
        :func:`crnt4sbml.MassConservationApproach.get_objective_fun_params`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_lambda_dch_matrix())
            <function _lambdifygenerated at 0x131a06ea0>
        """
        return self.__lambda_DCH_matrix

    def get_symbolic_objective_fun(self):
        """
        Returns SymPy expression for the objective function of the optimization problem. This is the determinant of the
        G matrix squared.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_symbolic_objective_fun())
            1.0*re1**2*re2**2*re3**2*re4**2*re5**2*re6**2*s1**2*s6**2*s7**2*((1.0*s2/s7 - 1.0*s2*(-re3r - re4)/
            (re3*s6*s7))/re4 + (1.0 + 1.0*re1r/(re1*s1))/re2 + 1.0/(re1*s1))**2*(-((1.0*s6*(-1.0*s1/s6 + 1.0)/s7 +
            1.0 - (-re3r - re4)*(-1.0*s1/s6 + 1.0)/(re3*s7))/re4 + 1.0/re2)*(-1.0*re5r*s2/(re5*re6*s1*s6) -
            1.0*s2*(1 + re5*s6/(re1*s2))/(re5*s1*s6) - (1.0 + 1.0*re1r/(re1*s1))/re2)/((1.0*s2/s7 - 1.0*s2*
            (-re3r - re4)/(re3*s6*s7))/re4 + (1.0 + 1.0*re1r/(re1*s1))/re2 + 1.0/(re1*s1)) + (2.0 + 1.0*re5r/(re5*s6))/
            re6 + 1.0*(1 + re5*s6/(re1*s2))/(re5*s6) - 1.0/re2 - 1.0/(re1*s2))**2
        """
        return self.__symbolic_objective_fun

    def get_lambda_objective_fun(self):
        """
        Returns a lambda function representation of the objective function of the optimization problem. Here the
        arguments of the lambda function are given by the values provided by
        :func:`crnt4sbml.MassConservationApproach.get_objective_fun_params`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_lambda_objective_fun())
            <function _lambdifygenerated at 0x12f6f7ea0>
        """
        return self.__lambda_objective_fun

    def get_concentration_vals(self):
        """
        Returns a list of SymPy expressions representing the species in terms of those variables present in the decision
        vector. The order is that established in :func:`crnt4sbml.Cgraph.get_species`. Note that if only a single
        species is provided as an element in the list, this means the species is a free variable.

        See also
        ---------
        crnt4sbml.MassConservationApproach.get_concentration_solutions


        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_concentration_vals())
            [s15*(re5r + re6)/(re5*s6), s2, re1*s15*s2*(re5r + re6)/(re5*s6*(re1r + re2)), s6,
            -s15*(re5*re5r*s6*(re1r + re2)*(re3r + re4) - (re5r + re6)*(-re1*re1r*re3r*s2 - re1*re1r*re4*s2 +
            re1*re3r*s2*(re1r + re2) + re1*re4*s2*(re1r + re2) + re5*s6*(re1r + re2)*(re3r + re4)))/(re3*re4*re5*s6**2*
            (re1r + re2)), s15*(re1*re2*re5r*s2 + re1*re2*re6*s2 + re1r*re5*re6*s6 + re2*re5*re6*s6)/(re4*re5*s6*(re1r + re2)), s15]
        """
        return self.__concentration_vals

    def get_decision_vector(self):
        """
        Returns a list of SymPy variables that represent the decision vector of the optimization problem.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_decision_vector())
            [re1, re1r, re2, re3, re3r, re4, re5, re5r, re6, s2, s6, s15]
        """
        return self.__decision_vector_x

    def get_concentration_bounds_species(self):
        """
        Returns a list of SymPy variables that represents the order of species for the concentration bounds provided
        to :func:`crnt4sbml.MassConservationApproach.run_optimization`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_concentration_bounds_species())
            [s1, s3, s7, s16]
        """
        return self.__concentration_bounds_species

    def get_concentration_funs(self):
        """
        Returns a list of lambda functions representing each of the species. Here the species are those expressions
        provided by :func:`crnt4sbml.MassConservationApproach.get_concentration_vals` where the arguments of each
        lambda function is provided by :func:`crnt4sbml.MassConservationApproach.get_decision_vector`.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_concentration_funs())
            [<function _lambdifygenerated at 0x135f8b4d0>, <function _lambdifygenerated at 0x135f72050>,
            <function _lambdifygenerated at 0x135f728c0>, <function _lambdifygenerated at 0x135f725f0>,
            <function _lambdifygenerated at 0x135f5f830>, <function _lambdifygenerated at 0x135fa0170>,
            <function _lambdifygenerated at 0x135fa04d0>]
        """
        return self.__concentration_funs

    def get_objective_fun_params(self):
        """
        Returns a list of SymPy variables that represent those variables that may be contained in the G matrix, Jacobian
        of the equilibrium manifold with respect to the species, or objective function.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_objective_fun_params())
            [re1, re1r, re2, re3, re3r, re4, re5, re5r, re6, s1, s2, s3, s6, s7, s16, s15]
        """
        return self.__objective_fun_params

    def get_conservation_laws(self):
        """
        Returns a string representation of the conservation laws. Here the values on the left hand side of each equation
        are the constants of the conservation laws.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_conservation_laws())
            C1 = 1.0*s16 + 1.0*s7
            C2 = 1.0*s2 + 1.0*s3
            C3 = 1.0*s1 + 2.0*s15 + 1.0*s16 + 1.0*s3 + 1.0*s6
        """
        rhs = self.__cgraph.get_b() * sympy.Matrix([self.__concentration_pars]).T
        laws = ""
        for i in range(rhs.shape[0]):
            laws += 'C' + str(i+1) + ' = ' + str(rhs[i]) + '\n'
        
        return laws

    def get_concentration_solutions(self):
        """
        Returns a more readable string representation of the species defined in terms of the decision vector.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> print(approach.get_concentration_solutions())
            s1 = s15*(re5r + re6)/(re5*s6)
            s2 = s2
            s3 = re1*s15*s2*(re5r + re6)/(re5*s6*(re1r + re2))
            s6 = s6
            s7 = -s15*(re5*re5r*s6*(re1r + re2)*(re3r + re4) - (re5r + re6)*(-re1*re1r*re3r*s2 - re1*re1r*re4*s2 + re1*re3r*s2*(re1r + re2) + re1*re4*s2*(re1r + re2) + re5*s6*(re1r + re2)*(re3r + re4)))/(re3*re4*re5*s6**2*(re1r + re2))
            s16 = s15*(re1*re2*re5r*s2 + re1*re2*re6*s2 + re1r*re5*re6*s6 + re2*re5*re6*s6)/(re4*re5*s6*(re1r + re2))
            s15 = s15
        """
        sols = ""
        for i in range(self.__N):
            sols += self.__species[i] + ' = ' + str(self.__concentration_vals[i]) + '\n'            

        return sols

    def get_independent_odes(self):
        """
        Returns a SymPy Matrix where the rows represent the independent ODEs used in the numerical continuation routine. Here
        the entries of the list correspond to the time derivatives of the corresponding species provided by
        :func:`crnt4sbml.MassConservationApproach.get_independent_species`. Note that the independent ODEs created are
        based on the species chosen for the numerical continuation. Thus, the continuation routine needs to be ran
        first. If this function is called before the numerical continuation routine then None will be returned.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> approach = network.get_mass_conservation_approach()
        >>> multistable_param_ind = approach.run_greedy_continuity_analysis(species="species", parameters=params_for_global_min,
                                                                            auto_parameters={'PrincipalContinuationParameter': "PCP"})
        >>> odes = approach.get_independent_odes()
        """

        return self.__independent_odes

    def get_independent_species(self):
        """
        Returns a list of SymPy representations of the independent species used in the numerical continuation routine.
        Note that the independent species created are based on the species chosen for the numerical continuation. Thus,
        the continuation routine needs to be ran first. If this function is called before the numerical continuation
        routine then None will be returned.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> approach = network.get_mass_conservation_approach()
        >>> multistable_param_ind = approach.run_greedy_continuity_analysis(species="species", parameters=params_for_global_min,
                                                                            auto_parameters={'PrincipalContinuationParameter': "PCP"})
        >>> species = approach.get_independent_species()
        """

        return self.__independent_species

    def get_optimization_bounds(self):
        """
        Builds all of the necessary physiological bounds for the optimization routine.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Returns
        --------
        bounds: list of tuples
            List of tuples defining the upper and lower bounds for the decision vector variables based on physiological
            ranges.
        concentration_bounds: list of tuples
            List of tuples defining the upper and lower bounds for those concentrations not in the decision vector
            based on physiological ranges.

        Examples
        ---------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_mass_conservation_approach()
            Creating Equilibrium Manifold ...
            Elapsed time for creating Equilibrium Manifold: 2.060944

        >>> bounds, concentration_bounds = approach.get_optimization_bounds()
        >>> print(bounds)
            [(1e-08, 0.0001), (1e-05, 0.001), (0.001, 1.0), (1e-08, 0.0001), (1e-05, 0.001), (0.001, 1.0),
            (1e-08, 0.0001), (1e-05, 0.001), (0.001, 1.0), (0.5, 500000.0), (0.5, 500000.0), (0.5, 500000.0)]

        >>> print(concentration_bounds)
            [(0.5, 500000.0), (0.5, 500000.0), (0.5, 500000.0), (0.5, 500000.0)]
        """

        graph_edges = self.__cgraph.get_g_edges()
        dec_vec_var_def = []
        for i in self.get_decision_vector():

            if i in self.__concentration_pars:
                dec_vec_var_def.append("concentration")
            elif i in self.__reaction_pars:

                ind = self.__reaction_pars.index(i)
                reaction = graph_edges[ind]

                reaction_type = self.__cgraph.get_graph().edges[reaction]['type']
                dec_vec_var_def.append(reaction_type)

                if reaction_type is None:
                    output_statement = "The reaction type of reaction " + self.__cgraph.get_graph().edges[reaction]['label'] \
                                       + " could not be identified as it does not fit any biological criteria " +\
                                       "established. \n" + "You must enter bounds manually for this reaction! \n"
                    print(output_statement)

        concentration_bounds = [self.get_physiological_range("concentration")]*len(self.get_concentration_bounds_species())

        bounds = [self.get_physiological_range(i) for i in dec_vec_var_def]

        return bounds, concentration_bounds

    def run_optimization(self, bounds=None, iterations=10, sys_min_val=numpy.finfo(float).eps, seed=0, print_flag=False,
                         numpy_dtype=numpy.float64, concentration_bounds=None):
        """
        Function for running the optimization problem for the mass conservation approach.

        Parameters
        -----------
            bounds: list of tuples
                A list defining the lower and upper bounds for each variable in the decision vector. Here the reactions
                are allowed to be set to a single value.
            iterations: int
                The number of iterations to run the feasible point method.
            sys_min_val: float
                The value that should be considered zero for the optimization problem.
            seed: int
                Seed for the random number generator. None should be used if a random generation is desired.
            print_flag: bool
                Should be set to True if the user wants the objective function values found in the optimization problem
                and False otherwise.
            numpy_dtype:
                The numpy data type used within the optimization routine. All variables in the optimization routine will
                be converted to this data type.
            concentration_bounds: list of tuples
                A list defining the lower and upper bounds for those species' concentrations not in the decision vector.
                The user is not allowed to set the species' concentration to a single value. See also:
                :func:`crnt4sbml.MassConservationApproach.get_concentration_bounds_species`.
        Returns
        --------
        params_for_global_min: list of numpy arrays
            A list of numpy arrays that correspond to the decision vectors of the problem.
        obj_fun_val_for_params: list of floats
            A list of objective function values produced by the corresponding decision vectors in params_for_global_min.

        Examples
        ---------
        See :ref:`quickstart-deficiency-label` and :ref:`my-deficiency-label`.
        """
        self.__numpy_dtype = numpy_dtype
        temp_c = numpy.zeros(self.__N, dtype=self.__numpy_dtype)

        # testing to see if there are any equalities in bounds
        equality_bounds_indices = []
        for i in range(len(bounds)):
            if not isinstance(bounds[i], tuple):
                equality_bounds_indices.append(i)

        # recasting user provided input to numpy_dtype
        for i in range(len(bounds)):
            bounds[i] = self.__numpy_dtype(bounds[i])
            
        for i in range(len(concentration_bounds)):
            concentration_bounds[i] = self.__numpy_dtype(concentration_bounds[i])

        if len(concentration_bounds) != len(self.__concentration_bounds_species):
            print("Concentration bounds is the incorrect length!")
            sys.exit()

        full_concentration_bounds = []
        for i in range(self.__N):
            if self.__concentration_pars[i] in self.__decision_vector_x:
                indx = self.__decision_vector_x.index(self.__concentration_pars[i])
                full_concentration_bounds.append(bounds[indx])
            else:
                indx = self.__concentration_bounds_species.index(self.__concentration_pars[i])
                full_concentration_bounds.append(concentration_bounds[indx])

        sys_min_val = self.__numpy_dtype(sys_min_val)

        params_for_global_min, obj_fun_val_for_params, self.__important_info = BistabilityFinder.run_optimization(
            bounds, iterations, sys_min_val, temp_c, self.__penalty_objective_func, self.__feasible_point_check,
            self.__objective_function_to_optimize, self.__final_constraint_check, seed, equality_bounds_indices,
            print_flag, numpy_dtype, full_concentration_bounds)

        return params_for_global_min, obj_fun_val_for_params

    def run_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                print_lbls_flag=False, auto_parameters=None):
        """
        Function for running the numerical continuation and bistability analysis portions of the mass conservation
        approach.

        Parameters
        ------------
            species: string
                A string stating the species that is the y-axis of the bifurcation diagram.
            parameters: list of numpy arrays
                A list of numpy arrays corresponding to the decision vectors that produce a small objective function
                value.
            dir_path: string
                A string stating the path where the bifurcation diagrams should be saved.
            print_lbls_flag: bool
                If True the routine will print the special points found by AUTO 2000 and False will not print any
                special points.
            auto_parameters: dict
                Dictionary defining the parameters for the AUTO 2000 run. Please note that one should **not** set
                'SBML' or 'ScanDirection' in these parameters as these are automatically assigned. It is absolutely
                necessary to set PrincipalContinuationParameter in this dictionary. For more information on these
                parameters refer to :download:`AUTO parameters <../auto2000_input.pdf>`. 'NMX' will default to
                10000 and 'ITMX' to 100.
        Returns
        ---------
            multistable_param_ind: list of integers
                A list of those indices in 'parameters' that produce multistable plots.
            plot_specifications: list of lists
                A list whose elements correspond to the plot specifications of each element in multistable_param_ind.
                Each element is a list where the first element specifies the range used for the x-axis, the second
                element is the range for the y-axis, and the last element provides the x-y values and special point label
                for each special point in the plot.

        Example
        ---------
        See :ref:`quickstart-deficiency-label` and :ref:`my-deficiency-label`.
        """
        # setting default values for AUTO
        if 'NMX' not in auto_parameters.keys():
            auto_parameters['NMX'] = 10000

        if 'ITMX' not in auto_parameters.keys():
            auto_parameters['ITMX'] = 100

        # making the directory if it doesn't exist
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

        species_num = self.__species.index(species) + 1

        species_y = str(self.__concentration_pars[species_num-1]) 

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_continuity_analysis(species_num, parameters,
                                                                                                               self.__initialize_ant_string,
                                                                                                               self.__finalize_ant_string,
                                                                                                               species_y, dir_path,
                                                                                                               print_lbls_flag,
                                                                                                               auto_parameters)

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    def run_greedy_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                       print_lbls_flag=False, auto_parameters=None):
        """
        Function for running the greedy numerical continuation and bistability analysis portions of the mass conservation
        approach. This routine uses the initial value of the principal continuation parameter to construct AUTO
        parameters and then tests varying fixed step sizes for the continuation problem. Note that this routine may
        produce jagged or missing sections in the plots provided. To produce better plots one should use the information
        provided by this routine to run :func:`crnt4sbml.MassConservationApproach.run_continuity_analysis`.

        Parameters
        ------------
            species: string
                A string stating the species that is the y-axis of the bifurcation diagram.
            parameters: list of numpy arrays
                A list of numpy arrays corresponding to the decision vectors that produce a small objective function
                value.
            dir_path: string
                A string stating the path where the bifurcation diagrams should be saved.
            print_lbls_flag: bool
                If True the routine will print the special points found by AUTO 2000 and False will not print any
                special points.
            auto_parameters: dict
                Dictionary defining the parameters for the AUTO 2000 run. Please note that only the
                PrincipalContinuationParameter in this dictionary should be defined, no other AUTO parameters should
                be set. For more information on these parameters refer to :download:`AUTO parameters <../auto2000_input.pdf>`.
        Returns
        ---------
            multistable_param_ind: list of integers
                A list of those indices in 'parameters' that produce multistable plots.
            plot_specifications: list of lists
                A list whose elements correspond to the plot specifications of each element in multistable_param_ind.
                Each element is a list where the first element specifies the range used for the x-axis, the second
                element is the range for the y-axis, and the last element provides the x-y values and special point label
                for each special point in the plot.

        Example
        ---------
        See :ref:`my-deficiency-label`.
        """
        # setting default values for AUTO
        if 'NMX' not in auto_parameters.keys():
            auto_parameters['NMX'] = 10000

        if 'ITMX' not in auto_parameters.keys():
            auto_parameters['ITMX'] = 100

        # making the directory if it doesn't exist
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

        species_num = self.__species.index(species) + 1

        species_y = str(self.__concentration_pars[species_num-1])

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_greedy_continuity_analysis\
            (species_num, parameters, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters)

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    def generate_report(self):
        """
        Prints out helpful details constructed by :func:`crnt4sbml.MassConservationApproach.run_optimization` and
        :func:`crnt4sbml.MassConservationApproach.run_continuity_analysis`.

        Example
        --------
        See :ref:`quickstart-deficiency-label` and :ref:`my-deficiency-label`.
        """
        print(self.__important_info)

    def __initialize_ant_string(self, species_num, pcp_x):
        y = self.__cgraph.get_y()
        a = self.__cgraph.get_a()
        bt = self.__cgraph.get_b()
        psi = self.__cgraph.get_psi()

        # forming ya matrix
        ya = y*a

        # finding how many rows are indep in ya
        _, vals = ya.T.rref()
        num_indp_eqns = len(vals)
        num_dep_eqns = ya.shape[0] - num_indp_eqns

        # getting dimensions of bt
        bt_rows = bt.shape[0]
        bt_cols = bt.shape[1]

        bt_nonzero_ind = []
        for i in range(bt_rows):
            bt_nonzero_ind.append([j for j in range(bt_cols) if bt[i, j] != 0 and j != species_num-1])

        chosen_indp_indices, chosen_dep_indices = self.__get_indp_dep_species_indices(bt_nonzero_ind, num_dep_eqns,
                                                                                      num_indp_eqns, ya)

        replacements, ind_spec_conc_temp, indp_odes_temp = self.__construct_important_variables(chosen_indp_indices,
                                                                                                chosen_dep_indices, ya,
                                                                                                psi, bt)

        ode_str = self.__create_ode_str(replacements, ind_spec_conc_temp, indp_odes_temp, species_num)

        return ode_str, pcp_x

    def __get_indp_dep_species_indices(self, bt_nonzero_ind, num_dep_eqns, num_indp_eqns, ya):
        # getting all combinations of the list indices
        possible_dep_species = list(itertools.product(*bt_nonzero_ind))

        removed_entries = []
        # remove tuples that have duplicate entries
        for i in range(len(possible_dep_species)):
            if len(set(possible_dep_species[i])) != num_dep_eqns:
                removed_entries.append(i)
        for index in sorted(removed_entries, reverse=True):
            del possible_dep_species[index]

        # get corresponding possible dependent species
        possible_indp_species = [] 
        species_ind = [i for i in range(len(self.__concentration_pars))]

        for i in possible_dep_species:
            possible_indp_species.append([j for j in species_ind if j not in i])

        # using YA to pick one of the possible indices
        chosen_indp_indices = []
        chosen_dep_indices = []
        for i in range(len(possible_indp_species)):
            _, vals = ya[possible_indp_species[i], :].T.rref()

            if len(vals) == num_indp_eqns: 
                chosen_indp_indices = possible_indp_species[i]
                chosen_dep_indices = possible_dep_species[i]
                break

        return chosen_indp_indices, chosen_dep_indices

    def __construct_important_variables(self, chosen_indp_indices, chosen_dep_indices, ya, psi, bt):
        # getting independent concentrations
        ind_spec_conc_temp = [self.__concentration_pars[i] for i in chosen_indp_indices]

        # getting dependent concentrations
        dep_spec_conc = [self.__concentration_pars[i] for i in chosen_dep_indices]

        # constructing the independent ODEs
        indp_odes_temp = ya[chosen_indp_indices, :] * psi

        # creating conservation laws string
        self.__cons_laws_sympy = bt*sympy.Matrix([self.__concentration_pars]).T

        # Lambda function of conservation laws
        self.__cons_laws_lamb = [sympy.utilities.lambdify(self.__concentration_pars, self.__cons_laws_sympy[i])
                                 for i in range(len(self.__cons_laws_sympy))]

        cons_laws_sympy_eq = [sympy.Eq(sympy.Symbol('C' + str(i+1), real=True), self.__cons_laws_sympy[i])
                              for i in range(len(self.__cons_laws_sympy))]

        dep_conc_in_laws = self.__dependent_species_concentrations(self.__cons_laws_sympy, dep_spec_conc)
        
        replacements = self.__find_dep_concentration_replacements(dep_conc_in_laws, self.__cons_laws_sympy,
                                                                  dep_spec_conc, cons_laws_sympy_eq)

        return replacements, ind_spec_conc_temp, indp_odes_temp

    def __create_ode_str(self, replacements, ind_spec_conc_temp, indp_odes_temp, species_num):
        # rearrange ind_spec_conc and indp_odes to make species of
        # interest be the first ODE
        indx_species_num = ind_spec_conc_temp.index(self.__concentration_pars[species_num-1])

        self.__ind_spec_conc = [ind_spec_conc_temp[indx_species_num]]
        for i in ind_spec_conc_temp:
            if i != self.__concentration_pars[species_num-1]:
                self.__ind_spec_conc.append(i)

        indp_odes = sympy.zeros(indp_odes_temp.shape[0], indp_odes_temp.shape[1])
        indp_odes[0] = indp_odes_temp[indx_species_num]
        count = 1
        for i in range(indp_odes_temp.shape[0]):
            if i != indx_species_num:
                indp_odes[count] = indp_odes_temp[i]
                count += 1

        # bulding ODE string in Antimony format
        ode_str = self.__building_ode_str(replacements, self.__ind_spec_conc, indp_odes)

        return ode_str
        
    def __finalize_ant_string(self, x, ode_str):
        concentration_vals = [self.__concentration_funs[j](*tuple(x)) for j in range(self.__N)]

        kinetic_vals = [x[i] for i in range(self.__R)] 

        antstr = self.__initialize_variables_in_antimony_string(self.__cons_laws_sympy, ode_str,
                                                                self.__cons_laws_lamb, concentration_vals, kinetic_vals,
                                                                self.__reaction_pars)

        return antstr

    def __final_constraint_check(self, x_initial, penalty_bounds, sys_min_val, equality_bounds_indices,
                                 concentration_bounds):
        non_equality_bounds_indices = [i for i in range(len(penalty_bounds)) if i not in equality_bounds_indices]

        x = numpy.zeros(len(penalty_bounds), dtype=self.__numpy_dtype)
        for j in equality_bounds_indices:              
            x[j] = penalty_bounds[j]              
        count = 0                                      
        for j in non_equality_bounds_indices:          
            x[j] = x_initial[count]                
            count += 1 

        # concentration > 0 check
        con = numpy.asarray([self.__concentration_funs[j](*tuple(x)) for j in range(self.__N)],
                            dtype=self.__numpy_dtype)
        con_temp = []
        for i in range(self.__N):
            con_temp.append(con[i] >= concentration_bounds[i][0] and con[i] <= concentration_bounds[i][1])
        concs_chk = numpy.all(con_temp)

        # boundary check
        test = []
        for j in non_equality_bounds_indices: 
            test.append(x[j] >= penalty_bounds[j][0] and x[j] <= penalty_bounds[j][1])
        boundry_chk = numpy.all(test)

        # rank(G) = N + delta - 1 check
        # xx = numpy.concatenate((x[0:self.__R],con),axis=None)

        # must convert xx to numpy.float64 because higher
        # is not supported in linalg
        # xx = numpy.float64(xx)

        # rank_G = numpy.linalg.matrix_rank(self.__lambda_G_matrix(*tuple(xx)))
        # rank_G_chk = rank_G == (self.__N + self.__delta - 1)

        # rank(DCH) = min(N,M-ell) check
        # rank_DCH = numpy.linalg.matrix_rank(self.__lambda_DCH_matrix(*tuple(xx)))

        # rank_DCH_chk = rank_DCH == min(self.__N,self.__M - self.__ell)

        if concs_chk and boundry_chk:  # and rank_G_chk and rank_DCH_chk:
            
            return [True, x]
        else:
            
            return [False, []]

    def __concentration_violation_fun(self, g, len_g):
        temp = numpy.zeros(len_g, dtype=self.__numpy_dtype)
        for i in range(len_g):
            temp[i] = numpy.maximum(self.__numpy_dtype(0.0), -g[i])**2
        return temp

    def __x_violation_fun(self, x, b, len_x):
        temp = numpy.zeros(len_x, dtype=self.__numpy_dtype)
        for i in range(len_x):
            temp[i] = numpy.maximum(self.__numpy_dtype(0.0), self.__numpy_dtype(b)-x[i])**2
        return temp

    def __penalty_objective_func(self, x_initial, temp_c, penalty_bounds, equality_bounds_indices, x,
                                 non_equality_bounds_indices, concentration_bounds):

        for j in equality_bounds_indices:
            x[j] = penalty_bounds[j]
        count = 0 
        for j in non_equality_bounds_indices:
            x[j] = x_initial[count]
            count += 1 

        # evaluating the concentrations first
        for i in range(self.__N):
            temp_val = self.__concentration_funs[i](*tuple(x))
            if numpy.iscomplex(temp_val): 
                temp_c = numpy.array([numpy.Inf for i in range(self.__N)], dtype=self.__numpy_dtype)
                break
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", numpy.ComplexWarning)
                    temp_c[i] = temp_val 

        if numpy.all(numpy.isfinite(temp_c)):

            # obtaining the sum of the violation functions squared
            sumval = self.__numpy_dtype(0.0)
            for j in range(self.__N):
                sumval += numpy.maximum(self.__numpy_dtype(0.0), concentration_bounds[j][0]-temp_c[j])**2
                sumval += numpy.maximum(self.__numpy_dtype(0.0), temp_c[j] - concentration_bounds[j][1])**2

            sum0 = self.__numpy_dtype(0.0)
            for j in non_equality_bounds_indices: 
                sum0 += numpy.maximum(self.__numpy_dtype(0.0), penalty_bounds[j][0]-x[j])**2
                sum0 += numpy.maximum(self.__numpy_dtype(0.0), x[j] - penalty_bounds[j][1])**2
            sumval += sum0

            # obtaining the violation function values for
            # k's and concentrations in x
            # xx = numpy.concatenate((x[0:self.__R],x[self.__R + self.__alpha_end_ind:self.__d_len]),axis=None)
            # temp = self.__x_violation_fun(xx,self.__numpy_dtype(0.0),self.__R + (self.__d_len - (self.__R +
            # self.__alpha_end_ind)))
            # sumval += numpy.sum(temp)
            return sumval 

        else:
            return numpy.PINF

    def __feasible_point_check(self, x, result_fun, sys_min_val, equality_bounds_indices, non_equality_bounds_indices,
                               penalty_bounds, concentration_bounds):
        result_x = numpy.zeros(len(penalty_bounds), dtype=self.__numpy_dtype)

        for j in equality_bounds_indices:
            result_x[j] = penalty_bounds[j]
        count = 0
        for j in non_equality_bounds_indices:
            result_x[j] = x[count]
            count += 1

        # double checking the concentrations
        con = numpy.asarray([self.__concentration_funs[i](*tuple(result_x)) for i in range(self.__N)],
                            dtype=self.__numpy_dtype)
        con_temp = []
        for i in range(self.__N):
            con_temp.append(con[i] >= concentration_bounds[i][0] and con[i] <= concentration_bounds[i][1])
        concs_chk = numpy.all(con_temp)

        finite_chk = numpy.isfinite(con)
        if concs_chk and numpy.all(finite_chk):

            # putting the feasible points in x_candidates
            if abs(result_fun) <= sys_min_val and numpy.all(con > self.__numpy_dtype(0)):
                return True
            else:
                return False
        else:
            return False 
                
    def __objective_function_to_optimize(self, x_initial, temp_c, penalty_bounds, sys_min_val, equality_bounds_indices,
                                         x, non_equality_bounds_indices, concentration_bounds):
        for j in equality_bounds_indices:
            x[j] = penalty_bounds[j]
        count = 0
        for j in non_equality_bounds_indices:
            x[j] = x_initial[count]
            count += 1

        test = []
        for j in non_equality_bounds_indices:
            test.append(x[j] >= penalty_bounds[j][0] and x[j] <= penalty_bounds[j][1])
        boundry_chk = numpy.all(test)

        if boundry_chk:
            # calculating the concentration values
            for i in range(self.__N):
                temp_val = self.__concentration_funs[i](*tuple(x))
                if numpy.iscomplex(temp_val):
                    temp_c = numpy.array([numpy.Inf for i in range(self.__N)], dtype=self.__numpy_dtype)
                    break
                else:
                    with warnings.catch_warnings(): 
                        warnings.simplefilter("ignore", numpy.ComplexWarning)
                        temp_c[i] = temp_val 

            finite_chk = numpy.isfinite(temp_c)
            con_temp = []
            for i in range(self.__N):
                con_temp.append(temp_c[i] >= concentration_bounds[i][0] and temp_c[i] <= concentration_bounds[i][1])
            concs_chk = numpy.all(con_temp)
        # making sure our concentrations are finite
            if concs_chk and numpy.all(finite_chk): 
                temp = numpy.zeros(self.__N, dtype=self.__numpy_dtype)
                for i in range(self.__N):
                    temp[i] = numpy.maximum(self.__numpy_dtype(0.0), -temp_c[i])
                sumval = numpy.sum(temp)

                xx = numpy.concatenate((x[0:self.__R], temp_c), axis=None)
                return self.__lambda_objective_fun(*tuple(xx)) + sumval
            else:
                return numpy.PINF

        else:
            return numpy.PINF

    def __dependent_species_concentrations(self, cons_laws_sympy, dep_spec_conc):
        # finding those dep_spec_conc that occur in each conservation law
        dep_conc_in_laws = []
        for i in range(len(cons_laws_sympy)):
            temp = []
            for j in range(len(dep_spec_conc)):
                if cons_laws_sympy[i].count(dep_spec_conc[j]) > 0:
                    temp.append(dep_spec_conc[j])
            dep_conc_in_laws.append(temp)
        return dep_conc_in_laws

    def __is_list_empty(self, inlist):
        if isinstance(inlist, list):  # Is a list
            return all(map(self.__is_list_empty, inlist))
        return False  # Not a list

    def __find_dep_concentration_replacements(self, dep_conc_in_laws, cons_laws_sympy, dep_spec_conc,
                                              cons_laws_sympy_eq):
        replacements = []
        flag = True
        while flag:
            for i in range(len(cons_laws_sympy_eq)):
                if len(dep_conc_in_laws[i]) == 1:
                    temp = sympy.solve(cons_laws_sympy_eq[i], dep_conc_in_laws[i])
                    cons_laws_sympy = [cons_laws_sympy[j].subs(dep_conc_in_laws[i][0], temp[0])
                                       for j in range(len(cons_laws_sympy))]

                    cons_laws_sympy_eq = [sympy.Eq(sympy.Symbol('C' + str(i+1), real=True), cons_laws_sympy[i])
                                          for i in range(len(cons_laws_sympy))]

                    replacements.append([dep_conc_in_laws[i][0], '(' + str(temp[0]) + ')'])

            dep_conc_in_laws = self.__dependent_species_concentrations(cons_laws_sympy, dep_spec_conc)

            if self.__is_list_empty(dep_conc_in_laws):
                flag = False

        return replacements    

    def __building_ode_str(self, replacements, ind_spec_conc, indp_odes):
        indp_odes_str = []
        # making the replacements in the indep. ODEs
        for i in range(len(indp_odes)):
            for j in range(len(replacements)):
                indp_odes[i] = indp_odes[i].subs(replacements[j][0], replacements[j][1])
            indp_odes_str.append(str(indp_odes[i]))

        self.__independent_odes = indp_odes
        self.__independent_species = ind_spec_conc

        # replacing all powers with ^ instead of **
        for i in range(len(indp_odes_str)):
            indp_odes_str[i] = indp_odes_str[i].replace('**', '^')

        # building the string of ODEs in Antimony syntax
        ode_str = ''
        for i in range(len(ind_spec_conc)):
            ode_str += 'J' + str(i) + ': -> ' + str(ind_spec_conc[i]) + '; ' + indp_odes_str[i] + ';'

        return ode_str

    def __building_ant_str(self, ode_str, kinetic_con, lhs_cons_laws, var_vals):
        vars_to_initialize = kinetic_con+lhs_cons_laws + [str(self.__concentration_pars[i]) for i in range(self.__N)]

        ant_str = ode_str

        for i in range(len(vars_to_initialize)):
            ant_str += str(vars_to_initialize[i]) + ' = ' + str(var_vals[i]) + ';'

        return ant_str

    def __initialize_variables_in_antimony_string(self, cons_laws_sympy, ode_str, cons_laws_lamb, concentration_vals,
                                                  kinetic_vals, kinetic_con):

        # string representation of variables on lhs of mass cons laws
        lhs_cons_laws = ['C'+str(i+1) for i in range(len(cons_laws_sympy))]

        conservation_law_vals = [cons_laws_lamb[i](*tuple(concentration_vals)) for i in range(len(cons_laws_lamb))]

        var_vals = kinetic_vals + conservation_law_vals + concentration_vals

        # The full Antimony string of system of ODEs
        ant_str = self.__building_ant_str(ode_str, kinetic_con, lhs_cons_laws, var_vals)

        return ant_str
