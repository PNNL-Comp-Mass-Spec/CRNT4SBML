import sympy
import sys
import itertools
import numpy
import scipy.optimize
import scipy.linalg
import warnings
import time
import math
import os
from .bistability_finder import BistabilityFinder


class GeneralApproach:
    """
    Class for constructing a more general approach to bistability detection for systems with mass action kinetics.
    """
    def __init__(self, cgraph, signal, response, fix_reactions):
        """
        Initialization of GeneralApproach class using a network contructed by crnt4sbml.CRNT().

        See also
        ---------
        fill in this
        """

        self.__cgraph = cgraph

        self.__signal = signal
        self.__response = response
        self.__fix_reactions = fix_reactions

        if not self.__cgraph.get_dim_equilibrium_manifold() > 0:
            print("# of species - rank(S) is not greater than zero.")
            print("This implies that there are no mass conservation laws.")
            print("The general approach cannot be ran!")
            sys.exit()

        # todo: check that there is more than one independent ODE

        self.__full_system = self.__cgraph.get_ode_system()

        self.__sympy_species = [sympy.Symbol(i, positive=True) for i in self.__cgraph.get_species()]
        self.__sympy_reactions = [sympy.Symbol(i, positive=True) for i in self.__cgraph.get_reactions()]

        self.__N = len(self.__sympy_species)
        self.__R = len(self.__sympy_reactions)
        self.__M = len(self.__cgraph.get_complexes())
        self.__BT_mat = None

        self.__A_mat = self.__cgraph.get_a()
        self.__Y_mat = self.__cgraph.get_y()
        self.__psi_vec = self.__cgraph.get_psi()
        self.__S_mat = self.__cgraph.get_s()

        try:
            self.__BT_mat = self.__cgraph.get_b()
        except:
            self.__create_BT_matrix()

        self.__cons_laws_sympy = None
        self.__indp_species = None
        self.__indp_system = None
        self.__indp_system_subs = None
        self.__replacements = None
        self.__det_jac_indp_subs = None
        self.__D2_g_u_w = None
        self.__decision_vector = None
        self.__lambda_indp_odes = None
        self.__lambda_det_jac = None
        self.__lambda_jec = None
        self.__lambda_D2_g_u_w = None
        self.__lambda_jac_eps = None
        self.__important_info = ""
        self.__lambda_variables = None
        self.__indices_in_dec_vec = None
        self.__comm = None
        self.__my_rank = None
        self.__num_cores = None

        self.__create_conservation_laws()
        self.__construct_independent_system()
        self.__construct_independent_system_with_substitution()
        self.__create_variables_for_optimization()

    def __create_BT_matrix(self):
        # (bt = NullSpace[Transpose[Y.A]]) // MatrixForm
        the_null_space = (self.__Y_mat * self.__A_mat).T.nullspace()

        # handling an edge case where the null space calculation isn't fully simplified by SymPy
        if len(the_null_space) > len(self.__sympy_species) - self.__S_mat.rank():
            the_null_space = self.__simplify_nullspace(the_null_space)

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
        bt = bt[0:sizes, :]

        # putting in a check to make sure that the rows
        # of bt are linearly independent
        _, vals = bt.T.rref()
        if len(vals) != sizes:
            print(" ")
            print("Error: Not all rows of B^T are linearly independent!")
            print(" ")
            sys.exit()

        self.__BT_mat = bt

    def __simplify_nullspace(self, the_null_space):

        # Flag that says a species or reaction was found in one of the column vectors of the nullspace
        denom_flag = False

        # find those columns that have a common denominator
        denom_list = []
        for ii in the_null_space:

            # getting common denominator
            denom_in_column = [sympy.denom(i) for i in ii]
            denom_list.append(denom_in_column[0])

        # collecting nullspace vectors with common denominators
        all_denom = list(set(denom_list))
        temp_null = [sympy.zeros(self.__N, 1)] * (len(self.__sympy_species) - self.__S_mat.rank())
        for ii in range(len(all_denom)):

            indices = [i for i, x in enumerate(denom_list) if x == all_denom[ii]]
            for i in indices:
                temp_null[ii] += the_null_space[i]
                temp_null[ii] = sympy.simplify(temp_null[ii])

        # checking to see if a reaction or species is in the reduced nullspace
        for i in temp_null:
            atoms_of_col_vec = [ii.atoms() for ii in i[:, 0]]
            for ii in atoms_of_col_vec:
                if any([iii in self.__sympy_reactions + self.__sympy_species for iii in ii]):
                    denom_flag = True

        if denom_flag:
            raise Exception("Nullspace calculation from S contains SymPy variables and this could not be resolved.")

        return temp_null

    def __create_non_negative_b_matrix(self, the_null_space, sizes):

        a_null = numpy.zeros((self.__N, sizes))
        for i in range(self.__N):
            temp_vec = []
            for j in range(sizes):
                temp_vec.append(the_null_space[j][i])
            a_null[i, :] = temp_vec

        a_eq = numpy.array([numpy.sum(a_null, axis=0)])
        # must multiply by negative one because in optimization we have the inequality <= 0.0
        a_ub = -1.0 * a_null
        b_ub = numpy.zeros(self.__N)
        b_eq = numpy.array([1.0])

        # defining the number of solutions to simulate
        num_sols = ((a_ub.shape[0] + 1) * (a_ub.shape[1] + 1)) * 10

        # a matrix to hold the different solutions of the method
        sol = numpy.zeros((num_sols, self.__N))

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
            unique_vecs[i, :] = unique_vecs[i, :] / numpy.float64(minval)

        unique_vecs = numpy.unique(unique_vecs.round(decimals=10), axis=0)

        return unique_vecs

    def __create_conservation_laws(self):

        self.__cons_laws_sympy = self.__BT_mat*sympy.Matrix(self.__sympy_species)

    def __construct_independent_system(self):

        # forming ya matrix
        #ya = self.__Y_mat * self.__A_mat   # old approach
        ya, __ = sympy.linear_eq_to_matrix(self.__Y_mat * self.__A_mat*self.__psi_vec, self.__sympy_reactions) # new

        # finding how many rows are indep in ya
        __, vals = ya.T.rref()

        num_indp_eqns = len(vals)
        num_dep_eqns = self.__N - num_indp_eqns

        # getting dimensions of bt
        bt_rows = self.__BT_mat.shape[0]
        bt_cols = self.__BT_mat.shape[1]

        bt_nonzero_ind = []
        species_num = self.__sympy_species.index(sympy.Symbol(self.__response, positive=True))
        for i in range(bt_rows):
            bt_nonzero_ind.append([j for j in range(bt_cols) if self.__BT_mat[i, j] != 0 and j != species_num])

        chosen_indp_indices, chosen_dep_indices = self.__get_indp_dep_species_indices(bt_nonzero_ind, num_dep_eqns,
                                                                                      num_indp_eqns, ya)

        self.__replacements, self.__indp_species, self.__indp_system = self.__construct_important_variables(chosen_indp_indices,
                                                                                                            chosen_dep_indices, ya)

    def __construct_independent_system_with_substitution(self):

        self.__indp_system_subs = sympy.zeros(len(self.__indp_system), 1)
        self.__indp_system_subs = self.__indp_system[:, :]

        # making the replacements in the indep. ODEs
        for i in range(len(self.__indp_system_subs)):
            for j in range(len(self.__replacements)):
                self.__indp_system_subs[i] = self.__indp_system_subs[i].subs(self.__replacements[j][0], self.__replacements[j][1])

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
        species_ind = [i for i in range(self.__N)]

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

    def __construct_important_variables(self, chosen_indp_indices, chosen_dep_indices, ya):
        # getting independent concentrations
        ind_spec_conc_temp = [self.__sympy_species[i] for i in chosen_indp_indices]

        # getting dependent concentrations
        dep_spec_conc = [self.__sympy_species[i] for i in chosen_dep_indices]

        # constructing the independent ODEs
        indp_odes_temp = ya[chosen_indp_indices, :] * sympy.Matrix(self.__sympy_reactions) # new approach
        #indp_odes_temp = ya[chosen_indp_indices, :] * self.__psi_vec  # old approach

        # Lambda function of conservation laws
        self.__cons_laws_lamb = [sympy.utilities.lambdify(self.__sympy_species, self.__cons_laws_sympy[i])
                                 for i in range(len(self.__cons_laws_sympy))]

        cons_laws_sympy_eq = [sympy.Eq(sympy.Symbol('C' + str(i + 1), real=True), self.__cons_laws_sympy[i])
                              for i in range(len(self.__cons_laws_sympy))]

        dep_conc_in_laws = self.__dependent_species_concentrations(self.__cons_laws_sympy, dep_spec_conc)

        replacements = self.__find_dep_concentration_replacements(dep_conc_in_laws, self.__cons_laws_sympy,
                                                                  dep_spec_conc, cons_laws_sympy_eq)

        return replacements, ind_spec_conc_temp, indp_odes_temp

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

                    replacements.append([dep_conc_in_laws[i][0], temp[0]])

            dep_conc_in_laws = self.__dependent_species_concentrations(cons_laws_sympy, dep_spec_conc)

            if self.__is_list_empty(dep_conc_in_laws):
                flag = False

        return replacements

    def __create_variables_for_optimization(self):

        self.__jac_indp_subs = self.__indp_system_subs.jacobian(sympy.Matrix(self.__indp_species))

        u_variables = sympy.Matrix(self.__indp_species)
        w_vector = sympy.Matrix([sympy.Symbol('w' + str(i)) for i in range(len(u_variables))])

        self.__D2_g_u_w = self.second_derivative_operator(self.__indp_system_subs, u_variables, w_vector)

        self.__det_jac_indp_subs = self.__jac_indp_subs.det(method='lu')

        self.__lambda_variables = self.__sympy_reactions + self.__indp_species + \
                                 [sympy.Symbol("C" + str(i+1), real=True) for i in range(len(self.__cons_laws_sympy))]

        self.__jac_eps = self.__indp_system_subs.jacobian([sympy.Symbol(self.__signal, real=True)])

        self.__lambda_jac_eps = sympy.utilities.lambdify(self.__lambda_variables, self.__jac_eps)

        self.__lambda_indp_odes = [sympy.utilities.lambdify(self.__lambda_variables, i) for i in self.__indp_system_subs]

        self.__lambda_det_jac = sympy.utilities.lambdify(self.__lambda_variables, self.__det_jac_indp_subs)

        self.__lambda_jec = sympy.utilities.lambdify(self.__lambda_variables, self.__jac_indp_subs)

        self.__lambda_D2_g_u_w = sympy.utilities.lambdify(self.__lambda_variables + [sympy.Symbol('w' + str(i))
                                                                                    for i in range(len(u_variables))],
                                                          self.__D2_g_u_w)

        if self.__fix_reactions:
            self.__enforce_steady_state()
            # self.__check()
        else:
            self.__fixed_reactions = []
            self.__vars_for_lam_fixed_reactions = []
            self.__lambda_fixed_reactions = []
            self.__soln_to_fixed_reactions2 = []

        self.__decision_vector = [i for i in self.__sympy_reactions if i not in self.__fixed_reactions] + self.__sympy_species

    def __check(self):

        temp_jac_eps = self.__jac_eps[:, :]
        k_vec = sympy.Matrix([[sympy.Symbol('k' + str(i))]  for i in range(self.__jac_eps.shape[0])])
        temp_D2_g_u_w = self.__D2_g_u_w[:, :]*k_vec

        jac_eps_subs = []
        for i in range(self.__jac_eps.shape[0]):
            temp_jac_eps[i] = [temp_jac_eps[i].subs(self.get_fixed_reactions()[j],
                                                    self.get_solutions_to_fixed_reactions()[j]) for j in
                               range(len(self.get_fixed_reactions()))][0]

            temp_jac_eps[i] = sympy.simplify(sympy.expand(temp_jac_eps[i]))

            if sympy.S.Zero == temp_jac_eps[i]:
                jac_eps_subs.append(True)
            else:
                jac_eps_subs.append(False)

        D2_g_u_w_subs = []
        for i in range(self.__D2_g_u_w.shape[0]):
            temp_D2_g_u_w[i] = [temp_D2_g_u_w[i].subs(self.get_fixed_reactions()[j],
                                                      self.get_solutions_to_fixed_reactions()[j])
                                for j in range(len(self.get_fixed_reactions()))][0]

            temp_D2_g_u_w[i] = sympy.simplify(sympy.expand(temp_D2_g_u_w[i]))

            if sympy.S.Zero == sympy.simplify(temp_D2_g_u_w[i]):
                D2_g_u_w_subs.append(True)
            else:
                D2_g_u_w_subs.append(False)

        a, b = sympy.linear_eq_to_matrix(temp_D2_g_u_w, self.__sympy_reactions)

        D2_nullspace = a.nullspace()

        for i in range(len(D2_nullspace)):
            D2_nullspace[i] = sympy.simplify(D2_nullspace[i])

        a, b = sympy.linear_eq_to_matrix(temp_jac_eps, self.__sympy_reactions)

        jac_eps_nullspace = a.nullspace()

        for i in range(len(jac_eps_nullspace)):
            jac_eps_nullspace[i] = sympy.simplify(jac_eps_nullspace[i])


    @staticmethod
    def second_derivative_operator(g_functions, u_variables, w_vector):

        Dg_u = g_functions.jacobian(u_variables)

        # operating on w_vector
        Dg_u_w = Dg_u*w_vector

        D2_g_u_w = Dg_u_w.jacobian(u_variables)

        return D2_g_u_w

    @staticmethod
    def unique_everseen(iterable, key=None):
        "List unique elements, preserving order. Remember all elements ever seen."
        # unique_everseen('AAAABBBCCDAABBB') --> A B C D
        # unique_everseen('ABBCcAD', str.lower) --> A B C D
        seen = set()
        seen_add = seen.add
        if key is None:
            for element in itertools.filterfalse(seen.__contains__, iterable):
                seen_add(element)
                yield list(element)
        else:
            for element in iterable:
                k = key(element)
                if k not in seen:
                    seen_add(k)
                    yield list(element)

    def __construct_possible_fixed_reactions(self):

        # obtaining the reactions in each ODE
        reactions_lists = []
        for i in self.__indp_system_subs:
            atoms = list(i.atoms(sympy.Symbol))
            reactions = [j for j in self.__sympy_reactions if j in atoms]
            reactions_lists.append(reactions)

        # finding all combinations of the reactions for each ODE minus those with repeating elements
        temp = [p for p in itertools.product(*reactions_lists) if len(set(p)) == len(p)]

        # finds unique elements of the list where the elements are seen as unordered sets
        possible_fixed_reactions = list(self.unique_everseen(temp, key=frozenset))

        return possible_fixed_reactions

    def __enforce_steady_state(self):

        self.__fixed_reactions = []
        a, vals = self.__S_mat.rref()

        self.__fixed_reactions = [self.__sympy_reactions[i] for i in vals]
        a, b = sympy.linear_eq_to_matrix(self.__indp_system, self.__fixed_reactions)
        temp_solution_tuple = list(sympy.linsolve((a, b), self.__fixed_reactions))[0]

        no_zero_values_flag = True

        # looking for other fixed reaction sets if no solution was found
        if sympy.S.Zero in temp_solution_tuple:
            print("")
            print("At least one fixed reaction found to be zero. \nAttempting to find nonzero values for all fixed reactions.")

            no_zero_values_flag = False
            reactions_found_to_be_zero = []

            zero_indices = [i for i in range(len(temp_solution_tuple)) if temp_solution_tuple[i] == sympy.S.Zero]
            [reactions_found_to_be_zero.append(self.__fixed_reactions[i]) for i in zero_indices]

            comb = self.__construct_possible_fixed_reactions()
            zeros = sympy.zeros(len(self.__indp_system), 1)
            for i in comb:
                i.sort(key=lambda j: self.__sympy_reactions.index(j))
                a, b = sympy.linear_eq_to_matrix(self.__indp_system, i)
                rank_a = a.rank()
                aug = a.row_join(zeros)
                rank_aug = aug.rank()
                self.__fixed_reactions = i
                if rank_a == rank_aug and rank_a == a.shape[1]:
                    self.__fixed_reactions = i
                    temp_solution_tuple2 = list(sympy.linsolve((a, b), self.__fixed_reactions))[0]
                    if sympy.S.Zero not in temp_solution_tuple2:
                        no_zero_values_flag = True
                        break
                    else:
                        zero_indices = [ii for ii in range(len(temp_solution_tuple2)) if
                                        temp_solution_tuple2[ii] == sympy.S.Zero]
                        [reactions_found_to_be_zero.append(i[ii]) for ii in zero_indices if i[ii] not in reactions_found_to_be_zero]
            print("Done searching for nonzero values for all fixed reactions.")
            print("")

        if no_zero_values_flag:
            self.__soln_to_fixed_reactions2 = list(temp_solution_tuple)
        else:
            print("There was no solution which resulted in all fixed reaction values being nonzero.")
            print(f"Consider removing one or all of the following reactions {reactions_found_to_be_zero}. \n")
            self.__soln_to_fixed_reactions2 = list(temp_solution_tuple2)

        for i in range(len(self.__soln_to_fixed_reactions2)):
            for j in range(len(self.__replacements)):
                self.__soln_to_fixed_reactions2[i] = self.__soln_to_fixed_reactions2[i].subs(self.__replacements[j][0], self.__replacements[j][1])

        self.__vars_for_lam_fixed_reactions = [i for i in self.__lambda_variables if i not in self.__fixed_reactions]
        self.__lambda_fixed_reactions = [sympy.utilities.lambdify(self.__vars_for_lam_fixed_reactions, i) for i in
                                         self.__soln_to_fixed_reactions2]

    def __feasible_point_obj_func(self, x, boundz, full_set_of_values, reaction_ind, cons_laws_lamb, indp_spec_ind, inputs):

        # filling in values for non-fixed reactions
        inputs[0:self.__R - len(reaction_ind)] = x[0:self.__R - len(reaction_ind)]

        # filling in values for independent species
        count = self.__R - len(reaction_ind)
        for i in indp_spec_ind:
            inputs[count] = x[i]
            count += 1

        # filling in values for conservation laws
        for i in cons_laws_lamb:
            inputs[count] = i(*tuple(x[self.__R - len(reaction_ind):]))
            count += 1

        # filling in all species values
        full_set_of_values[self.__R:] = x[self.__R - len(reaction_ind):]

        # filling in non-fixed reaction values
        count = 0
        for i in range(self.__R):
            if i not in reaction_ind:
                full_set_of_values[i] = x[count]
                count += 1

        # filling in fixed reaction values
        full_set_of_values[reaction_ind] = numpy.array([i(*tuple(inputs)) for i in self.__lambda_fixed_reactions])

        sumval = numpy.float64(0.0)
        for i in range(len(full_set_of_values)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(boundz[i][0]) - full_set_of_values[i])
            sumval += numpy.maximum(numpy.float64(0.0), full_set_of_values[i] - numpy.float64(boundz[i][1]))

        #cons_law = 0
        #sumval += numpy.maximum(numpy.float64(0.0), cons_laws_lamb[cons_law](*tuple(x[self.__R - len(reaction_ind):]) + numpy.float64(1e-2)))      ############### carefullll
        #sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(-100.0) - cons_laws_lamb[cons_law](*tuple(x[self.__R - len(reaction_ind):])))

        return sumval

    def __obj_func(self, x, jac_det, boundz, full_set_of_values, reaction_ind, cons_laws_lamb, indp_spec_ind, inputs, input_values):

        # filling in values for non-fixed reactions
        inputs[0:self.__R - len(reaction_ind)] = x[0:self.__R - len(reaction_ind)]

        # filling in values for independent species
        count = self.__R - len(reaction_ind)
        for i in indp_spec_ind:
            inputs[count] = x[i]
            count += 1

        # filling in values for conservation laws
        for i in cons_laws_lamb:
            inputs[count] = i(*tuple(x[self.__R - len(reaction_ind):]))
            count += 1

        # filling in all species values
        full_set_of_values[self.__R:] = x[self.__R - len(reaction_ind):]

        # filling in non-fixed reaction values
        count = 0
        for i in range(self.__R):
            if i not in reaction_ind:
                full_set_of_values[i] = x[count]
                count += 1

        # filling in fixed reaction values
        full_set_of_values[reaction_ind] = numpy.array([i(*tuple(inputs)) for i in self.__lambda_fixed_reactions])

        # filling in reactions
        input_values[0: self.__R] = full_set_of_values[0:self.__R]

        # filling in independent species
        input_values[self.__R: self.__R + len(indp_spec_ind)] = x[indp_spec_ind]

        # filling in conservation laws
        input_values[self.__R + len(indp_spec_ind):] = inputs[self.__R - len(reaction_ind)+ len(indp_spec_ind):]

        fun = numpy.float64(0.0)

        determinant = (jac_det(*tuple(input_values)))**2

        if numpy.isfinite(determinant):
            fun += determinant
        else:
            fun += numpy.PINF

        sumval = numpy.float64(0.0)
        for i in range(len(full_set_of_values)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(boundz[i][0]) - full_set_of_values[i]) #** 2
            sumval += numpy.maximum(numpy.float64(0.0), full_set_of_values[i] - numpy.float64(boundz[i][1])) #** 2

        cons_law = 0
        #sumval += numpy.maximum(numpy.float64(0.0), cons_laws_lamb[cons_law](*tuple(x[self.__R - len(reaction_ind):]) + numpy.float64(1e-2)))  ############### carefullll
        #sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(-100.0) - cons_laws_lamb[cons_law](*tuple(x[self.__R - len(reaction_ind):])))

        fun += sumval

        if not self.__fix_reactions:

            for i in self.__lambda_indp_odes:
                fun += (i(*tuple(input_values)))**2

        return fun

    def __feasible_point_method(self, bounds, num_constraint_method_iters, seed, print_flag, confidence_level_flag):

        all_vars_with_bounds = self.__sympy_reactions + self.__sympy_species
        decision_vector_bounds_ind = [all_vars_with_bounds.index(i) for i in self.__decision_vector]

        decision_vector_bounds = [bounds[i] for i in decision_vector_bounds_ind]

        samples = numpy.random.rand(num_constraint_method_iters, len(self.__decision_vector))

        ranges = numpy.asarray(decision_vector_bounds, dtype=numpy.float64)
        samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

        # import math
        # ranges = numpy.asarray([(math.log10(i[0]), math.log10(i[1])) for i in decision_vector_bounds],
        #                        dtype=numpy.float64)
        # samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]
        # samples = numpy.power(10, samples)

        full_set_of_values = numpy.zeros(len(bounds), dtype=numpy.float64)

        fixed_reaction_ind_all = [all_vars_with_bounds.index(i) for i in self.__fixed_reactions]
        indp_spec_ind_dec = [self.__decision_vector.index(i) for i in self.__indp_species]

        if self.__fix_reactions:
            inputs = numpy.zeros(len(self.__vars_for_lam_fixed_reactions), dtype=numpy.float64)

        feasible_point_sets = []
        if self.__fix_reactions:
            for i in range(num_constraint_method_iters):

                with numpy.errstate(divide='ignore', invalid='ignore'):

                    result = scipy.optimize.dual_annealing(self.__feasible_point_obj_func, bounds=decision_vector_bounds,
                                                           args=(bounds, full_set_of_values, fixed_reaction_ind_all,
                                                                 self.__cons_laws_lamb, indp_spec_ind_dec, inputs),
                                                           x0=samples[i],
                                                           seed=seed,
                                                           local_search_options={'method': "Nelder-Mead"}, maxiter=100)

                    if print_flag:
                        print("Objective function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)

                    if abs(result.fun) > numpy.float64(1e-100):
                        result1 = scipy.optimize.minimize(self.__feasible_point_obj_func, result.x,
                                                          args=(bounds, full_set_of_values, fixed_reaction_ind_all,
                                                                self.__cons_laws_lamb, indp_spec_ind_dec, inputs),
                                                          method='Nelder-Mead', tol=1e-16)

                        if print_flag:
                            print("Objective function value: " + str(result1.fun))
                            print("Decision vector used: ")
                            print(result1.x)

                        if abs(result1.fun) <= numpy.finfo(float).eps or confidence_level_flag:
                            feasible_point_sets.append(result1.x)

                    else:
                        feasible_point_sets.append(result.x)

                    if print_flag:
                        print("")
        else:
            feasible_point_sets = [i for i in samples]

        return feasible_point_sets, fixed_reaction_ind_all, indp_spec_ind_dec, decision_vector_bounds

    def __mpi_feasible_point_method(self, bounds, seed, print_flag, confidence_level_flag,
                                    samples, all_vars_with_bounds, decision_vector_bounds):

        full_set_of_values = numpy.zeros(len(bounds), dtype=numpy.float64)

        fixed_reaction_ind_all = [all_vars_with_bounds.index(i) for i in self.__fixed_reactions]
        indp_spec_ind_dec = [self.__decision_vector.index(i) for i in self.__indp_species]
        inputs = numpy.zeros(len(self.__vars_for_lam_fixed_reactions), dtype=numpy.float64)

        feasible_point_sets = []
        for i in samples:

            with numpy.errstate(divide='ignore', invalid='ignore'):

                result = scipy.optimize.dual_annealing(self.__feasible_point_obj_func, bounds=decision_vector_bounds,
                                                       args=(bounds, full_set_of_values, fixed_reaction_ind_all,
                                                             self.__cons_laws_lamb, indp_spec_ind_dec, inputs),
                                                       x0=i,
                                                       seed=seed,
                                                       local_search_options={'method': "Nelder-Mead"}, maxiter=100)

                if print_flag:
                    print("Objective function value: " + str(result.fun))
                    print("Decision vector used: ")
                    print(result.x)

                if abs(result.fun) > numpy.float64(1e-100):
                    result1 = scipy.optimize.minimize(self.__feasible_point_obj_func, result.x,
                                                      args=(bounds, full_set_of_values, fixed_reaction_ind_all,
                                                            self.__cons_laws_lamb, indp_spec_ind_dec, inputs),
                                                      method='Nelder-Mead', tol=1e-16)

                    if print_flag:
                        print("Objective function value: " + str(result1.fun))
                        print("Decision vector used: ")
                        print(result1.x)

                    if abs(result1.fun) <= numpy.finfo(float).eps or confidence_level_flag:
                        feasible_point_sets.append(result1.x)

                else:
                    feasible_point_sets.append(result.x)

                if print_flag:
                    print("")

        return feasible_point_sets, fixed_reaction_ind_all, indp_spec_ind_dec

    def run_mpi_optimization(self, bounds=None, iterations=10, seed=0, print_flag=False, dual_annealing_iters=1000,
                             confidence_level_flag=False, change_in_rel_error=1e-2):

        """
        Function for running the optimization problem for the mass conservation approach.

        Parameters
        -----------
            bounds: list of tuples
                A list defining the lower and upper bounds for each variable in the decision vector.
            iterations: int
                The number of iterations to run the feasible point method.
            seed: int
                Seed for the random number generator. None should be used if a random generation is desired.
            print_flag: bool
                Should be set to True if the user wants the objective function values found in the optimization problem
                and False otherwise.
            dual_annealing_iters: integer
                The number of iterations that should be ran for dual annealing routine in optimization.
            confidence_level_flag: bool
                If True a confidence level for the objective function will be given.
            change_in_rel_error: float
                The maximum relative error that should be allowed to consider :math:`f_k` in the neighborhood
                of :math:`\widetilde{f}`.
        Returns
        --------
        params_for_global_min: list of numpy arrays
            A list of numpy arrays that correspond to the decision vectors of the problem.
        obj_fun_val_for_params: list of floats
            A list of objective function values produced by the corresponding decision vectors in params_for_global_min.

        Examples
        ---------
        Fill in
        """

        # initializing MPI proccess
        from mpi4py import MPI
        self.__comm = MPI.COMM_WORLD
        self.__my_rank = self.__comm.Get_rank()
        self.__num_cores = self.__comm.Get_size()

        self.__comm.Barrier()

        start_time = MPI.Wtime()

        if self.__my_rank == 0:
            print("Starting optimization ...")

        all_vars_with_bounds = self.__sympy_reactions + self.__sympy_species
        decision_vector_bounds_ind = [all_vars_with_bounds.index(i) for i in self.__decision_vector]

        decision_vector_bounds = [bounds[i] for i in decision_vector_bounds_ind]

        # creating initial decision vectors for feasible point method
        if self.__my_rank == 0:
            numpy.random.seed(seed)
            samples = numpy.random.rand(iterations, len(self.__decision_vector))
            ranges = numpy.asarray(decision_vector_bounds, dtype=numpy.float64)
            samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

            # import math
            # ranges = numpy.asarray([(math.log10(i[0]), math.log10(i[1])) for i in decision_vector_bounds],
            #                        dtype=numpy.float64)
            # samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]
            # samples = numpy.power(10, samples)
        else:
            samples = None

        sample_portion = self.__distribute_points(samples)

        if self.__my_rank == 0:
            print("Running feasible point method for " + str(iterations) + " iterations ...")

        feasible_point_sets, fixed_reaction_ind_all, indp_spec_ind_dec = self.__mpi_feasible_point_method(bounds, seed,
                                                                                                          print_flag,
                                                                                                          confidence_level_flag,
                                                                                                          sample_portion,
                                                                                                          all_vars_with_bounds,
                                                                                                          decision_vector_bounds)

        # checking number of elements of feasible_point_sets for each core to see if we need to redistribute them
        redistribute_flag = len(feasible_point_sets) == len(sample_portion)
        val = self.__comm.allreduce(redistribute_flag, op=MPI.LAND)
        if not val:
            array_of_feasibles = self.__gather_numpy_array_of_values(feasible_point_sets)
            feasible_point_sets = self.__distribute_points(array_of_feasibles)

        self.__comm.Barrier()

        if self.__my_rank == 0:
            print("Feasible point method has finished.")

        if self.__my_rank != 0:
            self.__important_info += f"The number of feasible points used in determinant optimization by core {self.__my_rank}: " \
                                    + str(len(feasible_point_sets))
        else:
            self.__important_info += f"The number of feasible points used in determinant optimization by core {self.__my_rank}: " \
                                     + str(len(feasible_point_sets)) + "\n"

        det_point_sets, det_point_sets_fun, obtained_minimums, full_set_of_values, inputs, input_values, smallest_value = \
            self.__main_optimization_routine(decision_vector_bounds, bounds, feasible_point_sets, fixed_reaction_ind_all
                                             , indp_spec_ind_dec, confidence_level_flag, seed, dual_annealing_iters,
                                             print_flag)

        self.__comm.Barrier()

        if self.__my_rank == 0:
            end_time = MPI.Wtime()
            elapsed = end_time - start_time
            print(f"Elapsed time for optimization in seconds: {elapsed}")

        params = self.__final_check(det_point_sets, bounds, full_set_of_values, fixed_reaction_ind_all,
                                    self.__cons_laws_lamb, indp_spec_ind_dec, inputs, input_values)

        if confidence_level_flag:

            full_obtained_minimums = self.__gather_single_value(obtained_minimums, iterations)

            if self.__my_rank == 0:
                self.__confidence_level(full_obtained_minimums, change_in_rel_error)

        else:
            smallest_values = self.__comm.gather(smallest_value, root=0)
            if self.__my_rank == 0:
                min_value = min(smallest_values)
                self.__important_info += "Smallest value achieved by objective function: " + str(min_value)

        list_params = self.__gather_list_of_values(params)
        list_det_point_sets_fun = self.__gather_list_of_values(det_point_sets_fun)

        self.__comm.Barrier()

        return list_params, list_det_point_sets_fun, self.__my_rank

    def __gather_single_value(self, value, number_of_values):

        temp_full_value = self.__comm.gather(value, root=0)

        if self.__my_rank == 0:
            full_value = numpy.zeros(number_of_values, dtype=numpy.float64)

            # putting all the obtained minimum values into a single array
            count = 0
            for i in temp_full_value:
                for j in i:
                    full_value[count] = j
                    count += 1
        else:
            full_value = None

        return full_value

    def __gather_list_of_values(self, values):

        full_values = self.__comm.gather(values, root=0)

        if self.__my_rank == 0:
            list_of_values = []
            for i in range(len(full_values)):
                list_of_values += full_values[i]
        else:
            list_of_values = []

        return list_of_values

    def __gather_numpy_array_of_values(self, values):

        full_values = self.__comm.gather(values, root=0)

        if self.__my_rank == 0:
            list_of_values = []
            for i in range(len(full_values)):
                list_of_values += full_values[i]
            array_of_values = numpy.zeros((len(list_of_values), list_of_values[0].shape[0]), dtype=numpy.float64)
            for i in range(len(list_of_values)):
                array_of_values[i, :] = list_of_values[i]
        else:
            array_of_values = None

        return array_of_values

    def __distribute_points(self, samples):

        if self.__my_rank == 0:

            # number of tasks per core
            tasks = len(samples) // self.__num_cores  # // calculates the floor

            # remainder
            r = len(samples) - self.__num_cores * tasks

            # array that holds how many tasks each core has
            tasks_core = numpy.zeros(self.__num_cores, dtype=numpy.int64)
            tasks_core.fill(tasks)

            # distributing in the remainder
            ii = 0
            while r > 0:
                tasks_core[ii] += 1
                r -= 1
                ii += 1

            sample_portion = samples[0:tasks_core[0], :]

            if self.__num_cores > 1:
                for i in range(1, self.__num_cores):
                    start = sum(tasks_core[0:i])
                    end = start + tasks_core[i]
                    self.__comm.send(samples[start:end, :], dest=i, tag=i * 11)

        else:
            if self.__num_cores > 1:
                sample_portion = self.__comm.recv(source=0, tag=self.__my_rank * 11)

        return sample_portion

    def run_optimization(self, bounds=None, iterations=10, seed=0, print_flag=False, dual_annealing_iters=1000,
                         confidence_level_flag=False, change_in_rel_error=1e-2):

        """
        Function for running the optimization problem for the mass conservation approach.

        Parameters
        -----------
            bounds: list of tuples
                A list defining the lower and upper bounds for each variable in the decision vector.
            iterations: int
                The number of iterations to run the feasible point method.
            seed: int
                Seed for the random number generator. None should be used if a random generation is desired.
            print_flag: bool
                Should be set to True if the user wants the objective function values found in the optimization problem
                and False otherwise.
            dual_annealing_iters: integer
                The number of iterations that should be ran for dual annealing routine in optimization.
            confidence_level_flag: bool
                If True a confidence level for the objective function will be given.
            change_in_rel_error: float
                The maximum relative error that should be allowed to consider :math:`f_k` in the neighborhood
                of :math:`\widetilde{f}`.
        Returns
        --------
        params_for_global_min: list of numpy arrays
            A list of numpy arrays that correspond to the decision vectors of the problem.
        obj_fun_val_for_params: list of floats
            A list of objective function values produced by the corresponding decision vectors in params_for_global_min.

        Examples
        ---------
        Fill in
        """

        print("Starting optimization ...")
        start_t = time.time()
        numpy.random.seed(seed)

        print("Running feasible point method for " + str(iterations) + " iterations ...")
        feasible_point_sets, fixed_reaction_ind_all, indp_spec_ind_dec, decision_vector_bounds = self.__feasible_point_method(bounds, iterations, seed, print_flag, confidence_level_flag)
        print("Feasible point method has finished.")

        self.__important_info += "\nThe number of feasible points used in determinant optimization: " + str(len(feasible_point_sets)) + "\n"

        bounds_2 = decision_vector_bounds
        full_set_of_values = numpy.zeros(len(bounds), dtype=numpy.float64)

        if self.__fix_reactions:
            inputs = numpy.zeros(len(self.__vars_for_lam_fixed_reactions), dtype=numpy.float64)
        else:
            inputs = numpy.zeros(len(self.__lambda_variables), dtype=numpy.float64)

        input_values = numpy.zeros(len(self.__lambda_variables), dtype=numpy.float64)

        det_point_sets = []
        det_point_sets_fun = []
        smallest_value = numpy.float(1e8)
        if confidence_level_flag:
            obtained_minimums = numpy.zeros(len(feasible_point_sets), dtype=numpy.float64)

        if len(feasible_point_sets) != 0:

            print("Starting determinant optimization ...")
            for i in range(len(feasible_point_sets)):

                with numpy.errstate(divide='ignore', invalid='ignore'):

                    result = scipy.optimize.dual_annealing(self.__obj_func, bounds=bounds_2, args=(self.__lambda_det_jac, bounds,
                                                                                    full_set_of_values, fixed_reaction_ind_all,
                                                                                    self.__cons_laws_lamb, indp_spec_ind_dec,
                                                                                    inputs, input_values), x0=feasible_point_sets[i], seed=seed,
                                                           local_search_options={'method': "Nelder-Mead"}, maxiter=dual_annealing_iters)

                    if print_flag:
                        print("Global function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)

                    if abs(result.fun) > numpy.float64(1e-100):

                        result1 = scipy.optimize.minimize(self.__obj_func, result.x, args=(self.__lambda_det_jac, bounds,
                                                                                           full_set_of_values, fixed_reaction_ind_all,
                                                                                           self.__cons_laws_lamb, indp_spec_ind_dec,
                                                                                           inputs, input_values),
                                                         method='Nelder-Mead', tol=1e-16)

                        if print_flag:
                            print("Local function value: " + str(result1.fun))
                            print("Decision vector used: ")
                            print(result1.x)

                        if smallest_value > result1.fun:
                            smallest_value = result1.fun

                        if confidence_level_flag:
                            obtained_minimums[i] = result1.fun

                        if abs(result1.fun) <= numpy.finfo(float).eps:
                            det_point_sets.append(result1.x)
                            det_point_sets_fun.append(result1.fun)

                    else:
                        if smallest_value > result.fun:
                            smallest_value = result.fun
                        det_point_sets.append(result.x)
                        det_point_sets_fun.append(result.fun)
                        if confidence_level_flag:
                            obtained_minimums[i] = result.fun

                    if print_flag:
                        print("")
            print("Determinant optimization has finished.")

            end_t = time.time()
            elapsed = end_t - start_t
            print("Elapsed time for optimization in seconds: " + str(elapsed))

            params = self.__final_check(det_point_sets, bounds, full_set_of_values, fixed_reaction_ind_all,
                                        self.__cons_laws_lamb, indp_spec_ind_dec, inputs, input_values)

            if confidence_level_flag:
                self.__confidence_level(obtained_minimums, change_in_rel_error)

            else:
                self.__important_info += "Smallest value achieved by objective function: " + str(smallest_value) + "\n"

            return params, det_point_sets_fun

        else:
            raise Exception("Optimization needs to be run with more iterations or different bounds.")


    def __main_optimization_routine(self, decision_vector_bounds, bounds, feasible_point_sets, fixed_reaction_ind_all,
                                    indp_spec_ind_dec, confidence_level_flag, seed, dual_annealing_iters, print_flag):

        self.__comm.Barrier()

        if self.__my_rank == 0:
            print("Starting determinant optimization ...")

        bounds_2 = decision_vector_bounds

        inputs = numpy.zeros(len(self.__vars_for_lam_fixed_reactions), dtype=numpy.float64)
        full_set_of_values = numpy.zeros(len(bounds), dtype=numpy.float64)

        input_values = numpy.zeros(len(self.__lambda_variables), dtype=numpy.float64)

        det_point_sets = []
        det_point_sets_fun = []
        smallest_value = numpy.float(1e8)

        obtained_minimums = numpy.zeros(len(feasible_point_sets), dtype=numpy.float64)

        if len(feasible_point_sets) != 0:

            for i in range(len(feasible_point_sets)):

                with numpy.errstate(divide='ignore', invalid='ignore'):

                    result = scipy.optimize.dual_annealing(self.__obj_func, bounds=bounds_2,
                                                           args=(self.__lambda_det_jac, bounds,
                                                                 full_set_of_values, fixed_reaction_ind_all,
                                                                 self.__cons_laws_lamb, indp_spec_ind_dec,
                                                                 inputs, input_values), x0=feasible_point_sets[i],
                                                           seed=seed,
                                                           local_search_options={'method': "Nelder-Mead"},
                                                           maxiter=dual_annealing_iters)

                    if print_flag:
                        print("Global function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)

                    if abs(result.fun) > numpy.float64(1e-100):

                        result1 = scipy.optimize.minimize(self.__obj_func, result.x,
                                                          args=(self.__lambda_det_jac, bounds,
                                                                full_set_of_values, fixed_reaction_ind_all,
                                                                self.__cons_laws_lamb, indp_spec_ind_dec,
                                                                inputs, input_values),
                                                          method='Nelder-Mead', tol=1e-16)

                        if print_flag:
                            print("Local function value: " + str(result1.fun))
                            print("Decision vector used: ")
                            print(result1.x)

                        if smallest_value > result1.fun:
                            smallest_value = result1.fun

                        if confidence_level_flag:
                            obtained_minimums[i] = result1.fun

                        if abs(result1.fun) <= numpy.finfo(float).eps:
                            det_point_sets.append(result1.x)
                            det_point_sets_fun.append(result1.fun)

                    else:
                        if smallest_value > result.fun:
                            smallest_value = result.fun
                        det_point_sets.append(result.x)
                        det_point_sets_fun.append(result.fun)
                        if confidence_level_flag:
                            obtained_minimums[i] = result.fun

                    if print_flag:
                        print("")

            self.__comm.Barrier()
            if self.__my_rank == 0:
                print("Determinant optimization has finished.")


            return det_point_sets, det_point_sets_fun, obtained_minimums, full_set_of_values, inputs, input_values, smallest_value

        else:
            raise Exception("Optimization needs to be run with more iterations or different bounds.")


    def __confidence_level(self, obtained_minimums, change_in_rel_error):

        a = 1
        b = 5

        unique_elements, counts_elements = numpy.unique(obtained_minimums, return_counts=True)
        min_val_index = numpy.nanargmin(unique_elements)

        f_til = unique_elements[min_val_index]

        numpy_float64_smalles_positive_value = numpy.nextafter(numpy.float64(0), numpy.float64(1))

        if f_til > numpy_float64_smalles_positive_value:

            r = numpy.count_nonzero(
                numpy.abs(f_til - obtained_minimums) / f_til < numpy.float64(change_in_rel_error))

            n_til = obtained_minimums.shape[0]
            a_bar = a + b - 1
            b_bar = b - r - 1

            prob = 1 - (math.factorial(n_til + a_bar) * math.factorial(2 * n_til + b_bar)) / (
                    math.factorial(2 * n_til + a_bar) * math.factorial(n_til + b_bar))

        else:

            prob = 1.0

        self.__important_info += f"It was found that {unique_elements[min_val_index]} is the minimum objective function value with a confidence level of {prob} .\n"

    def __final_check(self, det_point_sets, boundz, full_set_of_values, reaction_ind,
                      cons_laws_lamb, indp_spec_ind, inputs, input_values):

        params = []
        for x in det_point_sets:

            # filling in values for non-fixed reactions
            inputs[0:self.__R - len(reaction_ind)] = x[0:self.__R - len(reaction_ind)]

            # filling in values for independent species
            count = self.__R - len(reaction_ind)
            for i in indp_spec_ind:
                inputs[count] = x[i]
                count += 1

            # filling in values for conservation laws
            for i in cons_laws_lamb:
                inputs[count] = i(*tuple(x[self.__R - len(reaction_ind):]))
                count += 1

            # filling in all species values
            full_set_of_values[self.__R:] = x[self.__R - len(reaction_ind):]

            # filling in non-fixed reaction values
            count = 0
            for i in range(self.__R):
                if i not in reaction_ind:
                    full_set_of_values[i] = x[count]
                    count += 1

            # filling in fixed reaction values
            full_set_of_values[reaction_ind] = numpy.array([i(*tuple(inputs)) for i in self.__lambda_fixed_reactions])

            # filing in reactions
            input_values[0: self.__R] = full_set_of_values[0:self.__R]

            # filling in independent species
            input_values[self.__R: self.__R + len(indp_spec_ind)] = x[indp_spec_ind]

            # filling in conservation laws
            input_values[self.__R + len(indp_spec_ind):] = inputs[self.__R - len(reaction_ind) + len(indp_spec_ind):]

            check = self.__saddle_node_bifurcation(input_values)
            # print("saddle node check:")
            # print(check)
            # print("")
            if True: #check:
                params.append(input_values.flatten())

        #sys.exit()

        return params


    def __saddle_node_bifurcation(self, result_x):

        its_zero = numpy.float(1e-14) # necessary to compare against zero

        jac_w_vals = self.__lambda_jec(*tuple(result_x)).astype(numpy.float64)
        jac_nullspace = scipy.linalg.null_space(jac_w_vals).astype(numpy.float64)

        # check for one dimensional nullspace
        if jac_nullspace.shape[1] == 1:

            k_vec = jac_nullspace.astype(numpy.float64)

            # check for nonzero nullspace
            if any([abs(i) > its_zero for i in k_vec]):

                eigs = numpy.linalg.eigvals(jac_w_vals).real.astype(numpy.float64)
                ind_zero_sum = sum([int(1) for i in eigs if abs(i) <= its_zero])

                # check for only one zero eigenvalue
                if ind_zero_sum == 1:

                    # check for steady state
                    if all([i(*tuple(result_x)) <= its_zero for i in self.__lambda_indp_odes]):

                        jac_eps_val = self.__lambda_jac_eps(*tuple(result_x)).astype(numpy.float64)

                        k_vec_vals = [numpy.float64(i) for i in k_vec]
                        D2_g_u_w_val = self.__lambda_D2_g_u_w(*tuple(list(result_x) + k_vec_vals)).astype(numpy.float64)

                        Duu_kk = numpy.dot(D2_g_u_w_val, k_vec).astype(numpy.float64)

                        # check for nonzero f_eps and f_uu
                        if any([abs(i) > its_zero for i in jac_eps_val]) \
                                and any([abs(i) > its_zero for i in Duu_kk]):

                            aug_eps = numpy.column_stack((jac_w_vals, jac_eps_val))
                            aug_Duu = numpy.column_stack((jac_w_vals, Duu_kk))

                            jac_w_vals_rank = numpy.linalg.matrix_rank(jac_w_vals)
                            aug_eps_rank = numpy.linalg.matrix_rank(aug_eps)
                            aug_Duu_rank = numpy.linalg.matrix_rank(aug_Duu)

                            # check that f_eps and f_uu are not in the range of f_u
                            if (jac_w_vals_rank < aug_eps_rank) and (jac_w_vals_rank < aug_Duu_rank):
                                return True
                            else:
                                #print("hi 0")
                                return False
                        else:
                            #print("hi 1")
                            return False
                    else:
                        #print("hi 2")
                        return False
                else:
                    #print("hi 3")
                    return False
            else:
                #print("hi 4")
                return False
        else:
            #print("hi 5")
            return False

    def run_greedy_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                       print_lbls_flag=False, auto_parameters=None, plot_labels=None):

        """
        Function for running the greedy numerical continuation and bistability analysis portions of the mass conservation
        approach. This routine uses the initial value of the principal continuation parameter to construct AUTO
        parameters and then tests varying fixed step sizes for the continuation problem. Note that this routine may
        produce jagged or missing sections in the plots provided. To produce better plots one should use the information
        provided by this routine to run :func:`crnt4sbml.GeneralApproach.run_continuity_analysis`.

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
            plot_labels: list of strings
                A list of strings defining the labels for the x-axis, y-axis, and title. Where the first element
                is the label for x-axis, second is the y-axis label, and the last element is the title label. If
                you would like to use the default settings for some of the labels, simply provide None for that
                element.
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
        See .
        """

        # setting default values for AUTO
        if 'NMX' not in auto_parameters.keys():
            auto_parameters['NMX'] = 10000

        if 'ITMX' not in auto_parameters.keys():
            auto_parameters['ITMX'] = 100

        # making the directory if it doesn't exist
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

        species_num = [str(i) for i in self.__indp_species].index(species) + 1

        species_y = str(self.__indp_species[species_num-1])

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_greedy_continuity_analysis \
            (species_num, parameters, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters, plot_labels)

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    def run_mpi_greedy_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                           print_lbls_flag=False, auto_parameters=None, plot_labels=None):

        """
        Function for running the greedy numerical continuation and bistability analysis portions of the mass conservation
        approach. This routine uses the initial value of the principal continuation parameter to construct AUTO
        parameters and then tests varying fixed step sizes for the continuation problem. Note that this routine may
        produce jagged or missing sections in the plots provided. To produce better plots one should use the information
        provided by this routine to run :func:`crnt4sbml.GeneralApproach.run_continuity_analysis`.

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
            plot_labels: list of strings
                A list of strings defining the labels for the x-axis, y-axis, and title. Where the first element
                is the label for x-axis, second is the y-axis label, and the last element is the title label. If
                you would like to use the default settings for some of the labels, simply provide None for that
                element.
        Returns
        ---------
            multistable_param_ind: list of integers
                A list of those indices in 'parameters' that produce multistable plots.
            sample_portion: list of 1D numpy arrays
                A list of 1D numpy arrays corresponding to those values in the input variable parameters that was
                distributed to the core.
            plot_specifications: list of lists
                A list whose elements correspond to the plot specifications of each element in multistable_param_ind.
                Each element is a list where the first element specifies the range used for the x-axis, the second
                element is the range for the y-axis, and the last element provides the x-y values and special point label
                for each special point in the plot.

        Example
        ---------
        See .
        """

        # setting default values for AUTO
        if 'NMX' not in auto_parameters.keys():
            auto_parameters['NMX'] = 10000

        if 'ITMX' not in auto_parameters.keys():
            auto_parameters['ITMX'] = 100

        # making the directory if it doesn't exist
        if not os.path.isdir(dir_path) and self.__my_rank == 0:
            os.mkdir(dir_path)

        species_num = [str(i) for i in self.__indp_species].index(species) + 1

        species_y = str(self.__indp_species[species_num-1])

        sample_portion = self.__distribute_list_of_points(parameters)
        self.__comm.Barrier()

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_mpi_greedy_continuity_analysis \
            (species_num, sample_portion, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters, plot_labels, self.__my_rank)

        self.__important_info += important_info

        return multistable_param_ind, sample_portion, plot_specifications

    def __distribute_list_of_points(self, samples):

        if self.__my_rank == 0:

            # number of tasks per core
            tasks = len(samples) // self.__num_cores  # // calculates the floor

            # remainder
            r = len(samples) - self.__num_cores * tasks

            # array that holds how many tasks each core has
            tasks_core = numpy.zeros(self.__num_cores, dtype=numpy.int64)
            tasks_core.fill(tasks)

            # distributing in the remainder
            ii = 0
            while r > 0:
                tasks_core[ii] += 1
                r -= 1
                ii += 1

            sample_portion = samples[0:tasks_core[0]]

            if self.__num_cores > 1:
                for i in range(1, self.__num_cores):
                    start = sum(tasks_core[0:i])
                    end = start + tasks_core[i]
                    self.__comm.send(samples[start:end], dest=i, tag=i * 11)

        else:
            if self.__num_cores > 1:
                sample_portion = self.__comm.recv(source=0, tag=self.__my_rank * 11)

        return sample_portion

    def run_mpi_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                    print_lbls_flag=False, auto_parameters=None, plot_labels=None):

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
            plot_labels: list of strings
                A list of strings defining the labels for the x-axis, y-axis, and title. Where the first element
                is the label for x-axis, second is the y-axis label, and the last element is the title label. If
                you would like to use the default settings for some of the labels, simply provide None for that
                element.
        Returns
        ---------
            multistable_param_ind: list of integers
                A list of those indices in 'parameters' that produce multistable plots.
            sample_portion: list of 1D numpy arrays
                A list of 1D numpy arrays corresponding to those values in the input variable parameters that was
                distributed to the core.
            plot_specifications: list of lists
                A list whose elements correspond to the plot specifications of each element in multistable_param_ind.
                Each element is a list where the first element specifies the range used for the x-axis, the second
                element is the range for the y-axis, and the last element provides the x-y values and special point label
                for each special point in the plot.

        Example
        ---------
        See .
        """

        # setting default values for AUTO
        if 'NMX' not in auto_parameters.keys():
            auto_parameters['NMX'] = 10000

        if 'ITMX' not in auto_parameters.keys():
            auto_parameters['ITMX'] = 100

        # making the directory if it doesn't exist
        if not os.path.isdir(dir_path) and self.__my_rank == 0:
            os.mkdir(dir_path)

        species_num = [str(i) for i in self.__indp_species].index(species) + 1

        species_y = str(self.__indp_species[species_num - 1])

        sample_portion = self.__distribute_list_of_points(parameters)
        self.__comm.Barrier()

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_mpi_continuity_analysis \
            (species_num, sample_portion, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters, plot_labels, self.__my_rank)

        self.__important_info += important_info

        return multistable_param_ind, sample_portion, plot_specifications

    def run_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                       print_lbls_flag=False, auto_parameters=None, plot_labels=None):

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
            plot_labels: list of strings
                A list of strings defining the labels for the x-axis, y-axis, and title. Where the first element
                is the label for x-axis, second is the y-axis label, and the last element is the title label. If
                you would like to use the default settings for some of the labels, simply provide None for that
                element.
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
        See .
        """

        # setting default values for AUTO
        if 'NMX' not in auto_parameters.keys():
            auto_parameters['NMX'] = 10000

        if 'ITMX' not in auto_parameters.keys():
            auto_parameters['ITMX'] = 100

        # making the directory if it doesn't exist
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

        species_num = [str(i) for i in self.__indp_species].index(species) + 1

        species_y = str(self.__indp_species[species_num - 1])

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_continuity_analysis \
            (species_num, parameters, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters, plot_labels)

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    def __initialize_ant_string(self, species_num, pcp_x):

        indp_odes_str = []
        indp_odes_str.append(str(self.__indp_system_subs[species_num-1]))
        [indp_odes_str.append(str(self.__indp_system_subs[i])) for i in range(len(self.__indp_species)) if i != species_num-1]

        for i in range(len(indp_odes_str)):
            indp_odes_str[i] = indp_odes_str[i].replace('**', '^')

        ind_spec_conc = []
        ind_spec_conc.append(self.__indp_species[species_num-1])
        [ind_spec_conc.append(str(i)) for i in self.__indp_species if i != self.__indp_species[species_num-1]]

        # building the string of ODEs in Antimony syntax
        ode_str = ''
        for i in range(len(ind_spec_conc)):
            ode_str += 'J' + str(i) + ': -> ' + str(ind_spec_conc[i]) + '; ' + indp_odes_str[i] + ';'

        return ode_str, pcp_x

    def __finalize_ant_string(self, result_x, ode_str):

        vars_to_initialize = [str(i) for i in self.__lambda_variables]
        var_vals = result_x

        ant_str = ode_str
        for i in range(len(vars_to_initialize)):
            ant_str += str(vars_to_initialize[i]) + ' = ' + str(var_vals[i]) + ';'

        return ant_str

    def get_optimization_bounds(self):
        """
        Returns a list of tuples that corresponds to the determined physiological bounds chosen for the problem. Each
        entry corresponds to the list provided by :func:`crnt4sbml.GeneralApproach.get_input_vector`.

         Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_optimization_bounds()
        """

        graph_edges = self.__cgraph.get_g_edges()
        dec_vec_var_def = []
        for i in self.__sympy_reactions + self.__sympy_species:

            if i in self.__sympy_species:
                dec_vec_var_def.append("concentration")
            elif i in self.__sympy_reactions:

                ind = self.__sympy_reactions.index(i)
                reaction = graph_edges[ind]

                reaction_type = self.__cgraph.get_graph().edges[reaction]['type']
                dec_vec_var_def.append(reaction_type)

                if reaction_type is None:
                    output_statement = "The reaction type of reaction " + self.__cgraph.get_graph().edges[reaction][
                        'label'] \
                                       + " could not be identified as it does not fit any biological criteria " + \
                                       "established. \n" + "You must enter bounds manually for this reaction! \n"
                    print(output_statement)

        bounds = [self.get_physiological_range(i) for i in dec_vec_var_def]

        return bounds

    @staticmethod
    def get_physiological_range(for_what=None):
        """
        Molecule concentrations: 5e-13 ... 5e-7 M (basically from 1 per cell to 1e6 per cell)
        Complex formation (bimolecular association): 1e4 ... 1e8  M^-1s^-1
        Example A + B -> [AB]  (or first step in enzymatic catalysis mechanism)
        Complex dissociation: 1e-5 ... 1e-3 s^-1
        Example [AB] -> A + B  (e.g enzyme substrat dissociation)
        Catalysis (equivalent to kcat in Michaelis-Menten): 1e-3 ... 1e0 s^-1
        Examples: [AB] -> A + C  or  X -> Y  (e.g. final step in enzymatic catalysis)
        """
        valid = {"concentration", "complex formation", "complex dissociation", "catalysis"}
        if for_what is None:
            warnings.warn("Please provide the argument what the range should be provided for.")
            return None
        if for_what not in valid:
            raise ValueError("for_what argument must be one of %r." % valid)
        if for_what == "concentration":
            return 5e-13, 5e-7
        if for_what == "complex formation":
            return 1e4, 1e8
        if for_what == "complex dissociation":
            return 1e-5, 1e-3
        if for_what == "catalysis":
            return 1e-3, 1e0

    def generate_report(self):
        """
        Prints out helpful details constructed by :func:`crnt4sbml.GeneralApproach.run_optimization` and
        :func:`crnt4sbml.GeneralApproach.run_continuity_analysis`.

        Example
        --------
        See .
        """
        if self.__comm == None:
            print(self.__important_info)
        else:

            all_important_info = self.__comm.gather(self.__important_info, root=0)
            self.__comm.Barrier()

            if self.__my_rank == 0:

                print("")

                for i in range(1, len(all_important_info)):
                    print(all_important_info[i])
                print(self.__important_info)

    def get_input_vector(self):
        """
        Returns a list of SymPy variables that specifies the ordering of the reactions and species for which bounds
        need to be provided.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_input_vector()
        """
        return self.__sympy_reactions + self.__sympy_species

    def get_decision_vector(self):
        """
        Returns a list of SymPy variables that specifies the ordering of the reactions and species of the decision
        vector to be used in optimization. Note: this method should not be used to create bounds for the optimization
        routine, rather :func:`crnt4sbml.GeneralApproach.get_input_vector` should be used.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_decision_vector()
        """
        return self.__decision_vector

    def get_independent_odes_subs(self):
        """
        Returns a Sympy Matrix representing the independent ODE system with conservation laws substituted in. Each row
        corresponds to the ODE for the species corresponding to the list provided by :func:`crnt4sbml.GeneralApproach.get_independent_species`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_independent_odes_subs()
        """
        return self.__indp_system_subs

    def get_independent_odes(self):
        """
        Returns a Sympy Matrix representing the independent ODE system without conservation laws substituted in. Each row
        corresponds to the ODE for the species corresponding to the list provided by :func:`crnt4sbml.GeneralApproach.get_independent_species`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_independent_odes()
        """
        return self.__indp_system

    def get_independent_species(self):
        """
        Returns a list of SymPy variables that reflects the independent species chosen for the general approach.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_independent_species()
        """

        return self.__indp_species

    def get_fixed_reactions(self):
        """
        Returns a list of SymPy variables that describe the reactions that were chosen to be fixed when ensuring a
        steady-state solution exists.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_fixed_reactions()
        """
        return self.__fixed_reactions

    def get_solutions_to_fixed_reactions(self):
        """
        Returns a list of SymPy expressions corresponding to the fixed reactions. The ordering of the elements
        corresponds to the list returned by :func:`crnt4sbml.GeneralApproach.get_fixed_reactions`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_solutions_to_fixed_reactions()
        """
        return self.__soln_to_fixed_reactions2

    def get_conservation_laws(self):
        """
        Returns a string representation of the conservation laws. Here the values on the left hand side of each equation
        are the constants of the conservation laws.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> print(approach.get_conservation_laws())
        """
        rhs = self.__BT_mat * sympy.Matrix([self.__sympy_species]).T
        laws = ""
        for i in range(rhs.shape[0]):
            laws += 'C' + str(i + 1) + ' = ' + str(rhs[i]) + '\n'

        return laws

    def get_determinant_of_jacobian(self):
        """
        Returns a Sympy expression of the determinant of the Jacobian, where the Jacobian is with respect to the
        independent species.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_determinant_of_jacobian()
        """
        return self.__det_jac_indp_subs

    def get_variables_for_lambda_functions(self):
        """
        Returns a list of all SymPy variables that are used in the independent ODE system with conservation laws
        substituted in, which is given by :func:`crnt4sbml.GeneralApproach.get_independent_odes_subs`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_variables_for_lambda_functions()
        """
        return self.__lambda_variables

    def get_indpendent_odes_lambda_function(self):
        """
        Returns a list of lamda functions where each element corresponds to the row of the independent ODE system with
        conservation laws substituted in, which is given by :func:`crnt4sbml.GeneralApproach.get_independent_odes_subs`.
        The variables for the lambda functions are given by :func:`crnt4sbml.GeneralApproach.get_variables_for_lambda_functions`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_indpendent_odes_lambda_function()
        """

        return self.__lambda_indp_odes



