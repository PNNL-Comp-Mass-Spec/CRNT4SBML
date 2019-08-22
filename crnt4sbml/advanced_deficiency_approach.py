import sympy
import sys
import numpy
import networkx

mach_eps = numpy.finfo(float).eps


class AdvancedDeficiencyApproach:

    def __init__(self, cgraph):
        self.__cgraph = cgraph

        # construct dictionary for those edges that are reversible
        temp = [[e[0].split('+'), e[1].split('+')] for e in self.__cgraph.get_g_edges()]
        edges = [e for e in self.__cgraph.get_g_edges()]
        for i in temp:
            for j in i:
                j.sort()

        self.__reversible_dict = {}
        for i in range(len(edges)):

            if [temp[i][1], temp[i][0]] in temp:
                self.__reversible_dict.update({edges[i]: True})
            else:
                self.__reversible_dict.update({edges[i]: False})

        self.__Y = self.__cgraph.get_y()

        self.__species = self.__cgraph.get_species()

        self.__dict_id_name = self.__cgraph.get_dict_id_name()

    def __create_s_matrix_nullspace(self):
        S = sympy.zeros(len(self.__species), len(self.__orientation_graph_edges))
        original_nodes = self.__cgraph.get_g_nodes()
        for i in range(len(self.__orientation_graph_edges)):
            edge = self.__orientation_graph_edges[i]
            target = original_nodes.index(edge[1])
            source = original_nodes.index(edge[0])
            S[:, i] = self.__Y[:, target] - self.__Y[:, source]


        print("len(species)")
        print(len(self.__species))
        print("S")
        sympy.pprint(S)

        the_temp_null = S.nullspace()

        self.__d = len(the_temp_null)

        SN_row = len(self.__orientation_graph_edges)
        SN_col = len(the_temp_null)

        the_S_null = []
        for j in range(SN_row):
            the_S_null.append([the_temp_null[i][j] for i in range(SN_col)])

        return the_S_null, SN_row, SN_col

    def __initialize_orientation(self):

        c_graph = self.__cgraph.get_graph()

        pruned_edges = []
        for edge in self.__cgraph.get_g_edges():
            if (edge[1], edge[0]) not in pruned_edges:
                pruned_edges.append(edge)

        self.__orientation_graph = c_graph.edge_subgraph(pruned_edges).copy()
        self.__orientation_graph_edges = pruned_edges

    def __is_scalar_multiple_greater_than_zero(self, a, b):

        # get entries in a and b that have zero entry
        a_ind = [i for i in range(len(a)) if abs(a[i]) <= mach_eps]
        b_ind = [i for i in range(len(b)) if abs(b[i]) <= mach_eps]

        if a_ind == b_ind:
            vals = [b[i] / a[i] for i in range(len(a)) if i not in a_ind]
            if vals[0] > mach_eps and vals.count(vals[0]) == abs(len(a) - len(a_ind)):
                return True
            else:
                return False
        else:
            return False

    def __find_equivalence_classes(self):

        self.__the_S_null, SN_row, SN_col = self.__create_s_matrix_nullspace()

        self.__equivalence_classes = []
        self.__equivalence_classes_edges = []
        temp_list = [0.0 for i in range(SN_col)]

        # forming P_0 equivalence class
        temp_P_0 = []
        for i in range(SN_row):
            if all([abs(self.__the_S_null[i][j] - temp_list[j]) <= mach_eps for j in range(SN_col)]):
                temp_P_0.append(self.__orientation_graph_edges[i])

        if temp_P_0:
            self.__equivalence_classes.append(self.__orientation_graph.edge_subgraph(temp_P_0).copy())
            self.__equivalence_classes_edges.append(temp_P_0)
        else:
            self.__equivalence_classes.append("Empty")
            self.__equivalence_classes_edges.append("Empty")

        # getting dependent rows of Ker(L_orientation)
        # to create nonzero equivalent classes
        self.__dep_rows_S_null = []
        loop_iter = [i for i in range(SN_row) if i not in temp_P_0]
        for i in loop_iter:
            temp = [self.__orientation_graph_edges[i]]
            for j in loop_iter:
                if j != i:
                    temp_mat = sympy.Matrix([self.__the_S_null[i], self.__the_S_null[j]])
                    _, vals = temp_mat.T.rref()

                    if len(vals) == 1:
                        temp.append(self.__orientation_graph_edges[j])

                temp.sort()
                if temp not in self.__dep_rows_S_null and len(temp) != 1:
                    self.__dep_rows_S_null.append(temp)

        # creating nonzero equivalent classes
        [self.__equivalence_classes.append(self.__orientation_graph.edge_subgraph(self.__dep_rows_S_null[i]).copy())
         for i in range(len(self.__dep_rows_S_null))]

        [self.__equivalence_classes_edges.append(self.__dep_rows_S_null[i]) for i in range(len(self.__dep_rows_S_null))]

        #print("equivalence_classes_edges")
        #for i in self.__equivalence_classes_edges:
        #    print(i)

        #print("")

    def __equivalence_classes_check(self):

        # check that all reactions in P_0 are reversible
        if self.__equivalence_classes[0] == "Empty":
            a_chk = True
        else:
            is_reversible = []
            for i in self.__equivalence_classes_edges[0]:
                is_reversible.append(self.__reversible_dict[i])
            a_chk = all(is_reversible)

        b_chk = []
        for i in range(1, len(self.__equivalence_classes)):

            temp = []
            for j in self.__equivalence_classes_edges[i]:
                temp.append(self.__reversible_dict[j])

            temp_f_ind = [ii for ii in range(len(temp)) if temp[ii] is False]

            if len(temp_f_ind) > 1:
                temp_scalar = []
                new_iter = [self.__orientation_graph_edges.index(self.__dep_rows_S_null[i - 1][k]) for k in
                            range(len(self.__dep_rows_S_null[i - 1])) if k in temp_f_ind]

                for ii in new_iter:
                    temp_scalar += [
                        self.__is_scalar_multiple_greater_than_zero(self.__the_S_null[ii], self.__the_S_null[jj]) for jj
                        in new_iter if jj != ii]

                b_chk.append(all(temp_scalar))

            else:
                b_chk.append(True)

        if all(b_chk) and a_chk:
            return True
        else:
            return False

    def __find_fundamental_classes(self):

        self.__fundamental_classes = []
        self.__fundamental_classes_edges = []
        for i in range(len(self.__equivalence_classes)):

            if self.__equivalence_classes[i] == "Empty":
                self.__fundamental_classes.append("Empty")
                self.__fundamental_classes_edges.append("Empty")
            else:
                self.__fundamental_classes.append(self.__equivalence_classes[i].copy())
                self.__fundamental_classes_edges.append(self.__equivalence_classes_edges[i][:])

                for j in range(len(self.__fundamental_classes_edges[i])):
                    if self.__reversible_dict[self.__fundamental_classes_edges[i][j]]:
                        self.__fundamental_classes[i].add_edge(self.__fundamental_classes_edges[i][j][1],
                                                               self.__fundamental_classes_edges[i][j][0])
                        self.__fundamental_classes_edges[i].append((self.__fundamental_classes_edges[i][j][1],
                                                                    self.__fundamental_classes_edges[i][j][0]))

    def __find_colinkage_sets(self):

        fund_len = len(self.__fundamental_classes)

        if self.__fundamental_classes[0] == "Empty":
            start = 1
        else:
            start = 0

        colinkage_sets = []
        strong_colinkage_sets = []
        for fund_i in range(start, fund_len):

            linkages = [self.__fundamental_classes[fund_i].subgraph(c) for c in
                        networkx.weakly_connected_components(self.__fundamental_classes[fund_i])]

            for lc_i in linkages:

                colinkage_sets.append([i for i in self.__cgraph.get_g_nodes() if i in lc_i.nodes])
                strong_lcs = [lc_i.subgraph(c).copy() for c in networkx.strongly_connected_components(lc_i)]

                temp2 = []
                for slc_j in strong_lcs:

                    slc_j_lbls = [i for i in self.__cgraph.get_g_nodes() if i in slc_j.nodes]

                    x = [lc_i.successors(v_i) for v_i in slc_j_lbls]
                    connected_out_vertices = [item for sublist in x for item in sublist]
                    d = set(connected_out_vertices) - set(slc_j_lbls)

                    if d == set():
                        temp2.append([slc_j_lbls, "T"])
                    else:
                        temp2.append([slc_j_lbls, "NT"])

                strong_colinkage_sets.append(temp2)

        return strong_colinkage_sets, colinkage_sets

    def __find_representative_set_w(self):

        self.__representative_set_W = []
        for P_i in self.__equivalence_classes_edges:
            if P_i != "Empty":
                P_i_reversibility = [self.__reversible_dict[i] for i in P_i]
                if all(P_i_reversibility):
                    self.__representative_set_W.append(P_i[0])
                else:
                    False_ind = P_i_reversibility.index(False)
                    self.__representative_set_W.append(P_i[False_ind])
            else:
                self.__representative_set_W.append("Empty")

    def __realign_orientation(self):

        #print("hi")

        #print(self.__orientation_graph_edges)

        #print("")
        realignment_flag = False
        for i in range(len(self.__representative_set_W)):
            if self.__representative_set_W[i] != "Empty":
                y_i = self.__representative_set_W[i]
                v_i = self.__orientation_graph_edges.index(y_i)
                for j in range(len(self.__equivalence_classes_edges[i])):
                    v_y = self.__orientation_graph_edges.index(self.__equivalence_classes_edges[i][j])
                    out = self.__is_scalar_multiple_greater_than_zero(self.__the_S_null[v_i], self.__the_S_null[v_y])
                    if not out:
                        print("attempting to realign orientation look into this!")
                        self.__orientation_graph.remove_edge(*self.__equivalence_classes_edges[i][j])
                        self.__orientation_graph_edges.remove(self.__equivalence_classes_edges[i][j])
                        self.__orientation_graph.add_edge(self.__equivalence_classes_edges[i][j][1],
                                                          self.__equivalence_classes_edges[i][j][0])
                        self.__orientation_graph_edges.append((self.__equivalence_classes_edges[i][j][1],
                                                               self.__equivalence_classes_edges[i][j][0]))

                        self.g_print(self.__cgraph.get_graph(), self.__orientation_graph_edges)

                        realignment_flag = True

                        sys.exit()

        if realignment_flag:
            self.__find_equivalence_classes()
            if not self.__equivalence_classes_check():
                "Network does not have the capacity for multiple steady states."
                sys.exit()

            # testing the realignment
            for i in range(len(self.__representative_set_W)):
                if self.__representative_set_W[i] != "Empty":
                    y_i = self.__representative_set_W[i]
                    v_i = self.__orientation_graph_edges.index(y_i)
                    for j in range(len(self.__equivalence_classes_edges[i])):
                        v_y = self.__orientation_graph_edges.index(self.__equivalence_classes_edges[i][j])
                        out = self.__is_scalar_multiple_greater_than_zero(self.__the_S_null[v_i], self.__the_S_null[v_y])
                        if not out:
                            print("After realignment criteria in step 5 not satisfied!")
                            sys.exit()

    def __create_basis_ker_L_O_intersect_Gam_W(self):

        print("self.__orientation_graph_edges")
        print(self.__orientation_graph_edges)

        temp_representative = self.__orientation_graph_edges

        indices = [i for i in range(len(self.__orientation_graph_edges)) if self.__orientation_graph_edges[i]
                   in self.__representative_set_W]

        print("indices")
        print(indices)

        basis_Gam_W = sympy.zeros(len(self.__orientation_graph_edges), len(indices))
        count = 0
        for i in indices:
            basis_Gam_W[i, count] = 1
            count += 1

        #print("basis_Gam_W")
        #sympy.pprint(basis_Gam_W)

        #sys.exit()

        #print("temp_representative")
        #print(temp_representative)
        #empty_indices = [i for i in range(len(temp_representative)) if temp_representative[i] == 'Empty']

        #remove empty entries
        #for i in reversed(empty_indices):
        #    del temp_representative[i]

        #del temp_representative[1]
        #print(temp_representative)

        S = sympy.zeros(len(self.__species), len(temp_representative))
        original_nodes = self.__cgraph.get_g_nodes()
        for i in range(len(temp_representative)):
            edge = temp_representative[i]
            target = original_nodes.index(edge[1])
            source = original_nodes.index(edge[0])
            S[:, i] = self.__Y[:, target] - self.__Y[:, source]


        #print("self.__the_S_null")
        #sympy.pprint(self.__the_S_null)


        #print("S")
        #sympy.pprint(S)

        #find a basis for R(L_O^T) i.e. column space of L_O^T
        S_T = S.T
        #print("S_T")
        #sympy.pprint(S_T)
        rref_matrix, rref_pivots = S_T.rref()

        #print("rref_pivots")
        #sympy.pprint(rref_pivots)

        basis_column_space_L_O_T = sympy.zeros(len(temp_representative), len(rref_pivots))
        count = 0
        for i in rref_pivots:
            #print(i)
            basis_column_space_L_O_T[:, count] = S_T[:, i]
            count += 1

        #print("basis_column_space_L_O_T")
        #sympy.pprint(basis_column_space_L_O_T)

        A = basis_column_space_L_O_T.row_join(-1*basis_Gam_W)

        #print("")

        #print("A")
        #sympy.pprint(A)

        #print("nullspace")
        A_null = A.nullspace()
        #print(A_null)

        if len(A_null) != ((len(self.__representative_set_W) - 1) - self.__d):
            print("q != w - d")
            sys.exit()

        #print("")
        #print(A_null[0])

        A_null_mat = sympy.zeros(len(rref_pivots), len(A_null))

        for i in range(len(A_null)):
            for j in range(len(rref_pivots)):
                A_null_mat[j, i] = A_null[i][j]

        #print("A_mat_null")
        #sympy.pprint(A_null_mat)

        #print("")
        #print("basis_column_space_L_O_T*A_null_mat")
        full_basis = basis_column_space_L_O_T*A_null_mat
        #sympy.pprint(full_basis)
        #print(full_basis.shape[1])
        full_basis_W = sympy.zeros(len(indices), full_basis.shape[1])
        #print("")
        count = 0
        for i in indices:
            full_basis_W[count, :] = full_basis[i, :]
            count += 1

        #temp_full_basis = full_basis_W[5,:].col_join(full_basis_W[0:5,:])

        #print("")
        #sympy.pprint(full_basis_W)
        #sympy.pprint(temp_full_basis)

        #sympy.pprint(sympy.Matrix([[0], [1], [-1], [0], [1], [-1]]))

        #temp = (temp_full_basis.row_join(sympy.Matrix([[1], [0], [0], [-1], [1], [0]]))).rref()

        #sympy.pprint(temp)


        return full_basis_W, indices


    def __forest_basis_exists(self, full_basis_W, indices):

        print("full_basis_W")
        sympy.pprint(full_basis_W)
        print("full basis W rows")
        basis_rows = [self.__orientation_graph_edges[i][0] + '->' + self.__orientation_graph_edges[i][1]
                      for i in indices]

        basis_columns = ['a' + str(i) for i in range(full_basis_W.shape[1])]

        print("basis_rows")
        print(basis_rows)
        print("basis_columns")
        print(basis_columns)
        edges = []
        for i in range(full_basis_W.shape[0]):
            for j in range(full_basis_W.shape[1]):
                if abs(full_basis_W[i,j]) > mach_eps:
                    edges.append((basis_rows[i], basis_columns[j]))

        basis_graph = networkx.DiGraph()

        print("edges")
        print(edges)

        basis_graph.add_edges_from(edges)

        print(networkx.is_forest(basis_graph))


    def change_name(self, edge):

        splits = edge.split('+')
        switched = [self.__dict_id_name[i] for i in splits]
        if len(switched) == 1:
            return switched[0]
        else:
            return switched[0] + '+' + switched[1]

    def g_print(self, graph, graph_edges):
        print("graph_edges")
        print(graph_edges)
        for e in graph_edges:
            if e == 'Empty':
                print(e)
            else:
                print(self.change_name(e[0]) + ' -> ' + self.change_name(e[1]) + '  --  ' + graph.edges[e]['label'])
        print("")

    def run_higher_deficiency_algorithm(self):

        #self.__cgraph.print()

        self.__initialize_orientation()
        self.g_print(self.__cgraph.get_graph(), self.__orientation_graph_edges)

        self.__find_equivalence_classes()

        #print("equivalence_classes")
        #for i in range(len(self.__equivalence_classes)):
        #    if i != 0:
        #        self.g_print(self.__equivalence_classes[i], self.__equivalence_classes_edges[i])
        #print("")

        if self.__equivalence_classes_check():
            self.__find_fundamental_classes()

            #print("")
            #print("fundamental classes")
            #for i in range(1, len(self.__fundamental_classes)):
            #    self.g_print(self.__cgraph.get_graph(), self.__fundamental_classes_edges[i])
            #print("")

            self.__strong_colinkage_sets, self.__colinkage_sets = self.__find_colinkage_sets()

            #print("strong_colinkage_sets")
            #for i in strong_colinkage_sets:
            #    for j in i:
            #        print(j)
            #    print("")
            #print("")

            #print("colinkage_sets")
            #for i in colinkage_sets:
            #    print(i)
            #print("")

            self.__find_representative_set_w()

            print("self.__representative_set_W")
            print(self.__representative_set_W)

            self.__realign_orientation()

            full_basis_W, indices = self.__create_basis_ker_L_O_intersect_Gam_W()

            self.__forest_basis_exists(full_basis_W, indices)

        else:
            "Network does not have the capacity for multiple steady states."
            sys.exit()

    def report(self):
        print("placeholder")
