import re


class LowDeficiencyApproach:
    """
    Class for testing the Deficiency Zero and One Theorems.
    """
    def __init__(self, cgraph):
        """
        Initialization of LowDeficiency Approach class.

        See also
        ---------
        crnt4sbml.CRNT.get_low_deficiency_approach()
        """
        self.__satisfies_deficiency_zero_theorem = False
        self.__satisfies_deficiency_one_theorem = False
        self.__satisfies_relaxed_deficiency_one_theorem = False

        # key conditions
        self.__cond1 = cgraph.get_if_cgraph_weakly_reversible()
        cond2 = cgraph.get_deficiency() == 0
        lcd = cgraph.get_linkage_classes_deficiencies()
        cond3 = all([i < 2 for i in lcd])
        cond4 = sum(lcd) == cgraph.get_deficiency()
        tlc = cgraph.get_number_of_terminal_strong_lc_per_lc()
        cond5 = all([i < 2 for i in tlc])

        # decision tree for low deficiency theorems
        if cond2:
            # deficiency zero theorem
            self.__satisfies_deficiency_zero_theorem = True
        elif self.__cond1 and cond3 and cond4:
            # deficiency one theorem
            self.__satisfies_deficiency_one_theorem = True
        elif cond5 and cond3 and cond4:
            # deficiency one theorem with relaxed criteria
            self.__satisfies_relaxed_deficiency_one_theorem = True

    def does_satisfy_any_low_deficiency_theorem(self):
        """
        Function to see if the network satisfies the Deficiency Zero or One Theorem. Returns True if the network
        satisfies the Deficiency Zero or One Theorem, False otherwise.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_low_deficiency_approach()
        >>> print(approach.does_satisfy_any_low_deficiency_theorem())
            False
        """
        return (self.__satisfies_deficiency_zero_theorem or
                self.__satisfies_deficiency_one_theorem or
                self.__satisfies_relaxed_deficiency_one_theorem)

    def does_satisfy_deficiency_zero_theorem(self):
        """
        Function to see if the network satisfies the Deficiency Zero Theorem. Returns True if the network
        satisfies the Deficiency Zero Theorem, False otherwise.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_low_deficiency_approach()
        >>> print(approach.does_satisfy_deficiency_zero_theorem())
            False
        """
        return self.__satisfies_deficiency_zero_theorem

    def does_satisfy_deficiency_one_theorem(self):
        """
        Function to see if the network satisfies the Deficiency One Theorem. Returns True if the network
        satisfies the Deficiency One Theorem, False otherwise.
        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_low_deficiency_approach()
        >>> print(approach.does_satisfy_deficiency_one_theorem())
            False
        """
        return (self.__satisfies_deficiency_one_theorem or
                self.__satisfies_relaxed_deficiency_one_theorem)

    def report_deficiency_zero_theorem(self):
        """
        Prints out the applicability of the Deficiency Zero Theorem for the provided network.
        Possible output:

        "By the Deficiency Zero Theorem, the differential equations
        cannot admit a positive equilibrium or a cyclic composition
        trajectory containing a positive composition. Thus, multiple
        equilibria cannot exist for the network."

        or

        "By the Deficiency Zero Theorem, there exists within each positive
        stoichiometric compatibility class precisely one equilibrium.
        Thus, multiple equilibria cannot exist for the network."

        or

        "The network does not satisfy the Deficiency Zero Theorem, multistability cannot be excluded."

        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_low_deficiency_approach()
        >>> print(approach.report_deficiency_zero_theorem())
            The network does not satisfy the Deficiency Zero Theorem, multistability cannot be excluded.
        """
        if self.__satisfies_deficiency_zero_theorem and not self.__cond1:
            report = ("""
            By the Deficiency Zero Theorem, the differential equations 
            cannot admit a positive equilibrium or a cyclic composition 
            trajectory containing a positive composition. Thus, multiple 
            equilibria cannot exist for the network. 
            """)
            print(re.sub(r"^\s+", "", report, flags=re.MULTILINE))
        elif self.__satisfies_deficiency_zero_theorem and self.__cond1:
            report = (""" 
            By the Deficiency Zero Theorem, there exists within each positive
            stoichiometric compatibility class precisely one equilibrium.
            Thus, multiple equilibria cannot exist for the network.
            """)
            print(re.sub(r"^\s+", "", report, flags=re.MULTILINE))
        else:
            print("The network does not satisfy the Deficiency Zero Theorem, multistability cannot be excluded.")

    def report_deficiency_one_theorem(self):
        """
        Prints out the applicability of the Deficiency One Theorem for the provided network.
        Possible output:

        "By the Deficiency One Theorem, the differential equations
        admit precisely one equilibrium in each positive stoichiometric
        compatibility class. Thus, multiple equilibria cannot exist
        for the network."

        or

        "The network satisfies relaxed Deficiency One Theorem. That is it
        is not weakly reversable, but each linkage class contains no more
        than one terminal linkage class. There can exist within a positive
        stoichiometric compatibility class at most one equilibrium.
        Thus, multiple equilibria cannot exist for the network."

        or

        "The network does not satisfy the Deficiency One Theorem, multistability cannot be excluded."

        :download:`Fig1Ci.xml <../../sbml_files/Fig1Ci.xml>` for the provided
        example.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")
        >>> approach = network.get_low_deficiency_approach()
        >>> print(approach.report_deficiency_zero_theorem())
            The network does not satisfy the Deficiency One Theorem, multistability cannot be excluded.
        """
        if self.__satisfies_deficiency_one_theorem:
            report = ("""
            By the Deficiency One Theorem, the differential equations
            admit precisely one equilibrium in each positive stoichiometric
            compatibility class. Thus, multiple equilibria cannot exist 
            for the network.
            """)
            print(re.sub(r"^\s+", "", report, flags=re.MULTILINE))
        elif self.__satisfies_relaxed_deficiency_one_theorem:
            report = ("""
            The network satisfies relaxed Deficiency One Theorem. That is it 
            is not weakly reversable, but each linkage class contains no more 
            than one terminal linkage class. There can exist within a positive 
            stoichiometric compatibility class at most one equilibrium.
            Thus, multiple equilibria cannot exist for the network.
            """)
            print(re.sub(r"^\s+", "", report, flags=re.MULTILINE))
        else:
            print("The network does not satisfy the Deficiency One Theorem, multistability cannot be excluded.")
