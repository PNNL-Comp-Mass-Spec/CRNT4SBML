.. _my-basic-label:

========================
Low Deficiency Approach
========================

Now that we have constructed the SBML file using the guidelines of :ref:`my-celldesigner-label`, we will proceed by
testing the Deficiency Zero and One Theorems of :cite:`fein_lecture`. We will complete this test for
:download:`Fig1Ci.xml <../sbml_files/Fig1Ci.xml>`. The first step we must
take is importing crnt4sbml. To do this open up a python script and add the following line:

.. code-block:: python

        import crnt4sbml

Next, we will take the SBML file created using CellDesigner and import it into the code. This is done by instantiating
CRNT with a string representation of the path to the SBML file. An example of this instantiation is as follows:

.. code-block:: python

        network = crnt4sbml.CRNT("/path/to/Fig1Ci.xml")

Once this line is ran the class CRNT takes the SBML file and parses it into a Python
`NetworkX <https://networkx.github.io/documentation/stable/>`_ object which is then used to
identify the basic Chemical Reaction Network Theory properties of the network. To obtain a full list of what is provided
by this instantiation, please see the getter methods of :meth:`crnt4sbml.CRNT`. To obtain a print out of the
number of species, complexes, reactions and deficiency of the network complete the following command:

.. code-block:: python

        network.basic_report()

For the closed portion of the C-graph the output should be as follows::

        Number of species: 7
        Number of complexes: 9
        Number of reactions: 9
        Network deficiency: 2

It is important for the user to verify that the number of species, complexes, reactions, and if possible deficiency
values are correct at this stage. To provide another check to make sure the parsing and CellDesigner model were
constructed correctly, one is encouraged to print the network constructed. To do this, add the following command
to the script:

.. code-block:: python

        network.print_c_graph()

After running this command for the constructed SBML file, the following output is obtained.

::

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

Notice that this output describes the reactions in terms of the species' id and not the species' name. Along with the
reactions, the reaction labels constructed during parsing are also returned. For this example the first reaction
s1+s2 -> s3 has a reaction label of 're1' and the reaction s15 -> s1+s6 has a reaction label of 're5r'.  Please note
that the species id and reaction labels may be different if the user has constructed the SBML file themselves. Further
information of the network can be found by analyzing the getter methods of :meth:`crnt4sbml.Cgraph`.

Once one has verified that the network and CellDesigner model were created correctly, we can begin to check the
properties of the network. If one is only interested in whether or not the network precludes bistability, it is best to
first check the Deficiency Zero and One Theorems of Chemical Reaction Network Theory. To do this add the following lines
to the script:

.. code-block:: python

        ldt = network.get_low_deficiency_approach()
        ldt.report_deficiency_zero_theorem()
        ldt.report_deficiency_one_theorem()

This provides the following output for the closed portion of the C-graph::

        The network does not satisfy the Deficiency Zero Theorem, multistability cannot be excluded.
        The network does not satisfy the Deficiency One Theorem, multistability cannot be excluded.

For information on the possible output for this run, please see :func:`crnt4sbml.LowDeficiencyApproach.report_deficiency_one_theorem`
and :func:`crnt4sbml.LowDeficiencyApproach.report_deficiency_zero_theorem`.
