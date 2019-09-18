.. _physio-bnds-label:

Creating Physiological Bounds
===============================

To make crnt4sbml more user friendly and make its search limited to physiological problems, we have constructed the
functions :meth:`crnt4sbml.MassConservationApproach.get_optimization_bounds` and
:meth:`crnt4sbml.SemiDiffusiveApproach.get_optimization_bounds`, which constructs the appropriate bounds that must be
provided to the mass conservation and semi-diffusive optimization routines, respectively. Although this feature can be
extremely useful especially if the user is continually changing the SBML file, it should be used with some amount of
caution.

+++++++++++++++++++++++++++
Assigning Reaction Types
+++++++++++++++++++++++++++

To provide these physiological bounds, crnt4sbml first identifies the reactions of the network. A reaction can be
identified as complex formation, complex dissociation, or catalysis, no other type of reaction is considered. To
make this assignment, the reactants and products of the reaction along with the associated stoichiometries are
found for the particular reaction. Using the sum of the stoichiometries for the reactants and products the decision
tree below is used to determine the type of reaction. If the reaction is not identified as complex formation, complex
dissociation, or catalysis, then an error message will be provided and the reaction type will be specified as "None".

.. image:: ./images_for_docs/decision_tree_for_reaction_type.png
   :width: 680px
   :align: center
   :height: 420px


++++++++++++++++++++++++++++
Mass Conservation Approach
++++++++++++++++++++++++++++

Under Development



+++++++++++++++++++++++++++
Semi-diffusive Approach
+++++++++++++++++++++++++++

Under Development