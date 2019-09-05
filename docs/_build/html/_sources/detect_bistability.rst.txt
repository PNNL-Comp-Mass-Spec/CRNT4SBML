.. _detect-bistability-label:

=================================
Steps for Detecting Bistability
=================================

The following are some simple steps to follow for detecting bistability using CRNT4SBML:

1. Construct an SBML file following the guidelines provided in :ref:`my-celldesigner-label`.
2. Check if the conditions for the Deficiency Zero or One Theorems are satisfied using the approach outlined in :ref:`my-basic-label`.
3. If the Deficiency Zero or One Theorems are not satisfied, then use :func:`crnt4sbml.Cgraph.get_dim_equilibrium_manifold` to find :math:`\lambda`, the number of mass conservation relationships.
4. If :math:`\lambda` is greater than zero one can use the details described in :ref:`my-deficiency-label` to conduct further analysis of the network.
5. If :math:`\lambda` is zero and there is a boundary species present in the SBML file then one can use the details described in :ref:`my-injectivity-label` to conduct further analysis of the network.

