===========
Quick Start
===========

To begin using CRNT4SBML, start by following the process outlined in :ref:`my-installation-label`. Once you have
correctly installed CRNT4SBML follow the steps below to obtain a general idea of how one can perform the mass conservation
and semi-diffusive approach of :cite:`irene`.

- If you are interested in running the Deficiency Zero and One theorems please consult :ref:`my-basic-label`.
- If one is interested in the general steps to follow in order to detect bistability, one should consult :ref:`detect-bistability-label`.

.. _quickstart-deficiency-label:

++++++++++++++++++++++++++++++++++++
Mass Conservation Approach Example
++++++++++++++++++++++++++++++++++++ 

In order to run the mass conservation approach one needs to first create an SBML file of the reaction network. The
SBML file representing the reaction network for this example is given by :download:`Fig1Ci.xml <../sbml_files/Fig1Ci.xml>`.
It is highly encouraged that the user consult :ref:`my-celldesigner-label` when considering their own individual network
as the format of the SBML file must follow a certain construction to be easily used by CRNT4SBML.

To run the mass conservation approach create the following python script:

.. code-block:: python

   import crnt4sbml

   network = crnt4sbml.CRNT("/path/to/Fig1Ci.xml")

   opt = network.get_mass_conservation_approach()

   bounds, concentration_bounds = opt.get_optimization_bounds()

   params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                        concentration_bounds=concentration_bounds)

   multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
                                                                                   auto_parameters={'PrincipalContinuationParameter': 'C3'})

   opt.generate_report()

This will provide the following output along with creating the directory "num\_cont\_graphs" in your current
directory that contains multistability plots. Please note that runtimes and the number of multistability plots produced
may vary among different operating systems. Please see :ref:`my-deficiency-label` for a more detailed explanation of
running the mass conservation approach and the provided output.

::

    Creating Equilibrium Manifold ...
    Elapsed time for creating Equilibrium Manifold: 1.9408779999999997

    Running feasible point method for 10 iterations ...
    Elapsed time for feasible point method: 2.790184

    Running the multistart optimization ...

    Smallest value achieved by objective function: 0.0

    Elapsed time for multistart method: 13.814762

    Running continuity analysis ...
    Elapsed time for continuity analysis: 22.570291996002197

    The number of feasible points that satisfy the constraints: 10
    Total feasible points that give F(x) = 0: 4
    Total number of points that passed final_check: 4
    Number of multistability plots found: 2
    Elements in params_for_global_min that produce multistability:
    [2, 3]

.. _`quickstart-injectivity-label`:

+++++++++++++++++++++++++++++++++++++
Semi-diffusive Approach Example
+++++++++++++++++++++++++++++++++++++

To run the semi-diffusive approach one needs to create the SBML file specific for semi-diffusive networks. The SBML file
representing the reaction network for this example is given by :download:`Fig1Cii.xml <../sbml_files/Fig1Cii.xml>`. It
is highly encouraged that the user consult :ref:`my-celldesigner-label` when considering their own individual network as
the format of the SBML file must follow a certain construction to be easily used by crnt4sbml.

To run the semi-diffusive approach create the following python script:

.. code-block:: python

    import crnt4sbml

    network = crnt4sbml.CRNT("path/to/Fig1Cii.xml")

    opt = network.get_semi_diffusive_approach()

    bounds = opt.get_optimization_bounds()

    params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds)

    multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s7", parameters=params_for_global_min,
                                                                                    auto_parameters={'PrincipalContinuationParameter': 're17'})

    opt.generate_report()


This will provide the following output along with creating the directory "num\_cont\_graphs" in your current
directory that contains multistability plots. Please note that runtimes and the number of multistability plots produced
may vary among different operating systems. Please see :ref:`my-injectivity-label` for a more detailed explanation of
running the semi-diffusive approach and the provided output.

::

    Running feasible point method for 10 iterations ...
    Elapsed time for feasible point method: 0.7991039999999998

    Running the multistart optimization ...

    Smallest value achieved by objective function: 0.0

    Elapsed time for multistart method: 45.470756

    Running continuity analysis ...
    Elapsed time for continuity analysis: 67.03238201141357

    The number of feasible points that satisfy the constraints: 10
    Total feasible points that give F(x) = 0: 9
    Total number of points that passed final_check: 9
    Number of multistability plots found: 9
    Elements in params_for_global_min that produce multistability:
    [0, 1, 2, 3, 4, 5, 6, 7, 8]


