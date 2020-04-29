.. _direct-simulation-label:

===========================================
Direct Simulation for the General Approach
===========================================

When using the general approach it is possible that the optimization routine finds kinetic constants that force the Jacobian
of the system to be ill-conditioned or even singular, even if species concentrations are varied. If this particular scenario
occurs, numerical continuation will not be able to continue as it relies on a well-conditioned Jacobian. To overcome this type
of situation we have constructed the function :func:`crnt4sbml.GeneralApproach.run_direct_simulation` for the general
approach. The direct simulation routine strategically chooses the initial conditions for the ODE system and then simulates
the ODEs until a steady state occurs. Then based on the user defined signal and optimization values provided, it will vary
the signal amount and simulate the ODE system again until a steady state occurs. By varying the signal for several values
and different initial conditions, direct simulation is able to construct a bifurcation diagram. Given the direct simulation
method is numerically integrating the system of ODEs, this method will often take longer than the numerical continuation
routine. Although this is the case, direct simulation may be able to provide a bifurcation diagram when numerical continuation
cannot.

+++++++++++++++++++++++++++++++++++
Failure of numerical continuation
+++++++++++++++++++++++++++++++++++

Under Development ...


++++++++++++++++++++++++++++++++++++++
Outline of direct simulation process
++++++++++++++++++++++++++++++++++++++

Under Development ...

------------------------------------------
Finding appropriate initial conditions
------------------------------------------

Under Development ...

------------------------------------------
Finding a steady state to the system
------------------------------------------

Under Development ...

-------------------------------------------
Conducting forward and reverse simulations
-------------------------------------------

Under Development ...

-------------------------------------------
Varying the signal of the ODE system
-------------------------------------------

Under Development ...

-------------------------------------------
Constructing the bifurcation diagram
-------------------------------------------

Under Development ...





