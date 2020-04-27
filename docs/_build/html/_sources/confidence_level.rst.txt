.. _confidence-level-label:

===========================
Confidence Level Routine
===========================

Although the general, mass conservation, and semi-diffusive approaches can quickly provide confirmation of bistability
for most examples, this may not always be the case. In fact, an important item of discussion is that these approaches
cannot exclude bistability, even if a large amount of random decision vectors are explored. It is this uncertainty that
we wish to address. This is done by assigning a probability that the minimum objective function value achieved is equal
to the true global minimum. We achieve this probability by considering a slightly modified version of the unified
Bayesian stopping rule in :cite:`BOLTON2004549` and Theorem 4.1 of :cite:`snyman_bayesian`, where the rule was first
established.

Let :math:`\alpha_k` and :math:`\alpha^*` denote the probability that the optimization routine has converged to the
local minimum objective function value, say :math:`f_k`, and global minimum objective function value, say :math:`f^*`.
Assuming that :math:`\alpha^* \geq \alpha_k` for all local minimum values :math:`f_k` we may then state that the
probability that :math:`\tilde{f} = f^*` is as follows:

:math:`Pr[\tilde{f} = f^*] \geq q(n, r) = 1 - \dfrac{(n + a + b - 1)! (2n + b - r - 1)!}{(2n + a + b - 1)! (n + b-r-1)!}`,

where :math:`n` is the number of initial decision vectors that are considered, :math:`\tilde{f} = min \{f_1, \dots, f_n \}`
, :math:`a` and :math:`b` are parameters of the Beta distribution :math:`\beta(a, b)`, and :math:`q(n, r)` is the
confidence level. We then let :math:`r` be the number of :math:`f_k` for :math:`k = 1, \dots, n` that are in the
neighborhood of :math:`\tilde{f}`.

Given our minimum objective function value is zero, for some networks it may be the case that the :math:`f_k` are nearly
zero with respect to machine precision. For this reason, we say that :math:`f_k` is in the neighborhood of :math:`\tilde{f}`
if

:math:`\dfrac{| \tilde{f} - f_k |}{\tilde{f}} \leq 10^{-2}`.

This means that :math:`f_k` is in the neighborhood of :math:`\tilde{f}` if the relative error of :math:`f_k` and :math:`\tilde{f}`
is less than 1%. If :math:`\tilde{f}` is considered zero with respect to the system's minimum positive normalized float, then we
consider this value zero and provide :math:`q(n,r) = 1.0`, skipping the computation of :math:`q(n,r)`. Thus, we can state that
the probability that the obtained :math:`\tilde{f}` is the global minimum (for the prescribed bounds of the decision vector)
is greater than or equal to the confidence level :math:`q(n, r)`. Using the standard practice in statistics, it should be
noted that :math:`q(n,r) \geq 0.95` is often considered an acceptable confidence level to make the conclusion that :math:`\tilde{f}`
is the global minimum of the objective function.

For information on how to enable the construction of a confidence level for each of the approaches, please refer to the
following for each approach:

- Mass conservation approach:
    - If using :func:`crnt4sbml.MassConservationApproach.run_optimization` set confidence_level_flag = True and and prescribe a value to change_in_rel_error (if applicable)
    - If using :func:`crnt4sbml.MassConservationApproach.run_mpi_optimization` set confidence_level_flag = True and and prescribe a value to change_in_rel_error (if applicable)
- Semi-diffusive approach:
    - If using :func:`crnt4sbml.SemiDiffusiveApproach.run_optimization` set confidence_level_flag = True and prescribe a value to change_in_rel_error (if applicable)
    - If using :func:`crnt4sbml.SemiDiffusiveApproach.run_mpi_optimization` set confidence_level_flag = True and and prescribe a value to change_in_rel_error (if applicable)
- General approach:
    - If using :func:`crnt4sbml.GeneralApproach.run_optimization` set confidence_level_flag = True and prescribe a value to change_in_rel_error (if applicable)
