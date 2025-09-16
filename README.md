# Analytical Electrostatics
---------------------------
Copyright (C) 2024-2025 Sergii V. Siryk

Copyright (C) 2024-2025 Vincenzo Di Florio

This software is distributed under the terms
the terms of the GNU/GPL licence v3

# Overview
----------

This repository contains a MATLAB script to compute the analytical solution of the Linearized Poisson-Boltzmann equation:

$$
-\mathrm{div} \left( \varepsilon_0 \varepsilon_r(\mathbf r) \nabla \varphi \right) + \kappa(\mathbf r)^2 \varphi = \rho^f
$$

for solutes modeled as $N_s$ non-overlapping dielectric spheres $\Omega_{m,i}$ ($i=1,\ldots,N_s$) with the same relative dielectric constant $\varepsilon_{r,m}$. Each sphere, centered at $\mathbf r_i\in\mathbb R^3$, has a radius $R_i$ and contains a fixed centrally-located point charge $q_i$.  The goal is to determine the total self-consistent potential  $\phi(\mathbf r)$  at a given point $\mathbf r\in\mathbb R^3$ in the form described in [1,2,3].

# Description of the script
---------------------------

- **`analytical_pb.m`** :   Main script (wrapper).
- **`multi_spheres_cg5.m`** :   Computes the expansion coefficients of the potentials, used for further calculations, and the total energy.
- **`calc_energy_components0.m`** :   Performs energy partitioning as described in [4].
- **`calculate_potentials1.m`** :   Computes the potential at specified points.
- **`array_cg15.mat`** :   Pre-computed Clebsch-Gordan coefficients for calculations with n_max<=15

# Bibliography
--------------

[1] Siryk, S. V., Bendandi, A., Diaspro, A., & Rocchia, W. (2021). Charged dielectric spheres interacting in electrolytic solution: A linearized Poisson–Boltzmann equation model. The Journal of Chemical Physics, 155(11).

[2] Siryk, S. V., & Rocchia, W. (2022). Arbitrary-Shape Dielectric Particles Interacting in the Linearized Poisson–Boltzmann Framework: An Analytical Treatment. The Journal of Physical Chemistry B, 126(49), 10400-10426.

[3] Di Florio, V., Ansalone, P., Siryk, S. V., Decherchi, S., De Falco, C., & Rocchia, W. (2025). NextGenPB: An analytically-enabled super resolution tool for solving the Poisson-Boltzmann Equation featuring local (de) refinement. Computer Physics Communications, 109816.

[4] Rocchia, W., Alexov, E., & Honig, B. (2001). Extending the applicability of the nonlinear Poisson− Boltzmann equation: multiple dielectric constants and multivalent ions. The Journal of Physical Chemistry B, 105(28), 6507-6514.

