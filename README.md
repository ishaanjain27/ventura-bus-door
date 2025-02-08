# ventura-bus-door
CSE Final Minor Project, Group 1: Mathematical Modeling and Numerical Simulation of Driving Busses and Vibrating Doors

## Folder: point_mass_vibration
This folder explores the vibration analysis of a single-degree-of-freedom (SDOF) mass-damper system. The goal is to compare numerical and analytical approaches while experimenting with the Julia programming language to understand its importance in our study. Different time integration schemes are also analyzed.
Files:
* mass_damper_numerical.py – Solves the SDOF mass-damper system numerically using Python.
* mass_damper_numerical.jl – Implements the same numerical approach in Julia.
* mass_damper_analytical.py – Derives the analytical solution symbolically using SymPy.
* analytical_solution_single_mass.jl – Provides the analytical solution using Julia’s symbolic capabilities.
* Compare ODE and SecondOrderODE.ipynb – A Jupyter notebook in Julia that compares different time marching schemes, including first- and second-order methods, to assess accuracy and stability.

## Folder: two_mass_vibration
This folder focuses on the vibration analysis of a two-degree-of-freedom (2DOF) mass-damper system. Both numerical and analytical solutions are considered, with an emphasis on time integration and error analysis.
Files:
* 2mass_damped_forced.jl – Solves the 2DOF mass-damper system both numerically and analytically, investigating errors and system behavior.
* double_msd.jl – Implements a custom time-marching numerical solution, avoiding built-in Julia solvers to better understand numerical integration.

## Folder: n_mass_vibration
This folder extends the vibration analysis to an n-degree-of-freedom (nDOF) mass-damper system, generalizing the approach from previous cases. The focus is on solving the system numerically using both Julia's built-in solvers and a custom time integration implementation.
Files:
* ndimensional_msd.ipynb – Implements a custom time-marching scheme for solving the nDOF system, providing deeper insights into numerical integration.
* nmasses.ipynb – Uses Julia’s built-in solvers to solve the same problem efficiently, allowing for performance comparisons.

## Folder: euler_beam_equation
This folder contains one of the core components of our study: the numerical solution of the Euler-Bernoulli beam equation for both static and dynamic cases. The models allow for different boundary conditions (spring, free, clamped, and simply supported) and assume uniform loading, with flexibility to modify it for point loading.
Main Files:
* Euler_Bernoulli_static.jl – Computes the static deformation of a beam, including eigenshape analysis for different boundary conditions.
* Euler_Bernoulli_dynamic.jl – Solves the dynamic beam equation using time integration methods.
Subfolder: dynamic-euler-beam
This subfolder contains initial implementation attempts in Python, exploring different discretization and time-marching approaches.
* eulerbeambcpaper.pdf – A reference paper that guided the discretization process and the application of boundary conditions.

## Folder: kirchhoff_love_plate_equation
This folder contains the final and most significant results of our study, focusing on the numerical solution of the Kirchhoff-Love plate equation. The models analyze both static and dynamic behavior for different boundary conditions. All implementations are in Julia.
Main Files:
* plate_static.jl – Solves the static plate equation for a simply supported plate.
* plate_dynamic.jl – Implements time integration to solve the dynamic plate equation for a simply supportedplate.
* plate_static_spring.jl – Extends the static analysis to a cantilever plate with spring supports on one edge.
* plate_dynamic_spring.jl – Computes the dynamic response of a cantilever plate with spring boundary conditions.