Cell Migration Simulator
========================

This tool allows to simulate the chemotactic migration of cells 
degrading the chemoattractant with a very simple model. Self-driven chemotactic behavior has been 
extensively studied, both experimentally and numerically, by Tweedy et al., who 
showed how this mechanism allows cells to find their way in complex environments such as mazes [Tweedy2020]_.
The program supports systems consisting of a 2D discrete environment, where the cells jump
from node to node until they reach a chemostat. Some nodes can be defined as "walls" to 
chemoattractant diffusion and cell migration, which allows for the definition of various system landscapes.
Reaction and diffusion are handled with the tau-leap algorithm [Gillespie2001]_ using
Bernstein's approach to model diffusion [Bernstein2005]_, as it is done with STReNGTHS [STReNGTHS2024]_.

How to use it 
-------------

First, one has to compile "interface.cpp" as a shared library with "sl.dll" as name.
Using the -fopenmp option to benefit from OpenMP acceleration is advised.

Next, one can import the sdchemotaxis module (sdchemotaxis.py) and run simulations with the "simulate" functions.
The example.py file has a few examples of how to run simulations and visualize the results. 
Functions and classes of the sdchemotaxis module are documented through docstrings.

References
----------

.. [Tweedy2020] Tweedy, L., Thomason, P. A., Paschke, P. I., Martin, K., Machesky, L. M., Zagnoni, M., & Insall, R. H. (2020). Seeing around corners: cells solve mazes and respond at a distance using attractant breakdown. Science, 369(6507), eaay9792. doi:10.1126/science.aay9792
.. [Gillespie2001] Gillespie, D. T. (2001). Approximate accelerated stochastic simulation of chemically reacting systems. The Journal of Chemical Physics, 115(4), 1716-1733. doi:10.1063/1.1378322
.. [Bernstein2005] Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm. Physical review. E, Statistical, nonlinear, and soft matter physics, 71(4), 041103. doi:10.1103/PhysRevE.71.041103
.. [STReNGTHS2024] Fillion, T. & Piazza, F. (2024). STReNGTHS, a Python package to model and simulate complex reaction-diffusion systems. The Journal of Open Source Software, 9(97), 6495. doi::10.21105/joss.06495

