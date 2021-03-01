---
layout: post
title:  "High Performance Solvers for Implicit Particle in Cell Simulation"
date:   2013-06-01 00:00:00
categories: publications
description: "High Performance Solvers for Implicit Particle in Cell Simulation"
published: true
show_sidebar: false
---

A three-dimensional implicit particle-in-cell (iPIC3D) method implemented by S. Markidis et. al. in [“Multiscale simulations of plasma with iPIC3D”, Mathematics and Computers in Simulation, 80(2010), 1509-1519] allows time steps at magnetohy- drodynamics time scale. The code requires the solution of two linear systems: a Poisson system related to divergence cleaning, and a system related to a second order formulation of Maxwell equation. In iPIC3D, the former is the most costly.

To reduce the cost of solving the Poisson system, a parallel matrix assembly and partitioning method are implemented, and conjugate gradient and algebraic multigrid (AMG) solvers from the Hypre library are called. The scalability of AMG as a solver is studied for 1D and 3D partitionings and compared to that of CG.

[https://doi.org/10.1016/j.procs.2013.05.396](https://doi.org/10.1016/j.procs.2013.05.396){:target="_blank"}

##### Authors
P. Kumar, S. Markidis, G. Lapenta, K. Meerbergen and D. Roose

##### Venue
_Procedia Computer Science_, Volume 18, 2013, Pages 2251-2258

##### Cite
```bibtex
@article{KUMAR20132251,
         title = "High Performance Solvers for Implicit Particle in Cell Simulation",
         journal = "Procedia Computer Science",
         volume = "18",
         pages = "2251 - 2258",
         year = "2013",
         note = "2013 International Conference on Computational Science",
         issn = "1877-0509",
         doi = "https://doi.org/10.1016/j.procs.2013.05.396",
         url = "http://www.sciencedirect.com/science/article/pii/S1877050913005395",
         author = "Pawan Kumar and Stefano Markidis and Giovanni Lapenta and Karl Meerbergen and Dirk Roose",
}
```
