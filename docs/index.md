---
title: iPIC3D
subtitle: Multi-scale simulations of plasma using Particle-In-Cell method
layout: page
hero_image: /iPIC3D/img/cover_pic.png
hero_height: is-fullheight
callouts: home_callouts
hero_link: /install/
hero_link_text: Download
show_sidebar: true
---

# Overview

iPIC3D is a three-dimensional parallel PIC code. On multiprocessor architectures, the domain decomposition technique is used to divide the computational workload among processors. For implicit PIC where the cost of particle moving and of field solving are of the same order (unlike explicit PIC where most of the cost resides with the particles), it is crucial that both field solving and particle moving be parallelized efficiently. An important aspect of efficiency is the need to retain the particles and cells belonging to a subdomain on the same processor. Large amounts of information is exchanged between grid and particles residing in the same physical domain and therefore it is important to avoid that this information exchange results in inter processor communication. The simulation box is divided among processors using a generic Cartesian virtual topology. Particles are divided among processors also depending on their location, and communicated to adjacent processors if exiting from the processor domain. The parallelization of the code is based on MPI libraries and blocking parallel communication has been chosen for the communication among processors.

# Applications


# Getting iPIC3D

iPIC3D source code and installation instruction is available [here](/install).

### System Requirement

- GNU Compilers or Intel Compilers
- A modern MPI installation
- Automake
- HDF5 (Optional)
