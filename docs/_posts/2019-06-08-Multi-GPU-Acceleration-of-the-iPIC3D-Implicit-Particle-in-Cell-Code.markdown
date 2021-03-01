---
layout: post
title:  "Multi-GPU Acceleration of the iPIC3D Implicit Particle-in-Cell Code"
date:   2019-06-08 00:00:00
categories: publications
description: "Multi-GPU Acceleration of the iPIC3D Implicit Particle-in-Cell Code"
published: true
show_sidebar: false
---

iPIC3D is a widely used massively parallel Particle-in-Cell code for the simulation of space plasmas. However, its current implementation does not support execution on multiple GPUs. In this paper, we describe the porting of iPIC3D particle mover to GPUs and the optimization steps to increase the performance and parallel scaling on multiple GPUs. We analyze the strong scaling of the mover on two GPU clusters and evaluate its performance and acceleration. The optimized GPU version which uses pinned memory and asynchronous data prefetching outperform their corresponding CPU versions by  5−10× on two different systems equipped with NVIDIA K80 and V100 GPUs.

[https://doi.org/10.1007/978-3-030-22750-0_58](https://doi.org/10.1007/978-3-030-22750-0_58){:target="_blank"}

##### Authors
C. P. Sishtla, S. W. D. Chien, V. Olshevsky, E. Laure and S. Markidis

##### Venue
International Coference on Computational Science 2019

##### Cite
```bibtex
@InProceedings{10.1007/978-3-030-22750-0_58,
               author="Sishtla, Chaitanya Prasad
               and Chien, Steven W. D.
               and Olshevsky, Vyacheslav
               and Laure, Erwin
               and Markidis, Stefano",
               title="Multi-GPU Acceleration of the iPIC3D Implicit Particle-in-Cell Code",
               booktitle="Computational Science -- ICCS 2019",
               year="2019",
               publisher="Springer International Publishing",
               address="Cham",
               pages="612--618",
               abstract="iPIC3D is a widely used massively parallel Particle-in-Cell code for the simulation of space plasmas. However, its current implementation does not support execution on multiple GPUs. In this paper, we describe the porting of iPIC3D particle mover to GPUs and the optimization steps to increase the performance and parallel scaling on multiple GPUs. We analyze the strong scaling of the mover on two GPU clusters and evaluate its performance and acceleration. The optimized GPU version which uses pinned memory and asynchronous data prefetching outperform their corresponding CPU versions by {\$}{\$}5-10{\backslash}times {\$}{\$}on two different systems equipped with NVIDIA K80 and V100 GPUs.",
               isbn="978-3-030-22750-0"
}
```
