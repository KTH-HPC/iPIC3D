---
title: Installing iPIC3D
subtitle: Download and compilation instruction
layout: page
show_sidebar: false
menubar: documentation_menu
hide_hero: true
hide_footer: true
---

## Install iPIC3D

Download and install iPIC3D using the following commands.

### Download

Clone the repository from Github:

#### Linux
```bash
git clone https://github.com/KTH-HPC/iPIC3D.git

```

### Build and install

Configure and build the software:
```bash
cd ipic3d-klm
mkdir build
cd build
CC=mpicc CXX=mpicxx cmake ..
make
```
Optionally install the software:
```bash
make install
```

[comment]: # (#### CUDA)
[comment]: # (Particle mover with CUDA can be compiled if enabled)
[comment]: # (```bash)
[comment]: # (--with-cuda=${CUDA_PATH})
[comment]: # (```)

## Detail user guide
A detailed user guide can be found at **[here](https://www.kth.se/files/view/markidis/5846cba0a4ba671c5d0bfb2d/ee_pdc_course_project_report_v.pdf)**.
