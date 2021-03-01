---
title: Run simulation
subtitle: Running example simulation with iPIC3D
layout: page
show_sidebar: false
menubar: documentation_menu
hide_hero: true
hide_footer: true
---

## Run simulation

Assume input file ```inputfiles/testGEM2Dsmall.inp``` is used and the number of MPI subdomains in each direction is set to (2,2,1):
```bash
mpirun -n 4 ./build/iPIC3D inputfiles/testGEM2Dsmall.inp
```
Ensure that the output directory exists as it will not be created automatically.

A detailed user guide can be found at **[here](https://www.kth.se/files/view/markidis/5846cba0a4ba671c5d0bfb2d/ee_pdc_course_project_report_v.pdf)**.
