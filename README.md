# Advection CMC42 Yin–Yang RK4 Solver

## Overview
This project provides a numerical solver for the 2D advection equation on the Yin–Yang grid using a compact fourth-order MacCormack scheme with a 4/2 spatial discretization formulation and fourth-order Runge–Kutta (RK4) time integration.

The implementation is designed for high-order accuracy in spherical geometry without coordinate singularities, achieved through the Yin–Yang grid approach.

---

## Numerical Method

- **Governing equation:** 2D linear advection equation  
- **Spatial discretization:** Compact fourth-order MacCormack scheme (4/2 formulation)  
- **Time integration:** Classical fourth-order Runge–Kutta (RK4)  
- **Grid system:** Yin–Yang overset grid for spherical domains  

---

## Features

- High-order accurate finite difference scheme  
- Singularity-free spherical grid (Yin–Yang grid)  
- Efficient time integration using RK4  
- Suitable for geophysical fluid dynamics applications  

---

## Output and Visualization

Simulation results are exported in a format compatible with **Tecplot** for post-processing and visualization.

---

## File Description

- `Adv_CMC42_RYY.f90` → Main Fortran source code implementing the solver

---

## Requirements

- Fortran compiler (e.g., gfortran or Intel Fortran)
- Tecplot (for visualization)

---

## How to Compile

Example using gfortran:

```bash
gfortran Adv_CMC42_RYY.f90 -o solver.exe
