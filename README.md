# pde-variational-method

This repository contains codes for solving partial differential equations (PDEs) manually, by writing the linear system derived from the variational formula without relying on any external libraries.

## Contents

- **1.L_2_projection.py**: Implements the L²-projection method to approximate a function \( f(x) \) (specifically \( f(x) = 2x\sin(2\pi x) + 3 \)) over an interval. The code constructs the global mass matrix and load vector by hand, solves the resulting linear system, and compares the original function with its L²-projection. Useful for understanding how the L²-projection works and how to assemble the matrices element by element.

- **2.Tdistribution.py**: Solves a steady-state temperature distribution problem along a rod of length \( L \) using the finite element method. The code manually assembles the global stiffness matrix and load vector, taking into account variable coefficients and boundary conditions. It then solves the resulting linear system to obtain the temperature profile along the rod. This script is useful for understanding the assembly of FEM matrices for diffusion-type problems and the treatment of boundary conditions.

- **3.mesh_refinement.py**: Implements a method to solve a Poisson problem with Dirichlet boundary conditions using a finite element approach. It includes functions to construct the global stiffness matrix and load vector, perform mesh refinement based on an error indicator, and solve the linear system derived from the variational formulation. The script also compares the effect of different refinement parameters.

- **4.L_2_porj_in_2D.py**: Extends the L²-projection method to two dimensions. The script defines a simple triangular mesh and projects the function \( f(x, y) = x \cdot y \) onto the finite element space. It manually assembles the global mass matrix and load vector for the 2D case, solves the resulting linear system, and visualizes both the exact function and its L²-projection for comparison. This code is useful for understanding the extension of variational methods and matrix assembly to two-dimensional domains.

## Purpose

These codes were developed for personal practice to understand the variational method, the construction of the matrices that make up the linear system, and how techniques like L2 norm and mesh refinement work.

Feel free to explore the code and adapt it for your own learning and experimentation with PDEs and finite