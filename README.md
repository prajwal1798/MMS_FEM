# MMS_FEM
Validating FEM for Elliptic Problems using Method of Manufactured Solutions

# **Finite Element Computational Analysis (MATLAB Project)**  

## **Project Overview**  
This repository contains MATLAB code developed as part of the **Finite Element Computational Analysis (EGM-23)** coursework at Swansea University. The project focuses on solving a **2D steady-state heat conduction problem** using the **Finite Element Method (FEM)**, with extensions to incorporate higher-order elements, internal heat generation, and validation techniques.  

## **Key Features**  
- **Basic Extensions**:  
  - Implementation of linearly varying internal heat sources within elements.  
  - Handling of non-zero boundary normal heat flux for Neumann boundary conditions.  
  - Validation of FEM results using the **Method of Manufactured Solutions (MMS)**.  

- **Mesh Generation and Accuracy Analysis**:  
  - Generation of FEM meshes with varying element sizes.  
  - L2-Norm error analysis to verify error propagation behavior.  

- **Higher-Order Finite Elements**:  
  - Support for higher-order polynomial basis functions in 4-node elements.  
  - Weak formulations for conductivity and heat source terms implemented with **Gaussian Quadrature**.  

- **Critical Results**:  
  - Verified temperature and flux distributions in line with prescribed boundary conditions.  
  - Achieved expected error convergence rates, validating the solver's stability and accuracy.  

## **Technical Details**  
- **Weak Formulations**:  
  - The FEM weak formulations for heat conduction include contributions from both Dirichlet (essential) and Neumann (natural) boundary conditions.  
  - Gaussian Quadrature was employed for numerical integration in higher-order element formulations.  

- **Error Analysis**:  
  - L2-Norm error convergence was calculated based on manufactured solutions, confirming the relationship between mesh refinement and error decay.  

## **How to Use**  
1. Clone the repository:  
   ```bash
   git clone https://github.com/prajwal_1798/finite-element-analysis.git
   

