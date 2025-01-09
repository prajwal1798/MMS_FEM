# MMS_FEM
Validating FEM for Elliptic Problems using Method of Manufactured Solutions
Finite Element Computational Analysis (MATLAB Project)
Project Overview
This project focuses on solving a 2D heat conduction problem using the Finite Element Method (FEM) implemented in MATLAB. The work extends existing FEM code by incorporating higher-order elements, internal heat generation, and verifying accuracy through analytical solutions. The project is part of the coursework for the Finite Element Computational Analysis (EGM-23) module at Swansea University.

Key Features
Basic Code Extensions:

Implementation of internal heat generation distributed linearly across elements.
Support for non-zero boundary normal heat flux (Neumann boundary conditions).
Validation of the extended code using the Method of Manufactured Solutions (MMS).
Mesh Generation and Accuracy Analysis:

Generated finite element meshes with varying element sizes and performed L2-Norm error analysis.
Verified error propagation behavior with respect to element sizes and polynomial orders, achieving the expected gradient of ~2.
Higher-Order Finite Element Formulation:

Extended the FEM code to support higher-order polynomial basis functions for 4-node elements.
Weak formulations for conductivity and heat source terms were derived and implemented using Gaussian Quadrature.
Critical Observations:

Demonstrated temperature and heat flux distributions in line with prescribed boundary conditions.
Showed reduced error convergence rates with finer meshes, validating code stability and accuracy.
Key Techniques
Weak formulation of FEM governing equations for both Dirichlet and Neumann boundary conditions.
Use of isoparametric formulation for heat source integration over elements.
Application of Gaussian Quadrature for higher-order basis functions.
Method of Manufactured Solutions for code verification and validation.
Results
Successfully verified code accuracy with L2-Norm error convergence.
Extended the solver to handle higher-order elements and linearly varying heat sources.
Visualized temperature contours and heat flux distributions for various test cases, validating against analytical solutions.
Technologies Used
Programming Language: MATLAB
Key Concepts: Finite Element Method (FEM), Weak Formulations, Gaussian Quadrature, Error Analysis
Motivation
This project was undertaken as part of my MSc in Computational Engineering at Swansea University. It highlights the integration of theoretical FEM concepts with practical coding implementations to solve real-world engineering problems.
