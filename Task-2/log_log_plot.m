%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODULE: EG-M23 Finite Element Computational Analysis
% Program for TASK 2 of Coursework by Group #3 
% with 3-noded triangular elements
%
% Prajwal Bharadwaj - 2337862
%
% Zienkiewicz Centre for Computational Engineering 
% College of Engineering
% Swansea University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear screen, workspace and close open figures
clc; clearvars; close all;


% Given error values for corresponding mesh sizes
error_values = [2.8422e-14, 1.5567e-13, 6.927e-13];
mesh_sizes = [0.1, 0.2, 0.4];

% Take the logarithm of error values and mesh sizes
log_error_values = log10(error_values);
log_mesh_sizes = log10(mesh_sizes);

% Plotting the log-log plot with a line
plot(log_mesh_sizes, log_error_values, 'r-*', 'LineWidth', 2);
title('Log-Log Plot of Error Values vs Mesh Sizes');
xlabel('Log(h)');
ylabel('Log(||e||)');
grid minor
grid on;


