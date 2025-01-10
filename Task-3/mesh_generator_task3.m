%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODULE: EG-M23 Finite Element Computational Analysis
% Program for TASK 2 of Coursework by Group #3 
% MESH-1 Generator (h = 0.4)
%
% Prajwal Bharadwaj - 2337862
%
% Zienkiewicz Centre for Computational Engineering 
% College of Engineering
% Swansea University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

% ELEMENTS Section
elements = 2; % Update the number of elements to 40

elementData = zeros(elements, 8);
idx = 1;

%==========================================================================
% Using Manufactured Data for getting the files written for Dirichlet Data 
%==========================================================================

x = linspace(0, 0.4, 100);
y = linspace(0, 0.6, 100);

[X, Y] = meshgrid(x, y);

T_manufactured = 200*X.^2 + 200*Y.^2 + 180;

%======================================================================
% Extracting the Dirichlet Data for the problem under consideration  
%======================================================================

T_Dirichlet_1 = T_manufactured(1,:)';
T_Dirichlet_2 = T_manufactured(end,:)';

% ===============================================================================
% Generate the first 4 rows based on the provided algorithm for columns 2, 3, and 4
% ===============================================================================
i= 1;
    for j = 0
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 2;
        elementData(idx, 4) = i + j + 8;
        elementData(idx, 5) = i + j + 1;
        elementData(idx, 6) = i + j + 5;
        elementData(idx, 7) = i + j + 4;
        elementData(idx, 8) = 1.5;
        idx = idx + 1;
    end

 i= 1;
    for j = 0
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 8;
        elementData(idx, 4) = i + j + 6;
        elementData(idx, 5) = i + j + 4;
        elementData(idx, 6) = i + j + 7;
        elementData(idx, 7) = i + j + 3;
        elementData(idx, 8) = 1.5;
        idx = idx + 1;
    end

% ============================================================================
%  i= 3;
%     for j = 0
%         elementData(idx, 1) = idx;
%         elementData(idx, 2) = i + j;
%         elementData(idx, 3) = i + j + 1;
%         elementData(idx, 4) = i + j + 3;
%         elementData(idx, 5) = 1.5;
%         idx = idx + 1;
%     end
% 
% i= 3;
%     for j = 0
%         elementData(idx, 1) = idx;
%         elementData(idx, 2) = i + j;
%         elementData(idx, 3) = i + j + 3;
%         elementData(idx, 4) = i + j + 2;
%         elementData(idx, 5) = 1.5;
%         idx = idx + 1;
%     end

% ============================================================================



% ============================================================================
% NODE_COORDINATES Section
% ============================================================================
nodes = 9; % Total number of nodes
nodeCoordinates = zeros(nodes, 4);

% Generate the node coordinates for the first 8 nodes
for i = 1:3
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 1) * 0.2;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.0;  
    nodeCoordinates(i, 4) = 0.0;
end


for i = 4:6
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 4) * 0.2;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.30;  
    nodeCoordinates(i, 4) = 0.0;
end
    

for i = 7:9
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 7) * 0.2;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.60;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end


%===========================================================================
% Writing the generated data to a text file
%===========================================================================
fileID = fopen('QuadMesh_task3.txt', 'w');

fprintf(fileID, 'TITLE = TestFile\n\n');
fprintf(fileID, 'ELEMENTS = %d\n\n', elements);

% Writing the element data with consistent spacing and alignment for the first 8 rows
for i = 1:idx-1
    fprintf(fileID, '%-4d %-3d %-3d %-3d %-3d %-3d %-3d %.1f \n', elementData(i, :));
end

fprintf(fileID, '\n');


fprintf(fileID, 'NODE_COORDINATES = %d \n', nodes);
fprintf(fileID, '%-3d   %.4f   %.4f   %.4f \n', nodeCoordinates.');

fprintf(fileID, '\n');


fprintf(fileID, 'NODES_WITH_PRESCRIBED_TEMPERATURE = 8\n');
for i_u = 1:3
    fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
end

for i_u = 4
    fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
end

for i_u = 6
    fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
end

for i_u = 7:9 
    fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
end
fprintf(fileID, '\n');


fprintf(fileID, 'EDGES_WITH_PRESCRIBED_CONVECTION = 0\n\n');

fprintf(fileID, '\n');

fprintf(fileID, 'EDGES_WITH_PRESCRIBED_NON_ZERO_HEAT_FLUX = 2\n');
fprintf(fileID, '3   6   -240\n6   9   -240\');

fclose(fileID);

