clear all
close all
clc

% ELEMENTS Section
elements = 16; % Update the number of elements to 40

elementData = zeros(elements, 5);
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

% Generate the first 4 rows based on the provided algorithm for columns 2, 3, and 4
% ===============================================================================
i= 1;
    for j = 0:1
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 4;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

 i= 1;
    for j = 0:1
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 4;
        elementData(idx, 4) = i + j + 3;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

% ============================================================================
 i= 4;
    for j = 0:1
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 4;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

i= 4;
    for j = 0:1
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 4;
        elementData(idx, 4) = i + j + 3;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

% ============================================================================
 i= 7;
    for j = 0:1
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 4;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end


i= 7;
    for j = 0:1
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 4;
        elementData(idx, 4) = i + j + 3;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end
% =========================================================================

i= 10;
    for j = 0:1
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 4;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

i= 10;
    for j = 0:1
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 4;
        elementData(idx, 4) = i + j + 3;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end



% ============================================================================
% NODE_COORDINATES Section
% ============================================================================
nodes = 15; % Total number of nodes
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
    nodeCoordinates(i, 3) = 0.15;  
    nodeCoordinates(i, 4) = 0.0;
end
    

for i = 7:9
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 7) * 0.2;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.30;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end

for i = 10:12
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 10) * 0.2;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.45;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end

for i = 13:15
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 13) * 0.2;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.60;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end


%===========================================================================
% Writing the generated data to a text file
%===========================================================================
fileID = fopen('Mesh2.txt', 'w');

fprintf(fileID, 'TITLE = TestFile\n\n');
fprintf(fileID, 'ELEMENTS = %d\n\n', elements);

% Writing the element data with consistent spacing and alignment for the first 8 rows
for i = 1:idx-1
    fprintf(fileID, '%-4d %-3d %-3d %-3d %.1f \n', elementData(i, :));
end

fprintf(fileID, '\n');


fprintf(fileID, 'NODE_COORDINATES = %d \n', nodes);
fprintf(fileID, '%-3d   %.4f   %.4f   %.4f \n', nodeCoordinates.');

fprintf(fileID, '\n');


fprintf(fileID, 'NODES_WITH_PRESCRIBED_TEMPERATURE = 6\n');
for i_u = 1:3
    fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
end

% for i_u = 4:3:10
%     fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
% end

for i_u = 13:15
    fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
end
fprintf(fileID, '\n');


fprintf(fileID, 'EDGES_WITH_PRESCRIBED_CONVECTION = 0\n\n');

fprintf(fileID, '\n');

fprintf(fileID, 'EDGES_WITH_PRESCRIBED_NON_ZERO_HEAT_FLUX = 4\n');
fprintf(fileID, '3     6   -240\n6     9   -240\n9    12   -240\n12   15   -240\');

fclose(fileID);

