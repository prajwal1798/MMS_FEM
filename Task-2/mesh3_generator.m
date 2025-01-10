clear all
close all
clc

% ELEMENTS Section
elements = 64; % Update the number of elements to 160

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

% Generate the first 8 rows based on the provided algorithm for columns 2, 3, and 4
% ===============================================================================
i= 1;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 6;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

 i= 1;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 6;
        elementData(idx, 4) = i + j + 5;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

% ============================================================================
 i= 6;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 6;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

i= 6;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 6;
        elementData(idx, 4) = i + j + 5;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

% ============================================================================
 i= 11;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 6;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

% ============================================================================
i= 11;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 6;
        elementData(idx, 4) = i + j + 5;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end
% =========================================================================
i= 16;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 6;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

i= 16;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 6;
        elementData(idx, 4) = i + j + 5;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

% ============================================================================

i= 21;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 6;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

i= 21;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 6;
        elementData(idx, 4) = i + j + 5;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

% ============================================================================

i= 26;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 6;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

i= 26;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 6;
        elementData(idx, 4) = i + j + 5;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

% ============================================================================

i= 31;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 6;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

i= 31;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 6;
        elementData(idx, 4) = i + j + 5;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

% ============================================================================

i= 36;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 1;
        elementData(idx, 4) = i + j + 6;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end

i= 36;
    for j = 0:3
        elementData(idx, 1) = idx;
        elementData(idx, 2) = i + j;
        elementData(idx, 3) = i + j + 6;
        elementData(idx, 4) = i + j + 5;
        elementData(idx, 5) = 1.5;
        idx = idx + 1;
    end
% % ============================================================================



% ============================================================================
% NODE_COORDINATES Section
% ============================================================================

nodes = 45; % Total number of nodes
nodeCoordinates = zeros(nodes, 4);
% Generate the node coordinates for the first 8 nodes
for i = 1:5
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 1) * 0.1;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.0;
    nodeCoordinates(i, 4) = 0.0;
end


for i = 6:10
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 6) * 0.1;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.075;  
    nodeCoordinates(i, 4) = 0.0;
end
    

for i = 11:15
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 11) * 0.1;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.150;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end

for i = 16:20
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 16) * 0.1;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.225;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end

for i = 21:25
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 21) * 0.1;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.300;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end

for i = 26:30
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 26) * 0.1;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.375;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end

for i = 31:35
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 31) * 0.1;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.450;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end

for i = 36:40
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 36) * 0.1;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.525;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end

for i = 41:45
    nodeCoordinates(i, 1) = i;  % First column: Node number i
    nodeCoordinates(i, 2) = (i - 41) * 0.1;  % Second column: i * 0.05
    nodeCoordinates(i, 3) = 0.600;
    nodeCoordinates(i, 4) = 0.0;% Third column: Fixed to 0.0
end


%===========================================================================
% Writing the generated data to a text file
%===========================================================================
fileID = fopen('Mesh3.txt', 'w');

fprintf(fileID, 'TITLE = TestFile\n\n');
fprintf(fileID, 'ELEMENTS = %d\n\n', elements);

% Writing the element data with consistent spacing and alignment for the first 8 rows
for i = 1:idx-1
    fprintf(fileID, '%-4d %-3d %-3d %-3d %.1f \n', elementData(i, :));
end

fprintf(fileID, '\n');


fprintf(fileID, 'NODE_COORDINATES =%d \n', nodes);
fprintf(fileID, '%-3d   %.4f   %.4f   %.4f \n', nodeCoordinates.');

fprintf(fileID, '\n');


fprintf(fileID, 'NODES_WITH_PRESCRIBED_TEMPERATURE = 10\n');
for i_u = 1:5
    fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
end

% for i_u = 6:5:36
%     fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
% end

for i_u = 41:45
    fprintf(fileID, '%2d  %7.3f \n', i_u, 200*(nodeCoordinates(i_u, 2)^2) + 200*nodeCoordinates(i_u, 3)^2 + 180);
end
fprintf(fileID, '\n');

fprintf(fileID, '\n');


fprintf(fileID, 'EDGES_WITH_PRESCRIBED_CONVECTION = 0\n\n');

fprintf(fileID, '\n');

fprintf(fileID, 'EDGES_WITH_PRESCRIBED_NON_ZERO_HEAT_FLUX = 8\n');
fprintf(fileID, '5    10    -240\n10   15    -240\n15   20    -240\n20   25    -240\n25   30    -240\n30   35    -240\n35   40    -240\n40   45    -240\');

fclose(fileID);

