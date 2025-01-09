%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for the finite element analysis of 2D heat flow problems
% with 3-noded triangular elements
%
% Dr. Rubén Sevilla
%
% Zienkiewicz Centre for Computational Engineering 
% College of Engineering
% Swansea University
%
% MODULE: EG-M23 Finite Element Computational Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear screen, workspace and close open figures
clc; clearvars; close all;
%==========================================================================
% Set some basic variables
%==========================================================================
nOfElemNodes=3;
%==========================================================================
% Data input phase. Read information from data file
%==========================================================================
fileName = input('Input data file name -------> ','s');
fileID = fopen(fileName, 'r');
fileTitle = fscanf(fileID, 'TITLE = %s',1);
% Total number of elements in the mesh
nOfElements = fscanf(fileID, '\nELEMENTS = %d', 1);

%--------------------------------------------------------------------------
% Read element information
%--------------------------------------------------------------------------
% Table of connectivities (also read conductivities and heat source data)

connect = fscanf(fileID, '\n%d %d %d %d %f', [nOfElemNodes+2,nOfElements]);
connect = connect';

% Extract element conductivity (isotropic conductivity assumed)
conductivity = connect(:,nOfElemNodes+2);
% Extract element heat source data
% Extract connectivity table
elementOrder = connect(:,1);
connect = connect(:,2:nOfElemNodes+1);
connect(elementOrder,:) = connect;
%--------------------------------------------------------------------------
% Read nodal coordinates
%--------------------------------------------------------------------------
nOfNodes = fscanf(fileID, '\nNODE_COORDINATES = %d', 1);
coord = fscanf(fileID, '\n%d %f %f %f', [4, nOfNodes]);
coord = coord';
% Extract nodal coordinates
nodeOrder = coord(:,1);
hs = coord(:,4);
coord = coord(:,2:3);
% First, second and third vertexes of all elements
x1 = coord(connect(:,1),1);
x2 = coord(connect(:,2),1);
x3 = coord(connect(:,3),1);
y1 = coord(connect(:,1),2);
y2 = coord(connect(:,2),2);
y3 = coord(connect(:,3),2);
% Calculate element areas
elemArea = zeros(nOfElements,1);
for iElem=1:nOfElements
    elemArea(iElem) = ( (x1(iElem)-x3(iElem))*(y2(iElem)-y3(iElem))-...
        (y1(iElem)-y3(iElem))*(x2(iElem)-x3(iElem)))/2;
    if  elemArea(iElem) < 0
        error('element nodes numbered clock-wise in data file')
    end
end
%--------------------------------------------------------------------------
% Read prescribed temperature
%--------------------------------------------------------------------------
nOfNodesDirichlet = fscanf(fileID,'\nNODES_WITH_PRESCRIBED_TEMPERATURE = %d',1);
dirichletData = fscanf(fileID, '\n%d %f', [2, nOfNodesDirichlet]);
dirichletData = dirichletData';
% Set fixed dof's information arrays
nodesDirichlet = dirichletData(:,1);
valueDirichlet = dirichletData(:,2);
% Create a vector of length nOfNodes and set 1 if the node has a imposed
% temperature or 0 otherwise
isNodeFixed = zeros(nOfNodes,1);
isNodeFixed(nodesDirichlet) = ones(nOfNodesDirichlet,1);
% Store a list of nodes with no imposed temperature (this is used in the
% elimintation process)
nodesFree = setdiff(1:nOfNodes, nodesDirichlet);
nOfNodesFree = length(nodesFree);
%--------------------------------------------------------------------------
% Read convection data
%--------------------------------------------------------------------------
nOfEdgesConvection = fscanf(fileID,'\nEDGES_WITH_PRESCRIBED_CONVECTION = %d',1);
%--------------------------------------------------------------------------
% Read Neumann Boundary data
%--------------------------------------------------------------------------
nOfEdgesNeumann = fscanf(fileID,'\nEDGES_WITH_PRESCRIBED_NON_ZERO_HEAT_FLUX = %d',1);
NeumannTable = fscanf(fileID,'\n%d %d %f', [3, nOfEdgesNeumann]);
NeumannTable = NeumannTable';
% Extract the connectivities and values of heat flux
NeumannNodes = NeumannTable(:,1:2);
q_edge = NeumannTable(:,3);
% First and second vertexes of all edges
xa = coord(NeumannNodes(:,1),1);
xb = coord(NeumannNodes(:,2),1);
ya = coord(NeumannNodes(:,1),2);
yb = coord(NeumannNodes(:,2),2);
% Compute the length of the edges on a Nuemann boundary
edgeLength = zeros(1, nOfEdgesNeumann);
for iEdge = 1:nOfEdgesNeumann
    edgeLength(iEdge) = sqrt( (xb(iEdge)-xa(iEdge))^2 + (yb(iEdge)-ya(iEdge))^2 );
end
%==========================================================================
% Solution phase. Assemble stiffness, forcing terms and solve system
%==========================================================================
%--------------------------------------------------------------------------
% Assemble global heat rate vector (external heat sources contribution)
%--------------------------------------------------------------------------
cr = connect';
RQ = zeros(nOfNodes,1);
for iElem = 1:nOfElements
    % Compute element heat source
    NTN = [1/12 1/24 1/24; 1/24 1/12 1/24; 1/24 1/24 1/12];
    Qe = [hs(cr(1,(iElem)),1); hs(cr(2,(iElem)),1); hs(cr(3,(iElem)),1)];
    rQ = 2*0.0060*NTN*Qe;

    % Assembly: Add element contribution to global vector
    globalElementNodes = connect(iElem,:);
    RQ(globalElementNodes) = RQ(globalElementNodes) + rQ;
end
%--------------------------------------------------------------------------
% Compute global stiffness matrix
%--------------------------------------------------------------------------
K = zeros(nOfNodes,nOfNodes);
for iElem = 1:nOfElements
    % Compute element conductivity stiffness
    Bmatx = 1/(2*elemArea(iElem))*...
            [y2(iElem)-y3(iElem)  y3(iElem)-y1(iElem)  y1(iElem)-y2(iElem);
             x3(iElem)-x2(iElem)  x1(iElem)-x3(iElem)  x2(iElem)-x1(iElem)];
    ke = conductivity(iElem)*elemArea(iElem)*(Bmatx'*Bmatx);
    
    % Assembly: Add element contribution to global stiffness
    globalElementNodes = connect(iElem,:);
    K(globalElementNodes,globalElementNodes) = K(globalElementNodes,globalElementNodes) + ke;
end
%--------------------------------------------------------------------------
% Add Neumann contributions to Forcing Vector 
%--------------------------------------------------------------------------

Rq = zeros(nOfNodes,1);
Rq(5)  = (edgeLength(1) * q_edge(1))/2 * 1;
Rq(10) = (edgeLength(1) * q_edge(1))/2 * 1 + (edgeLength(2) * q_edge(2))/2 * 1;
Rq(15) = (edgeLength(2) * q_edge(2))/2 * 1 + (edgeLength(3) * q_edge(3))/2 * 1;
Rq(20) = (edgeLength(3) * q_edge(3))/2 * 1 + (edgeLength(4) * q_edge(4))/2 * 1;
Rq(25) = (edgeLength(4) * q_edge(4))/2 * 1 + (edgeLength(5) * q_edge(5))/2 * 1;
Rq(30) = (edgeLength(5) * q_edge(5))/2 * 1;

%--------------------------------------------------------------------------
% Global right hand side of the system
%--------------------------------------------------------------------------
F = RQ - Rq;
%--------------------------------------------------------------------------
% Apply the elimination process to account for Dirichlet boundaries
%--------------------------------------------------------------------------
% Retain only equations for nodes that are NOT in a Dirichlet boundary
Kstar = K(nodesFree,nodesFree);
Fstar = F(nodesFree);
% Modify the right hand side with prescribed values
for iRow = 1:nOfNodesFree
    nodeRow = nodesFree(iRow);
    Krow  = K(nodeRow,nodesDirichlet);
    Fstar(iRow) = Fstar(iRow) - Krow*valueDirichlet;
end
%--------------------------------------------------------------------------
% Solve the reduced system
%--------------------------------------------------------------------------
Tstar = Kstar\Fstar;
%--------------------------------------------------------------------------
% Add to the global vector the solution at Dirichlet nodes (known values)
%--------------------------------------------------------------------------
ngdof = nOfNodes;
T = zeros(ngdof,1);
T(nodesFree)= Tstar;
T(nodesDirichlet) = valueDirichlet;
%==========================================================================
% Post-processing phase. 
%==========================================================================
% Compute reactions (heat flux) at nodes with prescribed temperature
R = zeros(ngdof,1);
for iRow = 1:nOfNodesDirichlet
    nodeRow = nodesDirichlet(iRow);
    Krow  = K(nodeRow,:);
    R(nodeRow) = Krow*T-F(nodeRow);
end
% Compute fluxes
% Plot mesh with heat flow vector distribution
xcentre = (x1 + x2 + x3)/3; 
ycentre = (y1 + y2 + y3)/3;
q = zeros(nOfElements, 2);
for iElem = 1:nOfElements
    % Compute element heat flow vector
    Bmatx = 1/(2*elemArea(iElem))*...
            [y2(iElem)-y3(iElem)  y3(iElem)-y1(iElem)  y1(iElem)-y2(iElem);
             x3(iElem)-x2(iElem)  x1(iElem)-x3(iElem)  x2(iElem)-x1(iElem)];
    q(iElem,:) = -conductivity(iElem)*Bmatx*T(connect(iElem,:));
end
% Plot mesh with temperature distribution
figure(1)
trisurf(connect,coord(:,1),coord(:,2),T,'facecolor','interp','facelight','phong','edgecolor','k');
set(gca,'FontSize',22)
xlabel('x')
ylabel('y')
axis equal
view(2)
colormap jet
colorbar
% Plot heat flow vectors with background mesh
figure(2)
quiver(xcentre,ycentre,q(:,1),q(:,2));
hold on
triplot(connect,coord(:,1),coord(:,2),'k');
set(gca,'FontSize',22)
xlabel('x')
ylabel('y')
axis equal
view(2)
colormap jet

%==========================================================================
% Close file(s) before terminating program
%==========================================================================
status = fclose(fileID);