%AMPLE 3D: A Material Point Learning Environment - 3D Electric Field Problem
%--------------------------------------------------------------------------
% Author: William Coombs (extended to 3D)
% Date:   27/02/2026
% Description:
% Steady-state 3D electric conduction material point method (MPM) code
% solving: -div(sigma * grad(phi)) = 0
% where phi = electric potential, sigma = electrical conductivity.
%
% Geometry: 3D rectangular domain with a spherical hole at the center.
% Elements: 8-noded tri-linear hexahedra.
% Boundary conditions:
%   Dirichlet: prescribed voltage (phi) at electrode faces
%     - Left face (x=0): phi = 0 (ground)
%     - Right face (x=lx): phi = 1 V (applied voltage)
%   Neumann: prescribed current source (fp) at material points
%
%--------------------------------------------------------------------------
% See also:
% SETUPGRID_ELECTRIC_3D   - 3D analysis specific information
% ELEMMPINFO_ELECTRIC     - material point-element information
% DETEXTFORCE_ELECTRIC    - external current sources
% DETFDOFS_ELECTRIC       - mesh unknown degrees of freedom
% LINSOLVE_ELECTRIC       - linear solver
% DETMPS_ELECTRIC         - material point conductance matrix and internal flux
% UPDATEMPS_ELECTRIC      - update material points
% POSTPRO_ELECTRIC_3D     - post processing function including VTK output
%--------------------------------------------------------------------------
clear; close all; clc;
addpath('constitutive','functions','plotting','setup');
[lstps,g,mpData,mesh] = setupGrid_electric_3d;                         % 3D setup information
NRitMax = 5; tol = 1e-9;                                              % Newton Raphson parameters
[nodes,nD] = size(mesh.coord);                                         % number of nodes and dimensions (3)
[nels,nen] = size(mesh.etpl);                                          % number of elements and nodes/element (8)
nDoF = nodes;                                                          % total number of degrees of freedom (scalar field)
nmp  = length(mpData);                                                 % number of material points
lstp = 0;                                                              % zero loadstep counter (for plotting function)

fprintf(1,'3D Electric Field MPM Simulation\n');
fprintf(1,'  Domain: %i x %i x %i\n', max(mesh.coord(:,1)), max(mesh.coord(:,2)), max(mesh.coord(:,3)));
fprintf(1,'  Nodes: %i, Elements: %i, Material Points: %i\n', nodes, nels, nmp);
fprintf(1,'  Dimensions: %i, Nodes/element: %i\n', nD, nen);

visualization_electric_3d(mpData, mesh);
run postPro_electric_3d;                                               % plotting initial state & mesh

phi_nodes = zeros(nDoF,1);                                             % initial electric potential (all zero)

for lstp=1:lstps                                                       % loadstep loop
  fprintf(1,'\n%s %4i %s %4i\n','loadstep ',lstp,' of ',lstps);
  [mesh,mpData] = elemMPinfo_Electric(mesh,mpData);                    % material point - element information
  fext = detExtForce_Electric(nodes,mpData);                           % external current source (Neumann BC)
  fext = fext*lstp/lstps;                                              % current external force value
  oobf = fext;                                                         % initial out-of-balance force
  fErr = 1;                                                            % initial error
  frct = zeros(nDoF,1);                                                % zero the reaction currents
  dphi_nodes = zeros(nDoF,1);                                          % zero the potential increment
  fd   = detFDoFs_Electric(mesh);                                      % free degrees of freedom
  NRit = 0;                                                            % zero the iteration counter
  Ke   = sparse(nDoF, nDoF);                                           % zero global conductance matrix
  while (fErr > tol) && (NRit < NRitMax) || (NRit < 2)                % global equilibrium loop
    [ddphi, drct] = linSolve_Electric(mesh.bc,Ke,oobf,NRit,fd);       % linear solver
    phi_nodes = phi_nodes + ddphi;                                     % update electric potential
    frct = frct+drct;                                                  % update reaction currents

    [fint,Ke,mpData] = detMPs_Electric(phi_nodes,mpData,nD);           % global conductance matrix & internal flux

    oobf = (fext-fint+frct);                                           % out-of-balance force
    fErr = norm(oobf)/norm(fext+frct+eps);                             % normalised oobf error
    NRit = NRit+1;                                                     % increment the NR counter
    fprintf(1,'%s %2i %s %8.3e\n','  iteration ',NRit,' NR error ',fErr);
  end
  mpData = updateMPs_Electric(phi_nodes,mpData);                       % update material point potentials
  run postPro_electric_3d;                                              % plotting and post processing
  visualization_electric_3d(mpData,mesh);
end
