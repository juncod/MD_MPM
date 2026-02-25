%AMPLE 1.1: A Material Point Learning Environment - Electric Field Problem
%--------------------------------------------------------------------------
% Author: William Coombs (adapted for electric field problem)
% Date:   27/08/2020
% Description:
% Steady-state electric conduction material point method (MPM) code
% solving: -div(sigma * grad(phi)) = 0
% where phi = electric potential, sigma = electrical conductivity.
% Boundary conditions:
%   Dirichlet: prescribed voltage (phi = V) at electrode nodes
%   Neumann:   prescribed current source (fp) at material points
%
%--------------------------------------------------------------------------
% See also:
% SETUPGRID_ELECTRIC    - analysis specific information
% ELEMMPINFO_ELECTRIC   - material point-element information
% DETEXTFORCE_ELECTRIC  - external current sources
% DETFDOFS_ELECTRIC     - mesh unknown degrees of freedom
% LINSOLVE_ELECTRIC     - linear solver
% DETMPS_ELECTRIC       - material point conductance matrix and internal flux
% UPDATEMPS_ELECTRIC    - update material points
% POSTPRO               - post processing function including vtk output
%--------------------------------------------------------------------------
clear; close all; clc;
addpath('constitutive','functions','plotting','setup');
[lstps,g,mpData,mesh] = setupGrid_electric;                            % setup information
NRitMax = 5; tol = 1e-9;                                              % Newton Raphson parameters
[nodes,nD] = size(mesh.coord);                                         % number of nodes and dimensions
[nels,nen] = size(mesh.etpl);                                          % number of elements and nodes/element
nDoF = nodes;                                                          % total number of degrees of freedom (scalar field)
nmp  = length(mpData);                                                 % number of material points
lstp = 0;                                                              % zero loadstep counter (for plotting function)

visualization_electric(mpData, mesh);
run postPro_electric;                                                  % plotting initial state & mesh

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
  run postPro_electric;                                                 % plotting and post processing
  visualization_electric(mpData,mesh);
end
