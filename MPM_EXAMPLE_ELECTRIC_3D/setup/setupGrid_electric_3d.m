function [lstps,g,mpData,mesh] = setupGrid_electric_3d
%% Analysis parameters - 3D Electric Conduction Problem
% Extends the 2D electric field problem to 3D with a spherical hole
% PDE: -div(sigma * grad(phi)) = 0
% Boundary conditions:
%   Dirichlet: phi=0 at x=0 (ground), phi=1 at x=lx (applied voltage)
%--------------------------------------------------------------------------
sigma = 1.0; beta = 0.0;                                              % electrical conductivity, potential-dependent coefficient
mCst  = [sigma beta];                                                  % material constants
g     = 0;                                                             % gravity (unused)
rho   = 1;                                                             % material density
lstps = 1;                                                             % number of loadsteps

%% 3D mesh parameters
nelsx = 20;                                                            % number of elements in x
nelsy = 20;                                                            % number of elements in y
nelsz = 20;                                                            % number of elements in z
lx = 20; ly = 20; lz = 20;                                            % domain dimensions
mp = 2;                                                                % material points per direction per element
mpType = 2;                                                            % 1 = MPM, 2 = GIMP
cmType = 1;                                                            % constitutive model (unused for electric)
mpSte  = 0;

%% Mesh generation
[etpl,coord] = formCoord3D(nelsx,nelsy,nelsz,lx,ly,lz);               % 3D background mesh
[~,nen]      = size(etpl);                                             % number of element nodes (8)
[nodes,nD]   = size(coord);                                            % number of nodes and dimensions (3)
h            = [lx ly lz]./[nelsx nelsy nelsz];                        % element lengths in each direction

%% Boundary conditions on background mesh
% Dirichlet BC: left face (x=0) grounded, right face (x=lx) applied voltage
bc = zeros(nodes, 2);
for node=1:nodes
  if coord(node,1)==0                                                  % left face: ground (phi=0)
    bc(node, :) = [node 0];
  elseif coord(node,1)==lx                                             % right face: applied voltage (phi=1V)
    bc(node, :) = [node 1];
  end
end
bc = bc(bc(:,1)>0,:);                                                  % remove empty rows

%% Mesh data structure
mesh.etpl  = etpl;
mesh.coord = coord;
mesh.bc    = bc;
mesh.h     = h;

%% Material point generation
ngp    = mp^nD;                                                        % material points per element (2^3 = 8)
GpLoc  = detMpPos(mp,nD);                                             % local MP positions (ngp,3)
N      = shapefunc(nen,GpLoc,nD);                                     % basis functions for MP positions (ngp,8)
[etplmp,coordmp] = formCoord3D(nelsx,nelsy,nelsz,lx,ly,lz);           % mesh for MP generation
nelsmp = size(etplmp,1);                                               % number of elements
nmp    = ngp*nelsmp;                                                   % total number of material points

mpC = zeros(nmp,nD);                                                   % zero MP coordinates
for nel=1:nelsmp
  indx = (nel-1)*ngp+1:nel*ngp;                                       % MP indices
  eC   = coordmp(etplmp(nel,:),:);                                    % element coordinates (8,3)
  mpPos = N*eC;                                                        % global MP coordinates (ngp,3)
  mpC(indx,:) = mpPos;
end

lp = zeros(nmp,3);                                                     % domain half-lengths
lp(:,1) = h(1)/(2*mp);
lp(:,2) = h(2)/(2*mp);
lp(:,3) = h(3)/(2*mp);
vp = 2^nD * lp(:,1).*lp(:,2).*lp(:,3);                               % volume: (2*lp_x)*(2*lp_y)*(2*lp_z)

mid = [lx/2, ly/2, lz/2];                                             % sphere center at domain center
Z = zeros(nmp,1);

%% Material point structure generation
for mpi = nmp:-1:1
  mpData(mpi).mpType = mpType;
  mpData(mpi).cmType = cmType;
  mpData(mpi).mpC    = mpC(mpi,:);                                     % 3D coordinates
  mpData(mpi).vp     = vp(mpi);
  mpData(mpi).vp0    = vp(mpi);

  mpData(mpi).phi     = zeros(1,1);                                    % electric potential
  mpData(mpi).gradPhi = zeros(nD,1);                                   % potential gradient (3D)
  mpData(mpi).J       = zeros(nD,1);                                   % current density (3D)

  mpData(mpi).mpM    = vp(mpi)*rho;
  mpData(mpi).nIN    = zeros(nen,1);                                   % nodes (8 for hex)
  mpData(mpi).eIN    = 0;
  mpData(mpi).Svp    = zeros(1,nen);                                   % basis functions (1,8)
  mpData(mpi).dSvp   = zeros(nD,nen);                                  % basis derivatives (3,8)
  mpData(mpi).Fn     = eye(3);
  mpData(mpi).F      = eye(3);
  mpData(mpi).sig    = zeros(6,1);
  mpData(mpi).epsEn  = zeros(6,1);
  mpData(mpi).epsE   = zeros(6,1);
  mpData(mpi).mCst   = mCst;
  mpData(mpi).fp     = 0;                                              % volumetric current source
  mpData(mpi).u      = zeros(nD,1);                                    % displacements (3D)
  if mpData(mpi).mpType == 2
    mpData(mpi).lp   = lp(mpi,:);
    mpData(mpi).lp0  = lp(mpi,:);
    mpData(mpi).R    = zeros(3);
  else
    mpData(mpi).lp   = zeros(1,nD);
    mpData(mpi).lp0  = zeros(1,nD);
  end
  mpData(mpi).mpSte = mpSte;
end

%% Remove material points inside spherical hole
sphereRadius = 3;                                                      % radius of spherical hole
for i=1:nmp
    r = norm(mpData(i).mpC - mid);
    if r < sphereRadius
        Z(i) = 1;
    end
end
mpData(find(Z)) = [];
end
