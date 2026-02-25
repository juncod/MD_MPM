function [lstps,g,mpData,mesh] = setupGrid_electric
%% Analysis parameters
sigma = 1.0; beta = 0.0;                                              % electrical conductivity, potential-dependent coefficient
mCst=[sigma beta];                                                     % material constants
g=0;                                                                   % gravity (unused)
rho=1;                                                                 % material density
J_source = 100;                                                        % applied current source (at right boundary MPs)
lstps=1;                                                               % number of loadsteps
a = 1;                                                                 % element multiplier
nelsx=42;                                                              % number of elements in the x direction
nelsy=40;                                                              % number of elements in the y direction
ly=40;  lx=42;                                                         % domain dimensions
d=20;  l=40;                                                           % conductor region dimensions
mp=2;                                                                  % number of material points in each direction per element
mpType = 2;                                                            % material point type: 1 = MPM, 2 = GIMP
cmType = 1;                                                            % constitutive model (unused for electric problem)
mpSte=0;
%% Mesh generation
[etpl,coord] = formCoord2D(nelsx,nelsy,lx,ly);                         % background mesh generation
[~,nen]      = size(etpl);                                             % number of element nodes
[nodes,nD]   = size(coord);                                            % number of nodes and dimensions
h            = [lx ly]./[nelsx nelsy];                                 % element lengths in each direction

%% Boundary conditions on background mesh
% Dirichlet BC: left edge (x=0) grounded (phi=0), right edge (x=lx) applied voltage (phi=1)
bc = zeros(nodes, 2);                                                  % generate empty bc matrix
for node=1:nodes                                                       % loop over nodes
  if coord(node,1)==0                                                  % left electrode: ground (phi=0)
    bc(node, :)=[node 0];
  elseif coord(node,1)==lx                                             % right electrode: applied voltage (phi=1 V)
    bc(node, :)=[node 1];
  end
end
bc = bc(bc(:,1)>0,:);                                                  % remove empty part of bc

%% Mesh data structure generation
mesh.etpl  = etpl;                                                     % element topology
mesh.coord = coord;                                                    % nodal coordinates
mesh.bc    = bc;                                                       % boundary conditions
mesh.h     = h;                                                        % mesh size

%% Material point generation
ngp    = mp^nD;                                                        % number of material points per element
GpLoc  = detMpPos(mp,nD);                                              % local MP locations
N      = shapefunc(nen,GpLoc,nD);                                      % basis functions for the material points
[etplmp,coordmp] = formCoord2D(40,20,l,d);                             % mesh for MP generation
coordmp(:,2)=coordmp(:,2)+(ly-d);                                      % adjust MP locations (vertical)
nelsmp = size(etplmp,1);                                               % no. elements populated with material points
nmp    = ngp*nelsmp;                                                   % total number of material points

mpC=zeros(nmp,nD);                                                     % zero MP coordinates
for nel=1:nelsmp
  indx=(nel-1)*ngp+1:nel*ngp;                                          % MP locations within mpC
  eC=coordmp(etplmp(nel,:),:);                                         % element coordinates
  mpPos=N*eC;                                                          % global MP coordinates
  mpC(indx,:)=mpPos;                                                   % store MP positions
end
lp = zeros(nmp,2);                                                     % zero domain lengths
lp(:,1) = h(1)/(2*mp);                                                 % domain half length x-direction
lp(:,2) = h(2)/(2*mp);                                                 % domain half length y-direction
vp      = 2^nD*lp(:,1).*lp(:,2);                                       % volume associated with each material point
count =0;
mid = [20,30];
Z = zeros(nmp,1);
%% Material point structure generation
for mp = nmp:-1:1                                                      % loop backwards over MPs so array doesn't change size
  mpData(mp).mpType = mpType;                                          % material point type: 1 = MPM, 2 = GIMP
  mpData(mp).cmType = cmType;                                          % constitutive model type (unused)
  mpData(mp).mpC    = mpC(mp,:);                                       % material point coordinates
  mpData(mp).vp     = vp(mp);                                          % material point volume
  mpData(mp).vp0    = vp(mp);                                          % material point initial volume

  mpData(mp).phi    = zeros(1,1);                                      % electric potential at MP
  mpData(mp).gradPhi = zeros(nD,1);                                    % potential gradient at MP
  mpData(mp).J      = zeros(nD,1);                                     % current density at MP (J = -sigma*gradPhi)

  mpData(mp).mpM    = vp(mp)*rho;                                      % material point mass
  mpData(mp).nIN    = zeros(nen,1);                                    % nodes associated with the material point
  mpData(mp).eIN    = 0;                                               % element associated with the material point
  mpData(mp).Svp    = zeros(1,nen);                                    % material point basis functions
  mpData(mp).dSvp   = zeros(nD,nen);                                   % derivative of the basis functions
  mpData(mp).Fn     = eye(3);                                          % previous deformation gradient (kept for compatibility)
  mpData(mp).F      = eye(3);                                          % deformation gradient (kept for compatibility)
  mpData(mp).sig    = zeros(6,1);                                      % Cauchy stress (kept for compatibility)
  mpData(mp).epsEn  = zeros(6,1);                                      % previous elastic strain (kept for compatibility)
  mpData(mp).epsE   = zeros(6,1);                                      % elastic strain (kept for compatibility)
  mpData(mp).mCst   = mCst;                                            % material constants: [sigma, beta]
  mpData(mp).fp     = 0;                                               % volumetric current source (zero: BCs drive the solution)
  mpData(mp).u      = zeros(nD,1);                                     % displacements (kept for compatibility)
  if mpData(mp).mpType == 2
    mpData(mp).lp     = lp(mp,:);                                      % material point domain lengths (GIMP)
    mpData(mp).lp0    = lp(mp,:);                                      % initial material point domain lengths (GIMP)
    mpData(mp).R = zeros(3);
  else
    mpData(mp).lp     = zeros(1,nD);                                   % material point domain lengths (MPM)
    mpData(mp).lp0    = zeros(1,nD);                                   % initial material point domain lengths (MPM)
  end
  mpData(mp).mpSte = mpSte;
end
for i=1:nmp
    r = norm(mpData(i).mpC-mid);
    if r<5
        Z(i) = 1;
    end
end
mpData(find(Z)) = [];
end
