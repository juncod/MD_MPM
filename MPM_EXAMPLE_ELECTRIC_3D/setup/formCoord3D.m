function [etpl,coord] = formCoord3D(nelsx,nelsy,nelsz,lx,ly,lz)

%Three dimensional finite element grid generation
%--------------------------------------------------------------------------
% Author: William Coombs (extended to 3D)
% Date:   27/02/2026
% Description:
% Function to generate a 3D finite element grid of linear hexahedral
% elements (8-noded bricks).
%
% Node numbering: loops over x first, then y, then z.
% Element topology: 8-noded hexahedral with node ordering consistent
% with the shapefunc.m convention for nen=8, nD=3.
%
%--------------------------------------------------------------------------
% [etpl,coord] = FORMCOORD3D(nelsx,nelsy,nelsz,lx,ly,lz)
%--------------------------------------------------------------------------
% Input(s):
% nelsx - number of elements in the x direction
% nelsy - number of elements in the y direction
% nelsz - number of elements in the z direction
% lx    - length in the x direction
% ly    - length in the y direction
% lz    - length in the z direction
%--------------------------------------------------------------------------
% Ouput(s);
% etpl  - element topology (nels,8)
% coord - nodal coordinates (nodes,3)
%--------------------------------------------------------------------------

nels  = nelsx*nelsy*nelsz;                                                 % number of elements
nodes = (nelsx+1)*(nelsy+1)*(nelsz+1);                                     % number of nodes

%% node generation
coord = zeros(nodes,3);                                                    % zero coordinates
node  = 0;                                                                 % zero node counter
for k=0:nelsz
  z = lz*k/nelsz;
  for j=0:nelsy
    y = ly*j/nelsy;
    for i=0:nelsx
      node = node+1;
      x = lx*i/nelsx;
      coord(node,:) = [x y z];
    end
  end
end

%% element generation
% Node ordering for 8-noded hexahedral element:
%   Nodes 1-4: bottom face (z=z_lo), nodes 5-8: top face (z=z_hi)
%   Within each face: same ordering as 2D quad
%
% shapefunc.m convention (nen=8, nD=3):
%   N1 = (1-xi)(1-eta)(1-zeta)/8   -> (xi-,eta-,zeta-)
%   N2 = (1-xi)(1-eta)(1+zeta)/8   -> (xi-,eta-,zeta+)
%   N3 = (1+xi)(1-eta)(1+zeta)/8   -> (xi+,eta-,zeta+)
%   N4 = (1+xi)(1-eta)(1-zeta)/8   -> (xi+,eta-,zeta-)
%   N5 = (1-xi)(1+eta)(1-zeta)/8   -> (xi-,eta+,zeta-)
%   N6 = (1-xi)(1+eta)(1+zeta)/8   -> (xi-,eta+,zeta+)
%   N7 = (1+xi)(1+eta)(1+zeta)/8   -> (xi+,eta+,zeta+)
%   N8 = (1+xi)(1+eta)(1-zeta)/8   -> (xi+,eta+,zeta-)
%
% With our node numbering (x fastest, then y, then z):
%   Given element at grid position (ix,iy,iz) (1-based):
%   n1 = base node at (ix-1, iy-1, iz-1)
%   base = (iz-1)*(nelsx+1)*(nelsy+1) + (iy-1)*(nelsx+1) + ix
%
%   Bottom face (z=z_lo, iz-1):
%     n_xm_ym = base
%     n_xp_ym = base + 1
%     n_xm_yp = base + (nelsx+1)
%     n_xp_yp = base + (nelsx+1) + 1
%   Top face (z=z_hi, iz):
%     offset = (nelsx+1)*(nelsy+1)
%     same pattern + offset
%
% Mapping to shapefunc convention:
%   N1(xi-,eta-,zeta-) -> n_xm_ym_zm = base
%   N2(xi-,eta-,zeta+) -> n_xm_ym_zp = base + planeOff
%   N3(xi+,eta-,zeta+) -> n_xp_ym_zp = base + 1 + planeOff
%   N4(xi+,eta-,zeta-) -> n_xp_ym_zm = base + 1
%   N5(xi-,eta+,zeta-) -> n_xm_yp_zm = base + rowOff
%   N6(xi-,eta+,zeta+) -> n_xm_yp_zp = base + rowOff + planeOff
%   N7(xi+,eta+,zeta+) -> n_xp_yp_zp = base + rowOff + 1 + planeOff
%   N8(xi+,eta+,zeta-) -> n_xp_yp_zm = base + rowOff + 1

etpl = zeros(nels,8);                                                      % zero element topology
nel  = 0;                                                                  % zero element counter
rowOff   = nelsx+1;                                                        % row offset (y direction)
planeOff = (nelsx+1)*(nelsy+1);                                            % plane offset (z direction)

for iz=1:nelsz
  for iy=1:nelsy
    for ix=1:nelsx
      nel  = nel+1;
      base = (iz-1)*planeOff + (iy-1)*rowOff + ix;
      etpl(nel,1) = base;                                                 % N1: (x-,y-,z-)
      etpl(nel,2) = base + planeOff;                                      % N2: (x-,y-,z+)
      etpl(nel,3) = base + 1 + planeOff;                                  % N3: (x+,y-,z+)
      etpl(nel,4) = base + 1;                                             % N4: (x+,y-,z-)
      etpl(nel,5) = base + rowOff;                                        % N5: (x-,y+,z-)
      etpl(nel,6) = base + rowOff + planeOff;                             % N6: (x-,y+,z+)
      etpl(nel,7) = base + rowOff + 1 + planeOff;                         % N7: (x+,y+,z+)
      etpl(nel,8) = base + rowOff + 1;                                    % N8: (x+,y+,z-)
    end
  end
end
end
