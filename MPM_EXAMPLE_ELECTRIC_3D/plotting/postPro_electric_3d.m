%Post processing script for 3D electric field MPM code
%--------------------------------------------------------------------------
% Author: William Coombs (adapted for 3D electric field problem)
% Description:
% The script produces VTK output files based on the background mesh and
% material point data (electric potential and current density) in 3D.
%--------------------------------------------------------------------------
% See also:
% MAKEVTK                  - VTK file for background mesh
% MAKEVTKMP_ELECTRIC_3D    - VTK file for 3D electric field MP data
%--------------------------------------------------------------------------

mpDataName = sprintf('output/mpData_%i.vtk',lstp);                    % MP output data file name
phi_mp = [mpData.phi]';                                                % electric potential at all MPs (nmp,1)
J_mp   = reshape([mpData.J],nD,nmp)';                                  % current density at all MPs (nmp,3)
mpC    = reshape([mpData.mpC],nD,nmp)';                                % material point coordinates (nmp,3)
makeVtkMP_Electric_3d(mpC, phi_mp, J_mp, mpDataName);                  % generate material point VTK file
if lstp==0
    makeVtk(mesh.coord,mesh.etpl,'output/mesh.vtk')                    % generate mesh VTK file
end
