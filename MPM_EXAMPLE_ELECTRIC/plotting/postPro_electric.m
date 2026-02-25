%Post processing script for electric field MPM code
%--------------------------------------------------------------------------
% Author: William Coombs (adapted for electric field problem)
% Date:   29/01/2019
% Description:
% The script produces VTK output files based on the background mesh and
% material point data (electric potential and current density).
% Background mesh is plotted for lstp = 0.
%
%--------------------------------------------------------------------------
% See also:
% MAKEVTK               - VTK file for background mesh
% MAKEVTKMP_ELECTRIC    - VTK file for electric field MP data
%--------------------------------------------------------------------------

mpDataName = sprintf('output/mpData_%i.vtk',lstp);                    % MP output data file name
phi_mp = [mpData.phi]';                                                % electric potential at all MPs (nmp,1)
J_mp   = reshape([mpData.J],nD,nmp)';                                  % current density at all MPs (nmp,nD)
mpC    = reshape([mpData.mpC],nD,nmp)';                                % material point coordinates (nmp,nD)
makeVtkMP_Electric(mpC, phi_mp, J_mp, mpDataName);                     % generate material point VTK file
if lstp==0
    makeVtk(mesh.coord,mesh.etpl,'output/mesh.vtk')                    % generate mesh VTK file
end
