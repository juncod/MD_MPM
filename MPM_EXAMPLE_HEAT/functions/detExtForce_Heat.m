function [fext] = detExtForce_Heat(nodes,mpData)

%Global external force determination  
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   23/01/2019
% Description:
% Function to determine the external forces at nodes based on body forces
% and point forces at material points.
%
%--------------------------------------------------------------------------
% [fbdy,mpData] = DETEXTFORCE(coord,etpl,g,eIN,mpData)
%--------------------------------------------------------------------------
% Input(s):
% nodes  - number of nodes (total in mesh)
% nD     - number of dimensions
% g      - gravity
% mpData - material point structured array. Function requires:
%           mpM   : material point mass
%           nIN   : nodes linked to the material point
%           Svp   : basis functions for the material point
%           fp    : point forces at material points
%--------------------------------------------------------------------------
% Ouput(s);
% fext   - external force vector (nodes*nD,1)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

nmp  = size(mpData,2);                                                      % number of material points & dimensions 
fext = zeros(nodes,1);                                                   % zero the external force vector
for mp = 1:nmp
   nIN = mpData(mp).nIN;                                                    % nodes associated with MP
   Svp = mpData(mp).Svp;                                                    % basis functions
   qp = mpData(mp).fp * Svp;
   fext(nIN) = fext(nIN) + qp';                                                % combine into external force vector
end
end