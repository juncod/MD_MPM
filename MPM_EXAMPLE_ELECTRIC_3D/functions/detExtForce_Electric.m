function [fext] = detExtForce_Electric(nodes,mpData)

%Global external current source determination
%--------------------------------------------------------------------------
% Author: William Coombs (adapted for electric field problem)
% Date:   23/01/2019
% Description:
% Function to determine the external current sources at nodes based on
% volumetric current sources at material points (Neumann BC contribution).
%
%--------------------------------------------------------------------------
% [fext] = DETEXTFORCE_ELECTRIC(nodes,mpData)
%--------------------------------------------------------------------------
% Input(s):
% nodes  - number of nodes (total in mesh)
% mpData - material point structured array. Function requires:
%           nIN   : nodes linked to the material point
%           Svp   : basis functions for the material point
%           fp    : volumetric current source at material point (A/m^3)
%--------------------------------------------------------------------------
% Ouput(s):
% fext   - external current source vector (nodes,1)
%--------------------------------------------------------------------------

nmp  = size(mpData,2);                                                 % number of material points
fext = zeros(nodes,1);                                                 % zero the external force vector
for mp = 1:nmp
   nIN = mpData(mp).nIN;                                               % nodes associated with MP
   Svp = mpData(mp).Svp;                                               % basis functions
   Jp = mpData(mp).fp * Svp;                                           % current source at nodes
   fext(nIN) = fext(nIN) + Jp';                                        % assemble into external force vector
end
end
