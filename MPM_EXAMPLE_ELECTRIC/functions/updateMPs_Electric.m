function [mpData] = updateMPs_Electric(phi_nodes, mpData)

%Material point update: electric potential
%--------------------------------------------------------------------------
% Author: William Coombs (adapted for electric field problem)
% Date:   29/01/2019
% Description:
% Function to update the electric potential at each material point by
% interpolating from the converged nodal potential values using the
% material point basis functions.
%
%--------------------------------------------------------------------------
% [mpData] = UPDATEMPS_ELECTRIC(phi_nodes,mpData)
%--------------------------------------------------------------------------
% Input(s):
% phi_nodes - nodal electric potential (nodes,1)
% mpData    - material point structured array.  The function requires:
%              - nIN : nodes associated with material point
%              - Svp : basis functions
%--------------------------------------------------------------------------
% Ouput(s);
% mpData    - material point structured array.  The function modifies:
%              - phi : electric potential at material point
%--------------------------------------------------------------------------

nmp = length(mpData);                                                  % number of material points
for mp=1:nmp
    nIN = mpData(mp).nIN;                                              % nodes associated with material point
    N   = mpData(mp).Svp;                                              % basis functions
    mpData(mp).phi = N * phi_nodes(nIN);                               % interpolate nodal potential to MP
end
end
