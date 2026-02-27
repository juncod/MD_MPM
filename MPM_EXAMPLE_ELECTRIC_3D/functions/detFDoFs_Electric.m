function [fd] = detFDoFs_Electric(mesh)
%Determine the free degrees of freedom on the background mesh
%--------------------------------------------------------------------------
% Author: William Coombs (adapted for electric field problem)
% Date:   17/12/2018
% Description:
% Function to determine the free degrees of freedom of the background mesh
% based on the elements that contain material points and the voltage
% boundary conditions (Dirichlet BCs for electric potential).
%
%--------------------------------------------------------------------------
% [fd] = DETFDOFS_ELECTRIC(mesh)
%--------------------------------------------------------------------------
% Input(s):
% mesh  - mesh structured array. Function requires:
%           - etpl : element topology (nels,nen)
%           - eInA : elements "active" in the analysis
%           - bc   : boundary conditions (*,2): node ID and prescribed voltage
%--------------------------------------------------------------------------
% Ouput(s):
% fd    - free degrees of freedom on the background mesh (*,1)
%--------------------------------------------------------------------------

incN   = unique(mesh.etpl(mesh.eInA>0,:));                             % unique active node list
incDoF = incN';                                                        % active degrees of freedom

[nodes,~] = size(mesh.coord);                                          % no. nodes and dimensions
fd     = (1:nodes);                                                    % all degrees of freedom
fd(mesh.bc(:,1)) = 0;                                                  % zero fixed voltage BCs

fd     = fd(incDoF);                                                   % only include active DoF
fd     = fd(fd>0);                                                     % remove fixed voltage BCs

end
