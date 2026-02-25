function [ddphi, drct] = linSolve_Electric(bc,Ke,oobf,NRit,fd)

%Linear solver for electric potential increment
%--------------------------------------------------------------------------
% Author: William Coombs (adapted for electric field problem)
% Date:   27/08/2020
% Description:
% Function to solve the linear system of equations for the increment in
% electric potential and reaction currents.  The linear system is only
% solved from the first Newton-Raphson iteration (NRit>0) onwards as the
% zeroth iteration is required to construct the conductance matrix based
% on the positions of the material points at the start of the loadstep.
%
% In the case of non-zero potential boundary conditions (e.g. applied
% voltage), these are applied when NRit = 1 and then held fixed for the
% remaining iterations.
%
%--------------------------------------------------------------------------
% [ddphi,drct] = LINSOLVE_ELECTRIC(bc,Ke,oobf,NRit,fd)
%--------------------------------------------------------------------------
% Input(s):
% bc    - boundary conditions (*,2): col1=node ID, col2=prescribed voltage
% Ke    - global conductance matrix (nDoF,nDoF)
% oobf  - out of balance force vector (nDoF,1)
% NRit  - Newton-Raphson iteration counter (1)
% fd    - free degrees of freedom (*,1)
%--------------------------------------------------------------------------
% Ouput(s):
% ddphi - electric potential increment (nDoF,1)
% drct  - reaction current increment (nDoF,1)
%--------------------------------------------------------------------------

nDoF = length(oobf);                                                   % number of degrees of freedom
ddphi = zeros(nDoF,1);                                                 % zero potential increment
drct = zeros(nDoF,1);                                                  % zero reaction current increment
if (NRit)>0
    fixedNodes = bc(:,1);
    fixedValues = bc(:,2);

    ddphi(fixedNodes) = (1+sign(1-NRit)) * fixedValues;
    ddphi(fd) = Ke(fd,fd) \ (oobf(fd) - Ke(fd, fixedNodes) * ddphi(fixedNodes));
    drct(fixedNodes) = Ke(fixedNodes,:) * ddphi - oobf(fixedNodes);
end
end
