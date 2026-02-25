function [ddT, drct] = linSolve_Heat(bc,Kt,oobf,NRit,fd)

%Linear solver
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   23/01/2019
% Description:
% Function to solve the linear system of equations for the increment in
% displacements and reaction forces.  The linear system is only solved for
% the first Newton-Raphson iteration (NRit>0) onwards as the zeroth
% iteration is required to construct the stiffness matrix based on the
% positions of the material points at the start of the loadstep.  This is
% different from the finite element method where the stiffness matrix from
% the last iteration from the previous loadstep can be used for the zeroth
% iteration. 
%
% In the case of non-zero displacement boundary conditions, these are
% applied when NRit = 1 and then the displacements for these degrees of
% freedom are fixed for the remaining iterations.
%
%--------------------------------------------------------------------------
% [duvw,drct] = LINSOLVE(bc,Kt,oobf,NRit,fd)
%--------------------------------------------------------------------------
% Input(s):
% bc    - boundary conditions (*,2)
% Kt    - global stiffness matrix (nDoF,nDoF)
% oobf  - out of balance force vector (nDoF,1)
% NRit  - Newton-Raphson iteration counter (1)
% fd    - free degrees of freedom (*,1)
%--------------------------------------------------------------------------
% Ouput(s);
% duvw  - displacement increment (nDoF,1)
% drct  - reaction force increment (nDoF,1)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

nDoF = length(oobf);                                                        % number of degrees of freedom 
ddT = zeros(nDoF,1);                                                       % zero displacement increment
drct = zeros(nDoF,1);                                                       % zero reaction increment
if (NRit)>0                                                             
    fixedNodes = bc(:,1);
    fixedValues = bc(:,2);
    
    ddT(fixedNodes) = (1+sign(1-NRit)) * fixedValues;
    ddT(fd) = Kt(fd,fd) \ (oobf(fd) - Kt(fd, fixedNodes) * ddT(fixedNodes));
    drct(fixedNodes) = Kt(fixedNodes,:) * ddT - oobf(fixedNodes);
end
end