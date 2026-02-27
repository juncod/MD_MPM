function [fint,Ke,mpData] = detMPs_Electric(phi_nodes,mpData,nD)
%Material point stiffness matrix and internal force for electric field problem
%--------------------------------------------------------------------------
% Based on: detMPs_Heat.m by William Coombs
% Governing equation: -div(sigma * grad(phi)) = 0
% Ohm's law: J = -sigma * grad(phi)
% Weak form stiffness: Ke = sum_mp (G' * D * G) * vol
%--------------------------------------------------------------------------
nmp   = length(mpData);
fint  = zeros(size(phi_nodes));
tnSMe = sum([mpData.nSMe]);
krow  = zeros(tnSMe,1); kcol=krow; kval=krow;
npCnt = 0;

for mp=1:nmp
    sigma0 = mpData(mp).mCst(1); beta = mpData(mp).mCst(2);
    nIN = mpData(mp).nIN;
    N = mpData(mp).Svp;
    dNx = mpData(mp).dSvp;
    nn  = size(dNx,2);
    ed = nIN;

    G = dNx;

    % Interpolate potential at MP (for potential-dependent conductivity)
    phiP = N * phi_nodes(ed);
    sigma_c = sigma0 + beta * phiP;             % sigma(phi) = sigma0 + beta*phi
    D = sigma_c * eye(nD);                       % isotropic conductivity tensor

    gradPhi = G * phi_nodes(ed);                 % potential gradient at MP
    J_flux  = D * gradPhi;                       % internal flux (D*gradPhi, used for weak-form residual)

    vol = mpData(mp).vp;
    kp = (G' * D * G) * vol;                     % conductance matrix contribution
    fp = (G' * J_flux) * vol;                    % internal force contribution

    mpData(mp).gradPhi = gradPhi;                % store potential gradient
    mpData(mp).J       = -D * gradPhi;           % current density: J = -sigma*gradPhi (Ohm's law)

    npDoF = nn^2;
    krow(npCnt+1:npCnt+npDoF) = repmat(ed', nn, 1);
    kcol(npCnt+1:npCnt+npDoF) = repmat(ed,  nn, 1);
    kval(npCnt+1:npCnt+npDoF) = kp(:);

    npCnt = npCnt + npDoF;
    fint(ed) = fint(ed) + fp;
end

nDoF = length(phi_nodes);
Ke = sparse(krow, kcol, kval, nDoF, nDoF);
end
