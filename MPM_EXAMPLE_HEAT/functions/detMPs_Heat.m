function [fint,Kt,mpData] = detMPs_Heat(T_nodes,mpData,nD)
nmp   = length(mpData);
fint  = zeros(size(T_nodes)); 
tnSMe = sum([mpData.nSMe]); 
krow  = zeros(tnSMe,1); kcol=krow; kval=krow;
npCnt = 0;

for mp=1:nmp
    k0 = mpData(mp).mCst(1); alpha = mpData(mp).mCst(2);
    nIN = mpData(mp).nIN;  
    N = mpData(mp).Svp;
    dNx = mpData(mp).dSvp;
    nn  = size(dNx,2);
    ed = nIN; 
    
    G = dNx;
   

    Tp = N * T_nodes(ed);
    k_current = k0 + alpha * Tp;
    D = k_current * eye(nD);
    
    gradT = G * T_nodes(ed);
    flux  = D * gradT;
    
    vol = mpData(mp).vp;
    kp = (G' * D * G) * vol;     
    fp = (G' * flux) * vol;
    
    mpData(mp).gradT = gradT;
    mpData(mp).flux  = flux;
    
    npDoF = nn^2;
    krow(npCnt+1:npCnt+npDoF) = repmat(ed', nn, 1);
    kcol(npCnt+1:npCnt+npDoF) = repmat(ed,  nn, 1);
    kval(npCnt+1:npCnt+npDoF) = kp(:);
    
    npCnt = npCnt + npDoF;
    fint(ed) = fint(ed) + fp;
end

nDoF = length(T_nodes);
Kt = sparse(krow, kcol, kval, nDoF, nDoF);
end