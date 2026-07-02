function A = getA(X, E, Z)
    % A = (J1*Z'+mu1*P'*X*Z')/(mu1*Z*Z');
    
    PXZ=(X-E)*Z'; % d'×c
    PXZ(isnan(PXZ))=0;
    PXZ(isinf(PXZ))=0;
    [U,~,V]=svd(PXZ,'econ');
    A = U*V';
end