function Z = get Z(X, E, A, W, S, H, J1, J2, J3, beta, mu1, mu2, mu3, Z0)
    [c, d] = size(Z0);
    [d, n] = size(X);
    es1 = ones(c, 1);
    es2 = ones(n,1);
    H(c,n)=0; % 扩展H
    % Z = pinv(2*eye(c)+mu1*A'*B'*B*A+mu2*es1*es1')*(A'*B'*J+mu1*A'*B'*X+mu2*es1*es2');
    Z = (mu1*A'*A+mu2*es1*es1'+mu3*W(1:c,1:c).*W(1:c,1:c)+beta*H)\(-es1*J2'+A'*J1+mu1*A'*X-mu1*A'*E+mu2*es1*es2'-J3.*W+mu3*S.*W);
    Z(isnan(Z))=0;
    Z(isinf(Z))=0;
end