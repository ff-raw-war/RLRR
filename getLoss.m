function loss = getLoss(Z, E, W, F, LZ, alpha, beta)
    % 目标函数是：

    % 在计算loss的时候，不用管s.t.约束项
    [c, n] = size(Z);
    WZ = W.*Z;
    [~, SigmaWZ, ~] = svd(WZ);
    SigmaWZ(n,n) = 0;
    loss = trace(SigmaWZ) + alpha*norm(E, 1) + beta*trace(F'*LZ*F);
    
end