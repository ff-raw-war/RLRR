function [A, Z, E] = train(X, alpha, beta, ep1, ep2, c, epsilon, maxIter)
    % tic
    lossVec = [];
    iterVec = [];
    errVec1 = [];
    errVec2 = [];
    errVec3 = [];

    rho = 1.01; 
	% rho = 3; 
    mu1 = 0.1; 
    mu2 = 0.1;
    mu3 = 0.1;
    maxMu = 10^5; 
    threshold = 2*10^-4; 
    [d, n] = size(X);
    [anchor, ind2, score] = VDA(X', c);
    
    A = anchor';
    %{
    scatter(X(:,1), X(:,2), 'b.') % Blue dots
    hold on
    scatter(anchor(:,1), anchor(:,2), 'r+') % Red plus signs
    %}
    Z = rand(c, n);
    W = ep2./(abs(Z)+ep1); 
    S = rand(c, n);
    F = rand(n, c);
    E = rand(d, n);
    LZ = rand(n, n);
    % A = rand(d, c);
    J1 = rand(d, n); 
    J2 = rand(n, 1);
    J3 = rand(c, n);
    es1 = ones(c, 1);
    es2 = ones(n, 1);

    err1 = 10*threshold; 
    errVec1 = [errVec1 err1];
    err2 = 10*threshold;
    errVec2 = [errVec2 err2];
    err3 = 10*threshold;
    errVec3 = [errVec3 err3];
    loss = getLoss(Z, E, W, F, LZ, alpha, beta);
    lossVec = [lossVec ; loss]; 
    iter = 1;
    iterVec = [iterVec iter];
	notConverged = 1;
    % initialTime = toc
    while notConverged
        iter = iter + 1;
        loss0 = loss;

        % tic;
        
        S = getS(W.*Z, J3, 1, mu3);
        E = getE(X, A, Z, J1, alpha, mu1);
        F = getF(Z);
        H = L2_distance_1(F',F');
        Z = getZ(X, E, A, W, S, H, J1, J2, J3, beta, mu1, mu2, mu3, Z);
        W = ep2./(abs(Z)+ep1); 
        A = getA(X, E, Z);
        LZ = getL(Z);
        % tic;
        
        J1 = J1 + mu1*(X-A*Z-E);
        J2 = J2 + mu2*(Z'*es1-es2);
        J3 = J3 + mu3*(W.*Z-S);
        % mu = ALF;
        mu1 = min(rho*mu1, maxMu); 
        mu2 = min(rho*mu2, maxMu);
        mu3 = min(rho*mu3, maxMu);
        
        
        err1 = norm(X-A*Z,'fro');
        err2 = norm(Z'*es1-es2, 'fro');
        err3 = norm(W.*Z-S, 'fro');
		% err = max(max(abs(P'*X-P'*X*C-E)));

		loss = getLoss(Z, E, W, F, LZ, alpha, beta);
        errVec1 = [errVec1 err1];
        errVec2 = [errVec2 err2];
        errVec3 = [errVec3 err3];
        lossVec = [lossVec loss];
        iterVec = [iterVec iter];
		
        
		if iter >= maxIter
			notConverged = 0;
        end
        
        %{
       
		if err1 < threshold && err2<threshold && err3<threshold
			notConverged = 0;
		end
        %}
		
        
        if ((loss0-loss)'*(loss0-loss) < epsilon) || (loss > loss0)
            notConverged = 0;
        end % 
        
        
        % afterTime = toc
    end % while iter
    
    
    
    %plot(hengzhou, zongzhou, 'd-');
    if iter < 10
        for i = 1:10-iter
            iterVec = [iterVec iter+i];
            lossVec = [lossVec min(lossVec)];
            % errVec = [errVec err];
        end
    end
    plot(iterVec, lossVec, 'd-');
    xlabel('Number of iterations');
    ylabel('Log of the objective value');
    

end
