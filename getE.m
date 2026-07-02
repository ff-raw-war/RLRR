function E = getE(X, A, Z, J1, alpha, mu1)
    E = soft_threshold1(X-A*Z+J1/mu1, alpha/mu1);
end