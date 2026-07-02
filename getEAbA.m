function E = getEAbA(X, Z, J1, alpha, mu1)
    E = soft_threshold1(X-X*Z+J1/mu1, alpha/mu1);
end