function F = getF(Z)
    
    [c, n]=size(Z);
    Z(n, n) = 0; % Z扩展为n×n
    %{
    LZ = (Z(1:c)+Z(1:c)')/2;
    %}
    LZ = (Z+Z')/2;
    LZ(isnan(LZ))=0;
    LZ(isinf(LZ))=0;
    [F, ~, ev] = eig1(LZ, c, 0);
end