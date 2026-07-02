function F = getF_(Z)
    [c, n]=size(Z);
    LZ = (Z(1:c)+Z(1:c)')/2;
    [F, ~, ev] = eig1(LZ, c, 0);
end