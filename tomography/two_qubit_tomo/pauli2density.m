function val = pauli2density(pauli)
%function val = pauli2density(initstate, rotangles)
%returns the density matrix corresponding to the given pauli vector
%pauli can be have length 3 or 4 (corresponding to 1 qubit)
%pauli can also have length 15 or 16 (2qubit)
% bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};

if size(pauli,1)>1
    pauli = pauli';
end

if length(pauli)==15 || length(pauli) == 3
    pauli = [1 pauli];
end

sx = [0 1;1 0]; sy = [0 -1i; 1i, 0]; sz = [1 0;0 -1];
Sigma = {eye(2), sx, sy, sz};

switch length(pauli)
    case 16
        rho = zeros(4);
        SSS = {};
        for j = 1:length(Sigma)
            for k = 1:length(Sigma)
                SSS{end+1} = kron(Sigma{k},Sigma{j});
            end
        end
        %SSS = {SSS{[1:5, 9, 13, 7:8, 10,12,14:15, 6,11, 16]}};
        SSS = {SSS{[1:5,9,13,10,14,7,15,8,12,6,11,16]}};
        for j = 1:length(SSS)
            rho = rho+SSS{j}*pauli(j);
        end
    case 4
        rho = [pauli(1) + pauli(4), pauli(2)-1i * pauli(3); ...
               pauli(2) + 1i*pauli(3), pauli(1)-pauli(4)];
    otherwise
        error('bad length of pauli vector!');
end
val = rho/length(rho);
end
