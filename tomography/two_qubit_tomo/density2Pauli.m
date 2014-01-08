function pauli = density2Pauli(density)
%function pauli = density2Pauli(density)
%returns the pauli plot for a given density matrix
%bs = {'XI','YI','ZI','IX','IY','IZ','XY', 'XZ','YX', 'YZ','ZX', 'ZY', 'XX', 'YY','ZZ'};

sigmaX = [0 1;1 0]; sigmaY = [0 -1i;1i 0]; sigmaZ = [1 0;0 -1];
XI = kron(sigmaX, eye(2));
YI = kron(sigmaY, eye(2));
ZI = kron(sigmaZ, eye(2));
IX = kron(eye(2), sigmaX);
IY = kron(eye(2), sigmaY);
IZ = kron(eye(2), sigmaZ);
XY = kron(sigmaX,sigmaY);
XZ = kron(sigmaX,sigmaZ);
YX = kron(sigmaY,sigmaX);
YZ = kron(sigmaY,sigmaZ);
ZX = kron(sigmaZ, sigmaX);
ZY = kron(sigmaZ, sigmaY);
XX = kron(sigmaX,sigmaX);
YY = kron(sigmaY,sigmaY);
ZZ = kron(sigmaZ, sigmaZ);

bs = {XI,YI,ZI,IX,IY,IZ,XY,XZ,YX, YZ,ZX, ZY, XX, YY,ZZ};

pauli = zeros(1,15);

for j = 1:length(bs)
    pauli(j) = trace(density*bs{j});
end

end