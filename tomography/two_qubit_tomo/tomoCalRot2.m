function val = tomoCalRot(Pauli, ML, MR, opts)
%function val = tomoCalRot(Pauli, transformation, opts)
%remakes pauli plots of variable Pauli (15-compenent vector) given the two
%compensation matricies ML and MR (left and right).
% for now, opts gets passed to PauliPlot (so 'samefig' is handy)
%data is to be organized in rows

if ~exist('opts','var')
    opts='';
end

if size(Pauli,1) ==15
    Pauli = Pauli';
end

%reshape entangled part of the matrix into
%  7  8  9 10 11 12 13 14 15
% XX XY XZ YX YY YZ ZX ZY ZZ
Pauli = Pauli([1:6 13,7,8,9,14,10,11,12,15]);

%       1    2    3    4    5    6    7     8    9    10    11    12    13   14   15
bs = {'XI','YI','ZI','IX','IY','IZ','XY', 'XZ','YX', 'YZ','ZX', 'ZY', 'XX', 'YY','ZZ'};


L = [ML zeros(3,12)];
R = [zeros(3,3) MR zeros(3,9)];
MM = kron(ML,MR);
M = [L; R; [zeros(9,6) MM]];

%data = M*Pauli';
data = Pauli*M;
%unshape it back into normal order:
data = data([1:6 8:10, 12:14, 7, 11,15]);

% newXX = ML*Pauli([13,9,11]).*MR*Pauli([13,7,8]);
% newXY = ML*Pauli([7 14 12]).*MR*Pauli([13 7 8]);
% newXZ = ML*Pauli([11 10 15]).*MR*Pauli([13 7 8]);
% 
% newYX = ML*Pauli([13 9 11]).*MR*Pauli([9 14 12]);
% newYY = ML*Pauli([7 14 12]).*MR*Pauli([9 14 10]);

%PauliPlot(data, opts);

val = data;
%val = Pauli*M';