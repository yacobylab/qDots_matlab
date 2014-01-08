function [bestAngles bestState bestFid statefun] = PartialBellStateFinder(initstate, data,opts)
% Set up all of the states and rotators
if ~exist('opts','var')
    opts='';    
end

sigmaZ = [1 0; 0 -1]; sigmaY = [0 -1i; 1i 0]; sigmaX = [0 1;1 0];
Sigma(:,:,1) = sigmaX;
Sigma(:,:,2) = sigmaY;
Sigma(:,:,3) = sigmaZ;

psi_p = (1/sqrt(2))*[1 0 0 1]';
psi_m = (1/sqrt(2))*[1 0 0 -1]';
phi_p = (1/sqrt(2))*[0 1 1 0]';
phi_m = (1/sqrt(2))*[0 1 -1 0]';

% Pick a Bell state, and make a function to generate it.

statefun = @(angles) rotatepv(initstate, angles);
% Test it, make some PauliPlots

% PauliPlot(statefun([0 0 0 0]));

% Compute expectations for a state.

bs = {'XI','YI','ZI','IX','IY','IZ','XY', 'XZ','YX', 'YZ','ZX', 'ZY', 'XX', 'YY','ZZ'};
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

% This makes a function that evaluates a pauli vector.
if 0
    fprintf('PauliVec = @(psi) [...\n');
    for j= 1:length(bs)
        fprintf([  'psi(:)''*',bs{j},'*psi(:), ...\n']);  % Complex conj
    end
    fprintf('  ]\n');
end
% The function from above
PauliVec = @(psi) [...
psi(:)'*XI*psi(:), ...
psi(:)'*YI*psi(:), ...
psi(:)'*ZI*psi(:), ...
psi(:)'*IX*psi(:), ...
psi(:)'*IY*psi(:), ...
psi(:)'*IZ*psi(:), ...
psi(:)'*XY*psi(:), ...
psi(:)'*XZ*psi(:), ...
psi(:)'*YX*psi(:), ...
psi(:)'*YZ*psi(:), ...
psi(:)'*ZX*psi(:), ...
psi(:)'*ZY*psi(:), ...
psi(:)'*XX*psi(:), ...
psi(:)'*YY*psi(:), ...
psi(:)'*ZZ*psi(:), ...
  ];

if length(data) == 4
    data = PauliVec(data);
end

bbs=15; % Nominal bell bank size
%persistent angles;
%persistent ebb;
angles = []; ebb = [];
if isempty(ebb) || 0
    angles=[];
    ebb=[];
    if ~isempty(strfind(opts,'stonly'))
        for p1=linspace(0,2*pi,bbs);
            for p2=linspace(0,2*pi,bbs);
                angles=[angles; pi/2 0 p1 p2];
                ebb=[ebb; PauliVec(statefun(angles(end,:)))];
            end
        end
    else
        tic
        for t1 = linspace(0,pi,bbs)
            for t2 = linspace(0,pi,bbs)
                for p1 = linspace(0,2*pi,abs(sin(t1))*bbs+1)
                    for p2 = linspace(0,2*pi,abs(sin(t2))*bbs+1)
                        angles=[angles ; t1 t2 p1 p2];
                        ebb=[ebb; (statefun([t1 t2 p1 p2]))];
                    end
                end
            end
        end
    end
    toc
    fprintf('Bell bank size is %d\n',size(ebb,1));
end

% Find the State from BellBank with most overlap with the current state
fidfunc = @(states) (1+sum(repmat(data(:)',size(states,1),1)'*states,1))/4;
fid = fidfunc(ebb);
[bestFid ind] = max(fid);
bestAngles = angles(ind,:);
if ~isempty(strfind(opts,'stonly'))    
  bestAngles(3:4)=fminsearch(@(a) 1-fidfunc(PauliVec(statefun([pi/2 pi/2 a(1:2)]))),bestAngles(3:4),optimset('TolX',1e-8));
else
  bestAngles=fminsearch(@(a) 1-fidfunc((statefun(a))),bestAngles,optimset('TolX',1e-8));
end
bestState = statefun(bestAngles);
%bestPauli=PauliVec(bestState);
bestFid = fidfunc(bestState);
if ~isempty(strfind(opts,'plot'))
  pauliPlot(bestState)
end
end

function f = rotatepv(pvec, rot)

%rot = theta theta phi phi
%opts = 'STonly';
% if strfind(opts, 'STonly')
%     ML = rotationmat3D(rot(1), [1 0 0]);
%     MR = rotationmat3D(rot(4), [1 0 0]);
% else
    ML = rotationmat3d(rot(3), [0 0 1])*rotationmat3D(rot(1), [0 1 0]);
    MR = rotationmat3d(rot(4), [0 0 1])*rotationmat3D(rot(2), [0 1 0]);
%end
L = [ML' zeros(3,12)];
R = [zeros(3,3) MR' zeros(3,9)];
MM = kron(ML',MR');
M = [L; R; [zeros(9,6) MM]];
if size(pvec,1)==1
    pvec = pvec';
end
f = M*pvec;
end