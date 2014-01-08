function out = unrotateEntangle(dataexp)

if size(dataexp,2) == 15
    dataexp = dataexp';
end

if 0
sigmaZ = [1 0; 0 -1]; sigmaY = [0 -1i; 1i 0]; sigmaX = [0 1;1 0];
Sigma(:,:,1) = sigmaX;
Sigma(:,:,2) = sigmaY;
Sigma(:,:,3) = sigmaZ;

psi_p = (1/sqrt(2))*[1 0 0 1]';
psi_m = (1/sqrt(2))*[1 0 0 -1]';
phi_p = (1/sqrt(2))*[0 1 1 0]';
phi_m = (1/sqrt(2))*[0 1 -1 0]';

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

% pvfunc=@(p,x) [cos((x(:)+p(2))*p(3)) * p(1), ...  % <Y,I> - like coefficent
%         cos((x(:)+p(2))*p(3)) * p(1), ...  % <I,Y> - like coefficent
%         repmat(p(1)^2,length(x),1),   ...  % <Y,Y> - like coefficent
%         sin((x(:)+p(2))*p(3)) * p(1), ... % <Y,ST> - like coefficent
%         sin((x(:)+p(2))*p(3)) * p(1)]; % <ST,Y> - like coefficent
%     
% p = [1 0 10]; %init amplitude, phase, w_ent
% v = [2 5 14 8 11]; 
% psetnodec = zeros(size(dataexp));
% wfid1 = []; wfid2 = [];
% stnrm = [];
% xv = .01*(1:64);
% for m = 1:64
%     paulivec = [0 1 0 0 1 0 0 0 0 0 0 0 0 1 0]; %initial state
%     rr =pvfunc(p, xv(m));
%     paulivec(v) = rr;
%     psetnodec(:,end+1) = paulivec;
% end % make the non-dephasing, non-rotated pauli set

% angs = zeros(6,1+size(dataexp,2));
% angs(:,1) = [pi/6 0.1 0.1 pi/6 0.1 0.1]';
% pfid = zeros(1,size(dataexp,2));
% for j = 1:size(dataexp,2)
%    angs(:,j+1) =findstate(dataexp(:,j),psetnodec(:,j),angs(:,j),'');
% end
end

%params = init amplitude, w_ent, init phase, T2_1, T2_2, 
%p(6+2i) = angle1_i, 
%p(7+2i) = angle2_i
ffunc = {'YI','IY','XZ','ZX','YY'};
xv = linspace(0,.64,64); %work in microseconds
p0 = [.8, .1, pi/6, .5, .5];
for j = 1:size(dataexp,2)
    mfdata(j).x = ones(size(dataexp,1),1)*xv(j);
    mfdata(j).y = dataexp(:,j);
    mmodel(j).fn = @paulifit;
    off=2*(j-1);
    mmod(j).pt=@(p) [ p(1:5), p(6+off) p(7+off)];
    p0 = [p0, pi/6, pi/6]; %make this more intelligent?
end

mask = p0./p0; mask(4:5) = [0 0];
[out.p out.chisq, out.cov]=mfitwrap(mfdata,mmodel,p0,'plinit plfit plotiter err',mask);
[out.p out.chisq, out.cov]=mfitwrap(mfdata,mmodel,out.p,'plinit plfit plotiter err');

out.dataexp = dataexp;

end

function val = findstate(pv,stateguess, rotguess,opts)
rotfid = @(data, initguess, rotations, popts) fidfunc(data,...
    rotatepv(initguess, rotations,popts));
val = fminsearch(@(rr)(1-rotfid(pv,stateguess,rr,opts)),rguess);

if ~isempty(strfind(opts, 'STonly'))
    val([2:3]) = [pi/2 0];
    val([5:6]) = [pi/2 0];
end
end

function paulivec = paulifit(p,x)
%bs = {'XI','YI','ZI','IX','IY','IZ','XY', 'XZ','YX', 'YZ','ZX', 'ZY', 'XX', 'YY','ZZ'};
x = x(1);
%params = init amplitude, w_ent, init phase, T2_1, T2_2, angle1, angle2
PP = [p(1)*cos(x*p(2)+p(3)).*exp(-(x./p(4)));... % YI
    p(1)*cos(x*p(2)+p(3)).*exp(-(x./p(5)));...  %IY
    p(1)*sin(x*p(2)+p(3)).*exp(-(x./p(4)));... %XZ
    p(1)*sin(x*p(2)+p(3)).*exp(-(x./p(5)));... %ZX
    p(1)*exp(-(x./p(4))).*exp(-(x./p(5)))]; % YY
paulivec = zeros(15,1);    
% PP = [p(1)*cos(x(:)*p(2)+p(3)).*exp(-(x(:)./p(4))),... % YI
%     p(1)*cos(x(:)*p(2)+p(3)).*exp(-(x(:)./p(5))),...  %IY
%     p(1)*sin(x(:)*p(2)+p(3)).*exp(-(x(:)./p(4))),... %XZ
%     p(1)*sin(x(:)*p(2)+p(3)).*exp(-(x(:)./p(5))),... %ZX
%     p(1)*exp(-(x(:)./p(4))).*exp(-(x(:)./p(5)))]; % YY

% PP =[cos((x+p(2))*p(3)) * p(1), ...  % <Y,I> - like coefficent
%         cos((x+p(2))*p(3)) * p(1),...  % <I,Y> - like coefficent
%         repmat(p(1)^2,length(x),1),   ...  % <Y,Y> - like coefficent
%         sin((x(:)+p(2))*p(3)) * p(1), ... % <Y,ST> - like coefficent
%         sin((x(:)+p(2))*p(3)) * p(1)]; % <ST,Y> - like coefficent
%paulivec = zeros(64,15);
 
v = [2 5 14 8 11]; 
paulivec(v) = PP;
% for j = 1:size(paulivec,1)
%   paulivec(j,:) = PauliRotator(paulivec(j,:), [p(6), 0, 0, p(7), 0, 0]);
% end
   
end

function val = YI(p,x)
  val = p(1)*cos(x*p(2)+p(3)).*exp(-(x./p(4)));
end
function val = IY(p,x)
  val = p(1)*cos(x*p(2)+p(3)).*exp(-(x./p(5)));
end
function val = XZ(p,x)
  val = p(1)*sin(x*p(2)+p(3)).*exp(-(x./p(4)));
end
function val = ZX(p,x)
  val = p(1)*sin(x*p(2)+p(3)).*exp(-(x./p(5)));
end
function val = YY(p,x)
  val = p(1)*exp(-(x./p(4))).*exp(-(x./p(5)));
end
function val = OTHS(p,x)
  val = zeros(size(x));
end

function val = paulifit2(p,x)
  %bs = {'XI','YI','ZI','IX','IY','IZ','XY', 'XZ','YX', 'YZ','ZX', 'ZY', 'XX', 'YY','ZZ'};
val = [OTHS(p,x), YI(p,x), OTHS(p,x), OTHS(p,x), IY(p,x), OTHS(p,x), OTHS(P,x),...
    XZ(p,x), OTHS(p,x), OTHS(p,x), ZX(p,x), OTHS(p,x), OTHS(p,x),OTHS(p,x),YY(p,x), OTHS(p,x)];
end
function val = paulifit3(p,x)
  
end

function out = PauliRotator(PauliVec, rotations)
%rotations = %rot = [rotation angle, axis_theta, axis_phi] for each qubit
scramble=[1:6 8:10 12:14 7 11 15];
PauliVec = PauliVec([1:6 13 7 8 9 14 10 11 12 15]);
ML = rotationmat3D(rotations(1), [sin(rotations(2))*cos(rotations(3)),...
    sin(rotations(2))*sin(rotations(3)), cos(rotations(2))]);

MR = rotationmat3D(rotations(4), [sin(rotations(5))*cos(rotations(6)),...
    sin(rotations(2))*sin(rotations(6)), cos(rotations(5))]);

L = [ML' zeros(3,12)];
R = [zeros(3,3) MR' zeros(3,9)];
MM = kron(ML',MR');
M = [L; R; [zeros(9,6) MM]];
if size(PauliVec,1)>1
    d = M*PauliVec;
    out = d(scramble);
else
    d= (M*PauliVec')';
    out = d(scramble);
end

end
