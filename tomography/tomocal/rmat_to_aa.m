function [axis angle values]=rmat_to_aa(m)
%function [axis angle values]=rmat_to_aa(m)
% Given a rotation matrix, convert it to axis angle.
%   note there is a axis->-axis, angle->-angle ambiguity here.
% The eigenvaluse of m are blown out to unit norm so dephased
% rotation matrices should give reasonable results.
[e v]=eig(m);
values=diag(v);
% Make all the eigenvectors unit
v=diag(v);
v=v./abs(v);
%
rv=find(imag(v) == 0 & real(v) ~= 0,1);
if length(rv) ~= 1
    v
    error('Pathological case; whine at Oliver\n');    
end
axis=e(:,rv);
m=e*diag(v)*e';

ct=(trace(m)-1)/2;
[t mi]=max(axis);  % Bigest element in axis
a=-m+m';
switch mi
    case 1
        st=a(2,3)/(2*axis(1));
    case 2
        st=-a(1,3)/(2*axis(2));
    case 3
        st=a(1,2)/(2*axis(3));
end
angle=atan2(st,ct);

if angle < 0
   angle = -1*angle;
   axis = -1*axis;
end
axis=axis';
if nargout == 0
    fprintf('Axis: (%g,%g,%g)  Angle: %g pi\m',axis(1),axis(2),axis(3),angle/pi)
end
return;
