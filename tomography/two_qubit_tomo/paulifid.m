function val = paulifid(p1,p2)
%function val = paulifid(p1,p2)
% return the fidelity with which p1 approximates p2.
% p1 and p2 may be a matrix with pauli vectors in columns
% ie [ 1 1 1 ; 0.5 0.5 0.5 ; 0 0 0 0 ; 0 0 0];
% fidelity is trace(sqrtm(r1)*r2*sqrtm(r1)).^2
if size(p1,1) == 1
    p1 = p1';
end
if size(p2,1) == 1
    p2 = p2';
end
val=zeros(1,size(p1,2));
for i=1:size(p1,2)
   m1=pauli2density(p1(:,i));
   sm2=sqrtm(pauli2density(p2(:,i)));   
   val(i)=real(trace(sqrtm(sm2*m1*sm2)).^2);
end
end
