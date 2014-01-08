function out = MtoChi(M)
%function out = MtoChi(M)
% return a chi matrix given an M matrix M.
persistent A;
if isempty(A)
    
    sx=[0 1;1 0];
    sy=[0 -1i;1i 0];
    sz=[1 0;0 -1];
    sI=eye(2);
    A=zeros(4,4,4,4);
    sigma = {sI,sx,sy,sz};
    
    for i = 1:4
        for j = 1:4
            for k = 1:4
                for l=1:4
                    A(i,j,k,l)= (1/2)*trace(sigma{i}*sigma{k}*sigma{j}*sigma{l});
                end
            end
        end
    end
end

if all(size(M)==3)
    d = zeros(4,4);
    d(1,1)=1;
    d(2:end,2:end)=M;
    M=d;
end

out = zeros(4,4);

for i = 1:4
    for j = 1:4
        r=squeeze(A(i,j,:,:)) .* M;
        out(i,j)=sum(r(:));
    end
end

out = (1/4)*out;


end