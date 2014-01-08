function [pt] = vzerocross(data)
% data is a cell-array of 3 3-dimensional tensors.
%  find the best point in space to approximate simultaneous zero crossing
%    of all tensors
%  Do this by finding voxles where a sign change occurs for each tensor
%   then interpolating.

for i=1:length(data)
    vdata{i} = data{i} > 0;
    sz=size(vdata);
    sz=sz-ones(size(sz));
    changed{i}=zeros(sz);
    switch(ndims(data{i}))
        case 1
            changed{i} = changed{i} | vdata{1}(1:end-1) ~= vdata{1}(2:end);
        case 2
            changed{i} = changed{i} | vdata{1}(1:end-1,:) ~= vdata{1}(2:end,:);
            changed{i} = changed{i} | vdata{1}(:,1:end-1) ~= vdata{1}(:,2:end);
        case 3
            changed{i} = changed{i} | vdata{1}(1:end-1,:,:) ~= vdata{1}(2:end,:,:);
            changed{i} = changed{i} | vdata{1}(:,1:end-1,:) ~= vdata{1}(:,2:end,:);            
            changed{i} = changed{i} | vdata{1}(:,:,1:end-1) ~= vdata{1}(:,:,2:end);
        otherwise
            error('Be more clever here');
    end
end
ints=ones(size(changed{1}));
for i=1:length(changed)
    ints = ints & changed{i};
end

% Ints is now 1 wherever there is an intersection
switch ndims(ints)
    case 1
        [i] = find(ints);
        for j=1:length(i)
            pt(j)=j+data{1}(j) / (data{1}(j+1)-data{1}(j));
        end    
    case 2
        [i j] = find(ints)
        for k=1:length(i)
        end
end
