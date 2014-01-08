%grad = tunedata.runs(28).grad([1 2 5 6 9:11], 1:2)
ng = length(tunedata.gatechan);
tunedata.basis = eye(ng);

%%
clear g;
g(1, :) = {'3b', '3a',  '4a'}; %Lead3
g(2, :) = {'4b', '3a', '4a'}; %Lead4
g(3, :) = {'T34', '3a', '4a'}; %T34
g(4, :) = {'N34', '3a', '4a'}; %N34
%g(5, :) = {'2a', '2b', '1a'};  %2a/2b
g(5, :) = {'3a', '4a',  [1 0]}; %X
g(6, :) = {'3a', '4a',  [0 1]}; %Y

m = (cellfun(@ischar, g));
[dummy, gind] = ismember(unique(smchanlookup(g(m(:)))), tunedata.gatechan);
if ~all(dummy) 
    error('Gate missing.');
end

tunedata.basenames = repmat({''}, 1, ng);
tunedata.basenames(gind)= {'Lead3', 'Lead4', 'T34', 'N34', 'X', 'Y'};
%% assumes size(g, 1) == length(gind), 
%i.e. number of directions = number  of gates varied.
grad = tunedata.runs(end).grad(:, 1:2);

tunedata.basis(gind, gind) = 0;
for i = 1:size(g, 1) 
    if ischar(g{i, 3})
        [dummy, g2] = ismember(smchanlookup(g(i, :)), tunedata.gatechan); 
        tunedata.basis(g2, gind(i)) = [1, - grad(g2(1), :)/grad(g2(2:3), :)]';
    else
        [dummy, g2] = ismember(smchanlookup(g(i, 1:2)), tunedata.gatechan);
        tunedata.basis(g2, gind(i)) = (g{i, 3}/grad(g2, :))';
    end
end

tunedata.runs(end).basis = tunedata.basis;
%% new gradient code including many dots
ng=length(tunedata.gatechan);
tunedata.basis=eye(ng);
clear g;
b={};
% If all entries are gates, try to control first gate w/out changing mu in
% dots.
% If last entry is a vector, try to use gates specified to change mu in
% that direction
g(1,:) ={'1b',  '1a','2a','3a','4a'}; b{end+1}='Lead1';
g(2,:) ={'2b',  '1a','2a','3a','4a'}; b{end+1}='Lead2';
g(3,:) ={'3b',  '1a','2a','3a','4a'}; b{end+1}='Lead3';
g(4,:) ={'4b',  '1a','2a','3a','4a'}; b{end+1}='Lead4';
g(5,:) ={'T12', '1a','2a','3a','4a'}; b{end+1}='T12';
g(6,:) ={'T34', '1a','2a','3a','4a'}; b{end+1}='T34';
g(7,:) ={'N12', '1a','2a','3a','4a'}; b{end+1}='N12';
g(8,:) ={'N34', '1a','2a','3a','4a'}; b{end+1}='N34';
g(9,:) ={'1a','2a','3a','4a',[1 0 0 0]}; b{end+1}='XL';
g(10,:)={'1a','2a','3a','4a',[0 1 0 0]}; b{end+1}='YL';
g(11,:)={'1a','2a','3a','4a',[0 0 1 0]}; b{end+1}='XR';
g(12,:)={'1a','2a','3a','4a',[0 0 0 1]}; b{end+1}='YR';
g(13,:)={'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4',[0 -1 0 0]}; b{end+1}='pYL';
g(14,:)={'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4',[-1 0 0 0]}; b{end+1}='pXL';
g(15,:)={'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4',[0 0 -1 0]}; b{end+1}='pXR';
g(16,:)={'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4',[0 0 0 -1]}; b{end+1}='pYR';

m = (cellfun(@ischar, g));
[dummy, gind] = ismember(unique(smchanlookup(g(m(:)))), tunedata.gatechan);
if ~all(dummy) 
    error('Gate missing.');
else
    fprintf('Gates ok\n');
end

tunedata.basenames = repmat({''}, 1, ng);
tunedata.basenames(gind)= b;
%% makes basis from gradient on both dots.  make sure autotune('grad') or 'grad copy' has been run.  (or resp).

ileft=find(strcmp({tunedata.sets.name},'left'));
iright=find(strcmp({tunedata.sets.name},'right'));
if isfield(tunedata.sets(ileft).runs(end),'grad') && ~isempty(tunedata.sets(ileft).runs(end).grad)
  lgrad=tunedata.sets(ileft).runs(end).grad(:,1:2);
else
  lgrad=zeros(length(tunedata.gatechan),2);
end

if isfield(tunedata.sets(iright).runs(end),'grad') && ~isempty(tunedata.sets(iright).runs(end).grad)
  rgrad=tunedata.sets(iright).runs(end).grad(:,1:2);
else
  rgrad=zeros(length(tunedata.gatechan),2);  
end

grad=[lgrad rgrad];
%%
tunedata.basis(gind, gind) = 0;
for i = 1:size(g, 1)     
    if ischar(g{i, end})   % Change gate
        [dummy, g2] = ismember(smchanlookup(g(i, :)), tunedata.gatechan);         
        tunedata.basis(g2, gind(i)) = [1, - grad(g2(1), :)* pinv(grad(g2(2:end), :))]';
    else  % Change mu
        [dummy, g2] = ismember(smchanlookup(g(i, 1:end-1)), tunedata.gatechan);
%        tunedata.basis(g2, gind(i)) = (g{i, end}/grad(g2, :))';
        tunedata.basis(g2, gind(i)) = (pinv(grad(g2, :)')*g{i,end}');
    end
end

tunedata.runs(end).basis = tunedata.basis;
for l=1:length(tunedata.sets)
    tunedata.sets(l).runs(end).basis = tunedata.basis;
end
%% check gradient in new basis
grad=[tunedata.sets(ileft).runs(end).grad(:,1:2) tunedata.sets(iright).runs(end).grad(:,1:2)]

tunedata.basis' * tunedata.sets(1).runs(end).grad
tunedata.basis' * tunedata.sets(2).runs(end).grad


%% compute compensation matrix
for i = 1:2
    ind(i) = find(~cellfun(@isempty, {tunedata.sets(i).runs.resp}), 1, 'last');
end


%strvcat(tunedata.sets(1).runs(ind(1)).resp.gx(15:16).name,tunedata.sets(2).runs(ind(2)).resp.gx(15:16).name)

gates = strvcat(tunedata.sets(1).runs(ind(1)).resp.gx(13:16).name, tunedata.sets(2).runs(ind(2)).resp.gx(13:16).name);
gateind = gates(:, end)-'0';
[gi1, gi1] = sort(gateind(1:4));
[gi2, gi2] = sort(gateind(5:8));

grad = [tunedata.sets(2).runs(ind(2)).resp.tpgrady(13:16);
    tunedata.sets(2).runs(ind(2)).resp.tpgradx(13:16);
    tunedata.sets(1).runs(ind(1)).resp.tpgradx(13:16);
    tunedata.sets(1).runs(ind(1)).resp.tpgrady(13:16)];
% row ordering hardcoded!


grad(1:2, :) = grad(1:2, gi2);
grad(3:4, :) = grad(3:4, gi1);

inv(-grad)

%% Make compensation matrix
ileft=find(strcmp({tunedata.sets.name},'left'));
iright=find(strcmp({tunedata.sets.name},'right'));
if isfield(tunedata.sets(ileft).runs(end),'grad') && ~isempty(tunedata.sets(ileft).runs(end).grad)
  lgrad=tunedata.sets(ileft).runs(end).grad(13:16,1:2);
else
  lgrad=zeros(length(tunedata.gatechan),2);
end

if isfield(tunedata.sets(iright).runs(end),'grad') && ~isempty(tunedata.sets(iright).runs(end).grad)
  rgrad=tunedata.sets(iright).runs(end).grad(13:16,1:2);
else
  rgrad=zeros(length(tunedata.gatechan),2);  
end

lgrad=lgrad';
rgrad=rgrad';
gradmat=-[lgrad(2,:); lgrad(1,:); rgrad];
compmatrix=pinv(gradmat);


%% This cell saves the old compmatrix with today's date and makes a new one
side=tunedata.name;
atswap('right'); atswap('left'); atswap(side);

load compmatrix;
compmatrix_old = compmatrix;
fname = ['compmatrix_' datestr(now,29)]; % see the help for why 29 
save(fname,'compmatrix');

ileft=find(strcmp({tunedata.sets.name},'left'));
iright=find(strcmp({tunedata.sets.name},'right'));
if isfield(tunedata.sets(ileft).runs(end),'grad') && ~isempty(tunedata.sets(ileft).runs(end).grad)
  lgrad=tunedata.sets(ileft).runs(end).grad(13:16,1:2); %rows 13-16 of grad refer to pls ramps 1-4
else
  lgrad=zeros(length(tunedata.gatechan),2);
  fprintf('Left gradient is empty');
end

if isfield(tunedata.sets(iright).runs(end),'grad') && ~isempty(tunedata.sets(iright).runs(end).grad)
  rgrad=tunedata.sets(iright).runs(end).grad(13:16,1:2);
else
  rgrad=zeros(length(tunedata.gatechan),2);  
  fprintf('Right gradient is empty');
end

lgrad=lgrad';
rgrad=rgrad';
gradmat=-[lgrad(2,:); lgrad(1,:); rgrad];
compmatrix=pinv(gradmat);
fprintf('change in compmatrix \n');
compmatrix-compmatrix_old
save compmatrix compmatrix
