function beta1 = fitwrap_o(ctrl, x, y, beta0, model, mask)
% beta1 = fitwrap(ctrl, x, y, beta0, model, mask)
% ctrl: plinit, plfit, woff, nofit, pause
%plinit plots the initial guess
%plfit plots the fit
%woff turns off warnings
%no fit does not fit the data
%pause will pause after each plot
% white : whiten fit parameters roughly

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

n = size(y, 1);
if size(x, 1) == 1
    x = repmat(x, n, 1);
end

if size(beta0, 1) == 1
    beta0 = repmat(beta0, n, 1);
end

if nargin < 6
    mask = true(1, size(beta0, 2));
end

if(strfind(ctrl, 'white'))
 bs = beta0;
 bs(bs < 1e-10) = 1;
else
 bs = beta0*0+1;
end

beta0=beta0./bs;
beta1 = beta0;

if strfind(ctrl, 'woff')
    ws(1) = warning('query', 'stats:nlinfit:IllConditionedJacobian');
    ws(2) = warning('query', 'stats:nlinfit:IterationLimitExceeded');
    ws2 = ws;
    [ws2.state] = deal('off');
    warning(ws2);
end

for i = 1:n
    if strfind(ctrl, 'pl')
        figure(500);
        clf;
        hold on;
    end

    if strfind(ctrl, 'plinit')
        plot(x(i, :), y(i, :), '.-', x(i, :), model(beta0(i, :).*bs(i,:), x(i, :)), 'r--');
    end
    
    if strfind(ctrl, 'nofit')
        continue
    end

    beta1(i, mask) = nlinfit(x(i, :), y(i, :), @(p,x)fitfn(p.*bs(i,:),x), beta0(i, mask));

    if strfind(ctrl, 'plfit')
        plot(x(i, :), y(i, :), '.-', x(i, :), model(beta1(i, :).*bs(i,:), x(i, :)), 'k');
    end
    
    if ~isempty(strfind(ctrl, 'pause')) && i < n
        pause
    end      
  beta1(i,:) = beta1(i,:).*bs(i,:);
end

if strfind(ctrl, 'woff')
    warning(ws);
end


    function y = fitfn(beta, x)
        beta([find(mask), find(~mask)]) = [beta, beta0(~mask)];
        y = model(beta, x);        
    end

end


