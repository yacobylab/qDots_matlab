function gradfeedbfn(x, data)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.



global fbdata;
ldi=length(fbdata.dataind);
newfig = ~ishandle(1034);
if ~newfig && ~ishandle(fbdata.ploth(1,1))
    close(1034);
    newfig=1;
end
    
xe=x(end);

for di=fbdata.dataind
    %cut = 10; use this if doing software average.
    %nframe = sum(isfinite(data{fbdata.dataind}(:, fbdata.pulses(1), cut)));
    %y = mean(data{fbdata.dataind}(nframe, fbdata.pulses, cut:end), 3);

    %nframe = sum(isfinite(data{fbdata.dataind}(:, fbdata.pulses(1))));
    % good for single group
    %y = data{fbdata.dataind}(x(2), fbdata.pulses);

    y = permute(mean(data{di}(xe, :, :, :), 4), [3 2 1]);
    y = y(fbdata.pulses);
    x = fbdata.x;

    fbdata.refval(di) = y(1);
    
    switch fbdata.method
        case 1 %good for long times
            y = (y-mean(y))./std(y);
            %x = 1:length(y);

            basis = [ones(size(x)); cos(2*pi*x*fbdata.freq); sin(2*pi*x*fbdata.freq)];

            coeff = y'\basis'; %matlab specific least square fit
            %coeff

            yfit = basis'*coeff';
            fbval = fbdata.freq + (mod(atan2(coeff(2), coeff(3)) - fbdata.setp + pi, 2*pi) - pi)/(2 * pi * fbdata.t0);

        case 2
            y = (y-mean(y))./std(y);
            %x = 0:length(y)-1;
            ft = ifft(y);
            ft = ft(1:round(end/2));
            [m, mi] = max(abs(ft(2:end)));


            fp(4) = 2* pi* mi/length(x);
            fp(2) = real(ft(mi+1));
            fp(3) = imag(ft(mi+1));
            %fitpar(i, 1) = y(1) - fitpar(i, 2);
            fp(1) = ft(1)+.1; % need this to condition fit.



            fp = nlinfit(x, y, @fitfn, fp);
            yfit = fitfn(fp, x);
            fbval = fp(4)/(2*pi);
            %plot(x, y, '-.', x0, fitfn(fp0, x0),':', x0, fitfn(fitpar(i, :), x0));

        case 3  % nlinfit with fft based initial  guess.
            y = y-y(1);
            %x = 0:length(y)-1;


            y = y - x * fbdata.bgslp;

            ft = ifft(y);
            ft = ft(1:round(end/2));
            [m, mi] = max(abs(ft(2:end)));


            fp(3) = 2* pi* mi/length(x);
            fp(1) = real(ft(mi+1));
            fp(2) = imag(ft(mi+1));


            fp = nlinfit(x(2:end), y(2:end), @fitfn3, fp);
            yfit = fitfn3(fp, x);
            fbval = fp(3)/(2*pi);
    end


    err = fbdata.setp - fbval;

    fbdata.fbval = [fbdata.fbval(2:end), fbval];
    fbdata.time = [fbdata.time(2:end), now*24*3600];

    if fbdata.fbon
        fbdata.intval = [fbdata.intval(2:end), fbdata.intval(end) + fbdata.igain * err * diff(fbdata.time(end-1:end))];
        fbdata.ctrlval = [fbdata.ctrlval(2:end), fbdata.pgain * err + fbdata.intval(end)];
    else
        fbdata.intval = [fbdata.intval(2:end), 0]; %[fbdata.intval(2:end), fbdata.ctrlval(end) - fbdata.pgain * err];
        fbdata.ctrlval = fbdata.ctrlval([2:end, end]);
    end


    plottime = (fbdata.time-fbdata.time(end))/60;

    if newfig   
        figure(1034);
        subplot(3,ldi,di);
        fbdata.ploth([1:2 7],di) = plot(x, y, '.', x, yfit, 'r',0,0,'rx');
        title('data');
        subplot(3,ldi,ldi+di)
        fbdata.ploth(3,di) = plot(plottime(2:ldi:end), fbdata.fbval(2:ldi:end));
        fbdata.ploth(6,di) = title(sprintf('fbval %d',fbdata.pulseind));
        subplot(3,ldi,2*ldi+di)
        fbdata.ploth(4:5,di) = plot(plottime, fbdata.intval, plottime, fbdata.ctrlval, 'r');
        title('polarization time');
    else
        set(fbdata.ploth(1,di), 'YData', y);
        set(fbdata.ploth(2,di), 'YData', yfit);
        if(~isfield(fbdata,'xval') || isempty(fbdata.xval))
            xv=-1;
        else
            xv=fbdata.xval;
        end
        set(fbdata.ploth(3,di), 'XData', plottime(2:ldi:end), 'YData', fbdata.fbval(2:ldi:end));
        sigma=std(fbdata.fbval(2:ldi:end))*2*pi;
        t2s=sqrt(2/(sigma*sigma));
        set(fbdata.ploth(4,di), 'XData', plottime, 'YData', fbdata.intval);
        set(fbdata.ploth(5,di), 'XData', plottime, 'YData', fbdata.ctrlval);
        a=fp(1)/sqrt(fp(1)^2+fp(2)^2);        
        per=2*pi/fp(3);
        tpi=(atan2(fp(2),-fp(1))+pi)/fp(3);
        tpi2=(atan2(fp(2),-fp(1))+pi/2)/fp(3);
        %tpi2=fzero(@(x) fitfn3(fp,x)-sqrt(fp(1)^2+fp(2)^2), tpi/2);
        %tpi2=min(acos(.5*(a+sqrt(3-3*a*a))),acos(.5*(a-sqrt(3-3*a*a))))/fp(3);

        if ~isfield(fbdata,'tpihist') || length(fbdata.tpihist) < di || any(isnan(fbdata.tpihist{di}))
            fbdata.tpihist{di}=[];
            fbdata.tpi2hist{di}=[];
        end
        npts=25;             
        
        if length(fbdata.tpihist{di}) >= npts  
          fbdata.tpihist{di}=[fbdata.tpihist{di}(2:npts) tpi];
          fbdata.tpi2hist{di}=[fbdata.tpi2hist{di}(2:npts) tpi2];
        else
          fbdata.tpihist{di}(end+1)=tpi;
          fbdata.tpi2hist{di}(end+1)=tpi2;
        end
        
        tpi=mean(fbdata.tpihist{di});
        tpi2=mean(fbdata.tpi2hist{di});
        set(fbdata.ploth(6,di), 'String', sprintf('Per: %.1f  \\pi: %.1f \\pi/2: %.1f T_2^* %.1f', per,tpi,tpi2,t2s));
        set(fbdata.ploth(7,di),'XData', xv);
    end
    drawnow;
end


return;

function y = fitfn(beta, x)
y = beta(1) + beta(2) * cos(beta(4) * x) + beta(3) * sin(beta(4) * x) ;

function y = fitfn3(beta, x)
y = sqrt(sum(beta(1:2).^2)) + beta(1) * cos(beta(3) * x) + beta(2) * sin(beta(3) * x) ;
