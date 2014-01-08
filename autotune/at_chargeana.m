function [trip, slp] = at_chargeana(scan, data,opts)
% function [trip, slp] = at_chargeana(scan, data)
% Analyze a charge scan to find triple points and slopes.

global tunedata;
d=data{1};
rng=vertcat(scan.loops.rng);

d = diff(d, [], 2);
d = d- median(d(:));
d = d .* sign(mean(d(:)));

thresh = 2 * std(d(:)); %changed from 3 sigma 09/20/10

ix = 1:size(d, 2);
jy = 1:size(d, 1);
x = linspace(rng(1, 1), rng(1, 2), size(d, 2)+1);
x = 0.5 * (x(2:end)+x(1:end-1));
y = linspace(rng(2, 1), rng(2, 2), size(d, 1));

plotderiv(rng,d);

% Automatic triple point identification
trip2=[];
autofit=0;
if isempty(strfind(opts,'man')) && isfield(tunedata.chrg,'imgl') && ~isempty(tunedata.chrg.imgl) && ~isempty(tunedata.chrg.imgr)
  manual = ~isempty(strfind(opts,'mnslp')); % Manual slope mode?
  try
      [tripx,tripy] = at_correlate(scan,data);
      autofit=1;
      trip2=[tripx', tripy'];
      if trip2(1,1) < trip2(2,1)
          fprintf('Correlator found bl triple point right of tr\n');
          trip2=[];  % Bad result; ignore
          autofit=0;
      end
      if trip2(1,2) < trip2(2,2)
          fprintf('Correlator found bl triple point above tr\n');
          trip2=[];  % Bad result; ignore
          autofit=0;
      end
  catch
      trip2=[];
      autofit=0;
  end  
end
if isempty(trip2)    
    while(size(trip2,1) < 2)
        plotderiv(rng,d);        
        manual=0;
        fprintf('Click on R and L triple point.  n to negate data.  m to begin manual fit.\n');
        while(size(trip2,1) < 2)
            [xp, yp, kp] = ginput(1);
            switch(kp)
                case 1
                    trip2 = [trip2 ; [xp, yp]];
                case 'n'
                    d=-d;
                    break;
                case 'm'
                    manual=~manual;
                    fprintf('Manual mode is %d\n',manual);
            end
        end
    end
end

if(trip2(1,1) < trip2(2,1)) % User is a doofus and clicked l-r rather than r-l
  trip2 = circshift(trip2,1);
  fprintf('You are a doofus\n');
end

if(manual)
   fprintf('Click on horizontal transition out of L near edge of plot\n');
   hl=ginput(1);
   fprintf('Click on vertical transition out of L near edge of plot\n');
   vl=ginput(1);
   fprintf('Click on horizontal transition out of R near edge of plot\n');
   hr=ginput(1);
   fprintf('Click on vertical transition out of R near edge of plot\n');
   vr=ginput(1);
   trip = trip2([1 3 2 4]);
   slp = [(hl(2)-trip2(2,2))/(hl(1)-trip2(2,1)), ...
          (hr(2)-trip2(1,2))/(hr(1)-trip2(1,1)), ...
          (vl(2)-trip2(2,2))/(vl(1)-trip2(2,1)), ...
          (vr(2)-trip2(1,2))/(vr(1)-trip2(1,1))]; 
    plot([trip2(1,1),vr(1)],[trip2(1,2),vr(2)],'k-',...
         [trip2(1,1),hr(1)],[trip2(1,2),hr(2)],'k-',...
         [trip2(2,1),vl(1)],[trip2(2,2),vl(2)],'k-',...
         [trip2(2,1),hl(1)],[trip2(2,2),hl(2)],'k-');
    return;
end
cntr = mean(trip2);
slp=diff(trip2,[],1);
slp_trip=slp(2)/slp(1);
slp_perp=-1/slp_trip;
[x2, y2] = meshgrid(x, y);
above=(y2 > (x2-cntr(1))*slp_trip + cntr(2));  % Above 11-02 transition
right=(x2 > cntr(1)+slp_perp*(y2-cntr(2)));    % right of 11-02 transition
above_r=(y2 > (x2-trip2(1,1))*slp_perp + trip2(1,2)); % above or right of tr
below_l=(y2 < (x2-trip2(2,1))*slp_perp + trip2(2,2)); % below or left of bl

if 0 % Debug image segmentation
    figure(7);
    subplot(221)
    scale=std(d(:));
    imagesc(rng(1, :), rng(2, :), d+scale*above);
    title('above');
    axis xy;
    subplot(222);
    imagesc(rng(1, :), rng(2, :), d+scale*right);
    title('right');
    axis xy;
    subplot(223)
    scale=std(d(:));
    imagesc(rng(1, :), rng(2, :), d+scale*above_r);
    title('above_r');
    axis xy;
    subplot(224);
    imagesc(rng(1, :), rng(2, :), d+scale*below_l);
    title('below_l');
    axis xy;
end

%UR sector
pmt=find((d(:) > thresh) & ( above(:)) & ( right(:)) & ( above_r(:)));
pmr=find((d(:) > thresh) & (~above(:)) & ( right(:)) & ( above_r(:)));
pml=find((d(:) > thresh) & ( above(:)) & (~right(:)) & ( below_l(:)));
pmb=find((d(:) > thresh) & (~above(:)) & (~right(:)) & ( below_l(:)));

figure(1);
plot(cntr(1), cntr(2), 'x' , ...
     x2(pmt),y2(pmt),'^k', ...
     x2(pmr),y2(pmr),'>k', ...
     x2(pml),y2(pml),'<k', ...
     x2(pmb),y2(pmb),'vk');

robust = 1;
if length(pmt) > 2 && robust
    fit1a = fliplr(robustfit(y2(pmt), x2(pmt))');
else
    fit1a = polyfit(y2(pmt),x2(pmt),1);
end
if length(pmr) > 2 && robust
    fit1b = fliplr(robustfit(y2(pmr), x2(pmr))');
else
    fit1b = polyfit(y2(pmr),x2(pmr),1);
end
if length(pml) > 2 && robust
    fit2a = fliplr(robustfit(y2(pml), x2(pml))');
else
    fit2a = polyfit(y2(pml),x2(pml),1);
end
if length(pmb) > 2 && robust
    fit2b = fliplr(robustfit(y2(pmb), x2(pmb))');
else
    fit2b = polyfit(y2(pmb),x2(pmb),1);
end
% find intersections;
df1 = fit1a-fit1b;
df2 = fit2a-fit2b;

if autofit % Use fit triple point
  trip = [[fit1a(1) , 1] * -df1(2)/df1(1) + [fit1a(2), 0], ...
    [fit2a(1) , 1] * -df2(2)/df2(1) + [fit2a(2), 0]];
else
  trip = trip2([1 3 2 4]);
end

slp = [fit1a(1), fit2b(1), fit1b(1), fit2a(1)];


plot(...
    polyval(fit1a, [y(jy(end)), trip2(1, 2)]), [y(jy(end)), trip2(1, 2)], 'k', ...
    polyval(fit1b, [trip2(1, 2), y(jy(1))]), [trip2(1, 2), y(jy(1))], 'k',...
    polyval(fit2a, [trip2(2, 2), y(jy(end))]), [trip2(2, 2), y(jy(end))], 'k', ...
    polyval(fit2b, [y(jy(1)), trip2(2, 2)]), [y(jy(1)), trip2(2, 2)], 'k', trip([1 3]), trip([2 4]), '.k');
return;
end

function plotderiv(rng,d)
global tunedata;
figure(1);
clf;
imagesc(rng(1, :), rng(2, :), d);
set(gca, 'ydir', 'normal')
axis image;
box on;
hold on;
% Add 0,2 and 1,1 labels
dict=pdload(tunedata.name);
if dict.sep.val(1) > dict.sep.val(2)
    l11=[max(rng(1,:)) min(rng(2,:))];
    l02=[min(rng(1,:)) max(rng(2,:))];
else
    l02=[max(rng(1,:)) min(rng(2,:))];
    l11=[min(rng(1,:)) max(rng(2,:))];
end
lambda=0.9;
l11=lambda*l11+(1-lambda)*l02;
l02=lambda*l02+(1-lambda)*l11;
text(l11(1),l11(2),'(1,1)','HorizontalAlignment','Center','FontWeight','Bold','FontName','Comic Sans MS');
text(l02(1),l02(2),'(0,2)','HorizontalAlignment','Center','FontWeight','Bold','FontName','Comic Sans MS');

end
