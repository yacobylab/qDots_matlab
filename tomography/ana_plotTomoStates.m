function ana_plotTomoStates(data,tvs,opts)
% function ana_plotTomoStates(data,tvs,opts)
% Data is a q-long cell array of 3xn data 
%   tvs is a q-long vector.  Cell entries with the same tv are put on the
%   same subplot.
if ~exist('opts','var')
    opts=struct();
end
opts=def(opts,'opts','');
opts=def(opts,'fignum',256);
tvv=unique(tvs);
spr=max(floor((length(tvv)^.45)),1);
spc=ceil(length(tvv)/spr);
figure(opts.fignum);
clf;
for j=1:length(tvv)
    subplot(spr,spc,j);
    bloch_sphere([],'hold');
    colors='rgbcmkgrgbcmykrgbcmyk';
    %colors='rgbcmrgbcmrgbcmrgbcmrgbcmrgbcm';   
    inds=find(tvs==tvv(j));
    for in=1:length(inds);
        if isopt(opts,'discrete')
            i=inds(in);
            hold on;
            for q=1:size(data{i},2)            
              p=plot3(data{i}(1,q),data{i}(2,q),data{i}(3,q),[colors(mod(q,end-1)+1) 'x']);
              set(p,'LineWidth',2);
            end
        else
            i=inds(in);
            hold on;
            q=1;
            p=plot3(data{i}(1,q),data{i}(2,q),data{i}(3,q),[colors(in) 'o']);
            set(p,'LineWidth',2);
            q=1:size(data{i},2);
            plot3(data{i}(1,q),data{i}(2,q),data{i}(3,q),[colors(in) '.-']);
            p=plot3(data{i}(1,end),data{i}(2,end),data{i}(3,end),[colors(in) 'x']);
            set(p,'LineWidth',2);
        end
    end
    h(j)=gca;
end


hlink=linkprop(h,{'CameraPosition','CameraUpVector'});
key = 'graphics_linkprop';
setappdata(h(1),key,hlink); 
return;



% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      s=setfield(s,f,v);
  end
return;

function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
return;