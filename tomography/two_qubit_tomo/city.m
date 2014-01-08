%function []=city(pdata,pideal)
%pdata=[1 0 0 -1 0 0 1 0 0 0 0 0 0 0 0 -1];
%pideal=[1 0 0 -1 0 0 1 0 0 0 0 0 0 0 0 -1];
%rdata=pauli2density(pdata);
%rideal=pauli2density(pideal);
rdata=[0 0 0 0; 0 0.05 0 0; 0 0 .95 0; 0 0 0 0];
rideal=[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];
figure(5)
clf
h=bar3(rdata,.5);
%hold on
shading interp
colormap(winter);
  for i = 1:length(h)
        zdata = get(h(i),'ZData');
        set(h(i),'CData',zdata)
        % Add back edge color removed by interpolating shading       
        set(h,'EdgeColor','k')
        set(h,'FaceAlpha',.5);
    end
    set(gca,'DataAspectRatio',[1 1 .5]);
    axis tight;      
set(gca,'XTickLabel',{'00';'01';'10';'11'})
set(gca,'YTickLabel',{'00';'01';'10';'11'})
set(gca,'ZLim',[-1 1]);
set(gca,'CLim',[-1 1]);

hold on
h=bar3(rideal,.50001);
        set(h,'EdgeColor','k')
        set(h,'FaceAlpha',0);
   
%end