function [fp,ff]=fitosc(x,y,opts,rng,style)
% function [fp,ff]=fitosc(x,y,opts,rng,style)
%range determines the time considered. 
%opts are: fitdecay, nocenter, plot
% q=fitosc(1:200,a.d{1}','fitdecay nocenter plotdata')
% fprintf('T2*=%g',1/q(6));
% fp are mean, amp x phase, amp y phase, freq, decay center, decay rate
   % initialization function
   fig=gcf;
   fifn.fn = @fioscill;
   fifn.args = {1};
 
   %cosfn and cosfn4 are the same...
   cosfn  = '@(y,x) y(1)+y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)'; %sin w/ phase and offset
   cosfn2 = '@(y,x) y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)'; %sin w/ phase, offset, gaussian decay
   cosfn3 = '@(y,x) y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*(y(5)*x)'; %linear decay
   cosfn4 = '@(y,x) y(1)+y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)'; %sin w/ phase and offset
   cosfn5 = '@(y,x) y(1)+y(2)*cos(y(4)*x+y(3)).*exp(-(x-y(5)).^2 * y(6).^2)'; %deals w/ phase differently, has offset decay, gaussian decay
      
   if exist('rng','var') && ~isempty(rng)
       pts=x>rng(1) & x <rng(2);
       x=x(pts);
       y=y(pts);
   end

   if isempty(strfind(opts,'fitdecay')) || (~isempty(strfind(opts,'afitdecay')) && std(y) < 2e-2)     
     fifn.args={2}; % No decay
     fp=fitwrap('fine',x,y,fifn, cosfn5, [1 0 1 1 0 0]);
     ig = [fp(1), fp(2)*cos(fp(3)), fp(2)*sin(fanap(3)), fp(4:6)];
     fp=fitwrap('fine',x,y,ig, cosfn, [1 1 1 1 0 0]);
     ff=str2func(cosfn);         
   elseif isempty(strfind(opts,'nocenter'))
     fifn.args={2}; % Decay and center
     fp=fitwrap('fine',x,y,fifn, cosfn5, [1 0 1 1 0 0]);
     ig = [fp(1), fp(2)*cos(fp(3)), fp(2)*sin(fp(3)), fp(4:6)];
     fp=fitwrap('fine',x,y,fifn,cosfn2, [1 1 1 1 0 0]);
     fp=fitwrap('fine',x,y,fp, cosfn2, [1 1 1 1 1 1]);
     ff=str2func(cosfn2);
   else      % Decay but no center
     fifn.args={2};
     fp=fitwrap('fine',x,y,fifn, cosfn5, [1 0 1 1 0 0]);
     fp = [fp(1), fp(2)*cos(fp(3)), -fp(2)*sin(fp(3)), fp(4:6)];
     fp=fitwrap('fine plfit',x,y,fp, cosfn2, [1 1 1 1 0 1]);
     ff=str2func(cosfn2);
   end
   if ~isempty(strfind(opts,'plot'))
       figure(fig);
       hold on;
       if ~isempty(strfind(opts,'interp'))
           x=linspace(x(1,1),x(1,end),512);
       end
       if ~exist('style','var') || isempty('style')           
         style='r';
       end
         plot(x,ff(fp,x),[style '-']);
       if ~isempty(strfind(opts,'plotdata'))
          plot(x,y,[style 'x']);
       end
   end
return;