function [basis basis2]= atxyfix2(stepsize)
%function val = atxyfix(stepsize)
% opts have not been added yet. can be 'repbasis' or 'avgbasis' to update the basis by replacing it
% with a new one, or adding to the running average. 
% updates the XL, YL, XR, YR part of the basis by recomputing in terms of
% itself. 

global tunedata;
global smdata;


gates = smget(1:18);

if ~exist('stepsize','var') || isempty(stepsize)
  vstep = .5e-3;  
else
    if stepsize > 1e-2
        vstep = 1e-3*stepsize;
       warning('detected huge step. Assuming you meant mV'); 
    else
       vstep = stepsize; % hard coded step size
    end
end


atswap('left');
atswap('right');
s=[{'left'},{'right'}];
order=[{'XL'},{'YL'},{'XR'},{'YR'}];

for j=1:2
    atswap(s{j})
    tunedata.chrg.scan.consts(2).val=awgseqind(tunedata.chrg.pls); %PulseLine
    smset({tunedata.chrg.scan.consts.setchan}, [tunedata.chrg.scan.consts.val]); %un-obsolete
    scan = smscanpar(tunedata.chrg.scan, tunedata.cntr);
    for i=0:length(order)
        if i>0
            atchg(order{i},vstep);            
            fprintf('atchg '); fprintf(order{i}); fprintf(' by %.2f mV \n' ,vstep*1e3);
        end
        data = smrun(scan);
        if any(isnan(data{1}(:))); smset(1:18, cell2mat(gates)); return; end

        figure(70); clf; hold on;
        imagesc(scan.loops(1).rng, scan.loops(2).rng, diff(data{1}));
        xlim([scan.loops(1).rng(1), scan.loops(1).rng(2)]);
        ylim([scan.loops(2).rng(1), scan.loops(2).rng(2)]);
        triple=ginput(2);
        trippt(i+1,1:2,j)=mean(triple,1);       
        dtrippoint=trippt(i+1,1:2,j)-trippt(1,1:2,j);        
        if i>0
            atchg(order{i},-vstep);
            fprintf('atchg '); fprintf(order{i}); fprintf(' by %.2f mV \n' ,-vstep*1e3);
            fprintf('triple point moves %.3f mV in x direction, %.3f mV in y direction for a change %.2f mV \n',dtrippoint(1)*1e3,dtrippoint(2)*1e3,vstep*1e3)
            
        end
    end
end

%Now we have all the changes, and need to create a new basis from them. 
%trip point has format (run,xy,side)

lgrad=trippt(:,:,1)-repmat(trippt(1,:,1),size(trippt,1),1); %this gives the gradient for the left trip. pt. in the format [dxl/dXL dyl/dXL; dxl/dYR dyl/dYR...]
lgrad=lgrad(2:end,:)./vstep;
rgrad=trippt(:,:,2)-repmat(trippt(1,:,2),size(trippt,1),1);
rgrad=rgrad(2:end,:)./vstep;

grad=[lgrad rgrad];

grad=grad';
basis=pinv(grad);


basisprint=basis';
basis2=tunedata.basis(1:4,9:12)*basis;
  
newbasis=[{'XLn'},{'YLn'},{'XRn'},{'YRn'}];
oldbasis=[{'XL'},{'YL'},{'XR'},{'YR'}];
basis2print=basis2';
fprintf('       ')
for i=1:4    
            fprintf(['%-4s',repmat('%-5s', 1, 5), '\n'], '',oldbasis{i});
end
               fprintf('\n');
         fprintf('-----------------------------------------------------\n')
        for i = 1:4
            fprintf(['%-9s:', repmat('%8.2g', 1, 5), '\n'], newbasis{i}, basisprint(i,1:4));
            fprintf('\n')
        end      

        fprintf('       ')
fprintf(['%-10s',repmat('%-8s', 1, 5), '\n'], '',smdata.channels(1:4).name);
         fprintf('\n');
         fprintf('-----------------------------------------------------\n')
        for i = 1:4
            fprintf(['%-9s:', repmat('%8.3g', 1, 5), '\n'], newbasis{i}, basis2print(i,1:4));
            fprintf('\n')
        end
        
               
   doit = strcmp(input('accept? (yes/[no]', 's'), 'yes');
    
    if doit
      tunedata.basis(1:4,9:12)=basis2;  
    end
