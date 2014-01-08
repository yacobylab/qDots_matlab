function center()
global tunedata
dv=35e-6;
for n=1:10
    autotune('stp tl tmp')
    chng(1:2,n)=abs(tunedata.measp(1:2,1)); 
    if mean(chng(:,n))>2*mean(chng(:,1))
        fprintf('This is not centering \n');
        sleep
        return
    end
    if chng(1,n)<dv && chng(2,n)<dv
        fprintf('Well centered after %01d runs \n',n);
        sleep
        return
    end
    atcenter('noconfirm'); 
end
    fprintf('Did not converge after 10 runs \n');
    sleep
end 