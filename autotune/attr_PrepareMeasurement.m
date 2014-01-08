function attr_PrepareMeasurement()
%%
% Recenter the dot and lock the gradient.
global tunedata;
pause(.5);
awgcntrl('on start wait err raw');
pause(.5);
done=0;
while ~done
    autotune('stp tl tmp quiet');
    done=all(tunedata.measp(1:2,1) < .05e-3);
    
    if all(tunedata.measp(1:2,1) < .5e-3)
        atcenter('noconfirm quiet');
    else
        fprintf('Warning: large change.  Will wait for user confirmation\n');
        atcenter;
    end
    pause(1);
end
 if ~sm_setgradient
     error('Unable to lock gradient.  Is fbdata.gradtgt correct?\n');
 end
end

function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
end
