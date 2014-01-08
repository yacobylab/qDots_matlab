%%
%chans = {'1a', '2a', '3a', '4a', '1b', '2b', '3b', '4b'};
chans={'SD1top','SD1mid','SD4top','SD4mid','SD1bot','SD4bot'};
counter = 0;
tic;
while 1
    smset(chans, -.01);
    val = cell2mat(smget(chans));
    if any(abs(val+.01) > 1e-3) %|| counter > 1000
        'failure'
        return
    else
        counter = counter+1;        
    end
    if ~mod(counter,100)
        fprintf('%2d %2f \n', counter,toc/60);
    end
end