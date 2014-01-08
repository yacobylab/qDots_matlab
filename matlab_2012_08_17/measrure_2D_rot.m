    epsvals = linspace(5,3,32);

    for j = 27:length(epsvals)
    smset('RFpow3', 10)
     fprintf('starting eps_offset = %.2f\n',epsvals(j));
    %find the freq
    scangen = fConfSeq2_v2([28+j],struct('nloop',75,'nrep',64, 'opts','pol','datachan','DAQ2'));
    
    scangen.loops(2).npoints = 25;
    if 1% j==1
    scangen.loops(1).rng = [190 260]*1e6;
    else
        scangen.loops(1).rng = goodfreq+25e6*[-1 1];
    end
    scangen.loops(1).setchan = {'RFfreq3'};
    scangen.loops(2).setchan  = 'count';
    scangen.loops(1).npoints = 100;

    if sm_setgradient;
    fprintf('dbz = %.2f\n',sm_getgradient);
    fprintf('Voffset =%.2f\n',epsvals(j));
    scanname = smnext(sprintf('ramsey_gen_rot_%.2f_R',epsvals(j)));
    d=smrun(scangen,scanname);
    if any(isnan(d{1}))
       break 
    end
    rehash path;
    goodfreq = find_rot_frame_freq(['sm_',scanname,'.mat']);
    %scangen.loops(1).rng = goodfreq + 25e6*[1 -1];
    smset('RFfreq3',goodfreq);
    fprintf('setting freq to %.2fMHz \n',1e-6*goodfreq);
    fprintf('dbz = %.2f \n',sm_getgradient);

    scangen2 = fConfSeq2_v2(28+j,struct('nloop',75,'nrep',100, 'opts','pol','datachan','DAQ2'));
    %scangen2 = fConfSeq2_v2([28, (28+j)*ones(1,99)],struct('nloop',100,'nrep',100, 'opts','pol','datachan','DAQ2'));
    scangen2.loops(2).npoints = 50;
    scangen2.loops(1).rng = [0 10];
    scangen2.loops(1).setchan = {'RFpow3'};
    scangen2.loops(2).setchan  = 'count';
    scangen2.loops(1).npoints = 100;

    scanname = smnext(sprintf('ramsey_gen_amp_%.2f_R',epsvals(j)));
    
    if sm_setgradient
       d = smrun(scangen2,scanname); 
       fprintf('dbz = %.2f \n',sm_getgradient);
       
    else
        fprintf('no lock gradient %i \n',j)
    end

    else
        sleep
       fprintf('no lock gradient for vrf = %.2f \n',j); 
    end



    end
    sleep