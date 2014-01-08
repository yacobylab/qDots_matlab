function spect = scopespect(spect_par, tds)
% spect = scopespect(spect_par, tds)
%
% Note that this program was originally written for use with a DAQ and that
% the adatation for use with an oscilloscope is a bit of a hack.
%
% tds is an open instrument connected to the Tek scope. The samplerate and
% record length are set according to spectpar, all other settings are kept.
% spectpar is a struct with the following optional entries. If an entry is
% not given, a default value is used.
%   
% entry                 default
%--------------------------------
% hardFiltOrd           1
% hardFiltFreq          30000
% aves                  10
% filterOrd             15
% amp_gain              1000
% roughMaxFreq          1e4
%       if roughMaxFreq > sampRate/2;
%            roughMaxFreq = sampRate/2;
%       end    
% roughMinFreq          1
%       if roughMinFreq > roughMaxFreq
%            rougMinFreq =  roughMaxFreq/10;
%       end
% samplerate            ai.SampleRate
% channel               1
% numFFTs               1 number of frequency window
% FFTshift              8  factor by which to shift frequency window
%       FFTshift = 2 ^ floor(log2(FFTshift));    
% filterData            1
% FFTOverlaps           1
% plotFFT               1
% plotTime              1
% plotTimeWindow        15   duration of time trace to be shown. Inf for all
% plotTimeSamples       5000    number of samples displaued in time window.
%                               Inf = no downsamplng

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


%Initialize variables for the program state event loop.
%This loop cycles through to check the values of the various User Interface elements.
%These variables are global so they can be accessed by the ParseEvents function.

persistent hTimeFig;
persistent hFFTFig;

if isfield(spect_par, 'samplerate')
   fprintf(tds,'HOR:MAI:SAMPLER %f', spect_par.samplerate)
end

sampRate = query(tds, 'HOR:MAI:SAMPLER?', '%s\n', '%f');

fprintf(tds,'HOR:DEL:MOD 0');
fprintf(tds,'HOR:DEL:POS 0');
fprintf(tds,'ACQ:MODE HIR');
%Defaults for settable parameters
if isfield(spect_par, 'hardFiltOrd')
    hardFiltOrd = spect_par.hardFiltOrd;
else
    hardFiltOrd = 1;
end

if isfield(spect_par, 'hardFiltFreq')
    hardFiltFreq = spect_par.hardFiltFreq;
else
    hardFiltFreq = 30000;
end

if isfield(spect_par, 'aves')
    aves = spect_par.aves;
else
    aves = 10;
end


if isfield(spect_par, 'filterOrd')
    filterOrd = spect_par.filterOrd;
else
    filterOrd = 15;
end


if isfield(spect_par, 'amp_gain')
    amp_gain = spect_par.amp_gain;
else
    amp_gain = 1000;
end

if isfield(spect_par, 'channel')
    ch = spect_par.channel;
else
    ch = 1;
end

if isfield(spect_par, 'roughMaxFreq')
    roughMaxFreq = spect_par.roughMaxFreq;
else
    roughMaxFreq = 1e4;
end

if roughMaxFreq > sampRate/2;
   roughMaxFreq = sampRate/2;
end

if isfield(spect_par, 'roughMinFreq')
    roughMinFreq = spect_par.roughMinFreq;
else
    roughMinFreq = 1;
end

if roughMinFreq > roughMaxFreq
    rougMinFreq =  roughMaxFreq/10;
end


% if isfield(spect_par, 'numDecades')
%     numDecades = min(spect_par.numDecades, 100*floor(log10(sampRate * 5 / roughMaxFreq)));
% else
%     numDecades = 1;
% end

if isfield(spect_par, 'FFTshift')
    FFTshift = 2 ^ floor(log2(spect_par.FFTshift));
else
    FFTshift = 8;
end

if isfield(spect_par, 'numFFTs')
    numFFTs = min(spect_par.numFFTs, ceil((log2(sampRate / roughMaxFreq))/log2(FFTshift)));
else
    numFFTs = 1;
end

if isfield(spect_par, 'filterData')
    filterDataFlg = spect_par.filterData;
else
    filterDataFlg = 1;
end

if 0&&isfield(spect_par, 'FFTOverlaps')
    FFTOverlaps = spect_par.FFTOverlaps;
else
    FFTOverlaps = 1;
end

if isfield(spect_par, 'plotFFT')
    plotFFT = spect_par.plotFFT;
else
    plotFFT = 1;
end

if isfield(spect_par, 'plotTime')
    plotTime = spect_par.plotTime;
else
    plotTime = 1;
end

if isfield(spect_par, 'plotTimeWindow')
    plotTimeWindow = spect_par.plotTimeWindow;
else
    plotTimeWindow = 15;
end

if isfield(spect_par, 'plotTimeSamples')
    plotTimeSamples = spect_par.plotTimeSamples;
else
    plotTimeSamples = 5000;
end


%what sbout those?
plotDataMinMax = 0;
cutoffFreq = 0.45;

%FFT window
if (isempty(hFFTFig) || ~ishandle(hFFTFig)) && plotFFT
    hFFTFig = figure('position',[5 55 815 450]);
end

%Time trace window
if (isempty(hTimeFig) || ~ishandle(hTimeFig)) && plotTime
    hTimeFig = figure('position',[5 580 815 213]);
end

% figure out filter parameters for a filter that has a cutoff frequency of
% 0.5 the Nyquist frequency of whatever array it's filtering.
%[n,Wn] = buttord(0.45,0.5,0.1,80);
[b,a] = butter(filterOrd,cutoffFreq);
[bNone,aNone] = butter(filterOrd, 2*cutoffFreq);
%[b,a] = ellip(10,1,80,0.4);

 %cutSize = 2^16; % cutSize is the number of data points to grab each time.
 %dataCut = zeros(cutSize,1);

 spect.FFTdata = [];
 spect.FFTfreq = [];

 
 for i = 1:numFFTs

     downsamples = max(floor(log2(sampRate / roughMaxFreq)) - 1 - log2(FFTshift) * (i-1), 0);
     
     sampFreq = sampRate / 2 ^ downsamples;
     FFTMaxFreq = sampFreq / 2;  
     %sampTime = 1/sampFreq;
     
     % number of bins = FFTMaxFreq / FFTMinFreq 
     % number of samples = 2 * number of bins
     
     %log2(FFTMaxFreq / roughMinFreq)
     
     FFTsize = max(4, 2^(ceil(log2(FFTMaxFreq / roughMinFreq) - log2(FFTshift) * (i-1)) + 1));
     FFTMinFreq = sampFreq / FFTsize;
     FFTfreq = (2:FFTsize/2-1) * FFTMinFreq;
     FFTstep = FFTsize / FFTOverlaps;  % each FFT of the data is offset from the previous one by this amount
     
     cutSize = FFTsize * 2^downsamples;
               
     fprintf(tds,'HOR:RECO %d', cutSize);
     if query(tds,'HOR:RECO?', '%s\n', '%f') < cutSize
         fprintf(tds,'HOR:RECO %d', 2*cutSize);
     end

     if query(tds,'HOR:RECO?', '%s\n', '%f') < cutSize
        fprintf('Unable to set record length. Aborting.\n');
        return;
     end

     fprintf(tds, 'DAT:STAR 1');
     fprintf(tds, 'DAT:STOP %d', cutSize);
     
     %FFTsize
     %numCuts = ceil(FFTsize * aves * 2 ^ downsamples / cutSize); 
     totSamples = aves * cutSize / 2 ^ downsamples;
  
     
     windowf = (1-cos(linspace(0, 2*pi, FFTsize)))/2;

      
     workData = repmat(NaN, 1, totSamples);
     sampTimes = (0:totSamples-1) / sampFreq ;
     linkconditions = [];  
     FFTdata = zeros(1, FFTsize/2-2);
  
     % number of samples to be plotted in time window 
     plotSamples = round(plotTimeWindow * sampFreq);
     plotDownSamples = max(1, floor(plotSamples/plotTimeSamples));
     
% normalisation: / count / n * 2 / sqrt(1/(time * n)) (from ffts.m)
%	n = FFTsize
%   time = sampTime = 1/sampFreq = 2 * fftMaxFreq    
%  =  2 * sqrt(time / n) =  2 / sqrt(sampFreq n) = 2 / sqrt(minFreq n^2)

    % changed 4->16/3 on 04/20/09 - see modelling.txt
     hardFiltCorrect = (1 + (FFTfreq./hardFiltFreq).^2).^ hardFiltOrd * 16/3./ ...
         (amp_gain^2 * FFTMinFreq * FFTsize^2) ; 
     %correction for hardware filtering
    
    disp(['Window ', num2str(i)]);
    disp(['Max freq: ', num2str(FFTMaxFreq) ' Hz'])
    disp(['Min freq : ', num2str(FFTMinFreq) ' Hz'])
    disp(['FFT size: ', num2str(FFTsize)]);    
    disp(['Downsampling: ', num2str(downsamples)]);
    
    insPlace = 1;
    newPlace = FFTsize - FFTstep; % newPlace = first new data point. Initialise to allow 
    %one full FFTsize to come in before starting to process
    nOld = 0;
       
    nacq = 0;
    for loop = 0:aves % discard first block of data
        %pause(1);
        %dataCut = getdata(ai, cutSize);
        while  nacq >= query(tds,'ACQ:NUMAC?', '%s\n', '%d');
            pause(.02);
        end

        dataCut = tdsgetcurve(ch, tds);
        nacq = query(tds,'ACQ:NUMAC?', '%s\n', '%d');
        
        %tic
        %dataCut = randn(1, cutSize) * sqrt(sampRate/2);

        if filterDataFlg 
            if downsamples == 0   %if not downsampling: just filter
                if isempty(linkconditions)
                    [dataCut,linkconditions] = filter(bNone,aNone,dataCut);
                else
                    [dataCut,linkconditions] = filter(bNone,aNone,dataCut,linkconditions);
                end
            else
                if isempty(linkconditions)
                    for downLoop=1:downsamples
                        [dataCut, linkconditions(:,downLoop)] = filter(b,a,dataCut);
                        dataCut = downsample(dataCut,2);
                        %figure;plot(dataCut);drawnow;
                    end
                else
                    for downLoop=1:downsamples
                        [dataCut, linkconditions(:,downLoop)] = filter(b,a,dataCut,linkconditions(:,downLoop));
                        dataCut = downsample(dataCut,2);
                    end
                end
            end
        else    
            dataCut = downsample(dataCut, 2^downsamples);
        end
        %      dataCut = filtfilt (b,a,dataCut);
        
        if loop < 1  %only initialise initial conditions
            continue;
        end;
        
        downsampledSize = length(dataCut);
        workData(insPlace:insPlace+downsampledSize-1) = dataCut; 

        insPlace = insPlace + downsampledSize;
        nNew = 0;
        newFFTdata = zeros(1, FFTsize);
        %size(windowf)
        %size(workData((newPlace + FFTstep - FFTsize + 1):(newPlace + FFTstep)))
        %FFTsize
        while insPlace - newPlace >= FFTstep %&& insPlace > FFTsize % ready for new FFT
            newFFTdata = newFFTdata + ...
                abs(fft(windowf .* workData((newPlace + FFTstep - FFTsize + 1):(newPlace + FFTstep)))).^2;
        
            newPlace = newPlace + FFTstep;
            nNew = nNew + 1;
        end    
        newFFTdata = newFFTdata(3:end/2); % exclude DC bin
        newFFTdata = newFFTdata .* hardFiltCorrect;

        if nNew > 0
            FFTdata = FFTdata * (nOld / (nOld + nNew)) + newFFTdata ./ (nOld + nNew); 
        end;
        nOld = nOld + nNew;
        
        if plotFFT && nNew > 0
          if double(get(hFFTFig, 'CurrentCharacter')) == 27
            set(hFFTFig, 'CurrentCharacter',' ');
            return
          end
            figure (hFFTFig);
            hold off;
            loglog(FFTfreq, sqrt(FFTdata));  
        end

        if plotTime
    	  if double(get(hTimeFig, 'CurrentCharacter')) == 27
            set(hTimeFig, 'CurrentCharacter',' ')
            return
          end
            figure(hTimeFig);
            hold off;
            plot(downsample(sampTimes(max(1, insPlace - plotSamples):insPlace-1), plotDownSamples),...
                downsample(workData(max(1, insPlace - plotSamples):insPlace-1), plotDownSamples),'r');
            %size(downsample(workData(max(1, insPlace - plotSamples):insPlace-1), plotDownSamples))
        end
    
        drawnow;
    
        %disp([num2str(ai.SamplesAvailable), ' samples available']);
    end
    %continue;
    
    spect.timestep(i) = 1/sampFreq;
    if isfield(spect_par, 'savedata') && spect_par.savedata
        spect.data{i} = workData;
    end
    
       
    if i == 1
        first = 1;
    else 
        first = ceil(spect.FFTfreq(end) / FFTMinFreq); % max freq of last window
    end
    
    if i == numFFTs
        last = floor(FFTsize * cutoffFreq); 
    else
        last = floor(FFTsize * cutoffFreq * .9); % stay away from cutoff
    end 
    
    spect.FFTdata = [spect.FFTdata, sqrt(FFTdata(first:last))];
    spect.FFTfreq = [spect.FFTfreq, FFTfreq(first:last)];
   
end

figure (hFFTFig);      
loglog(spect.FFTfreq, spect.FFTdata);
% a = axis;
%axes('color', 'none', 'XTick', [], 'yaxislocation', 'right', ...
%    'yscale', 'log', 'tickdir', 'out', 'ylim', a(3:4) * sqrt(FFTMinFreq) * 2);


spect.aves = aves;
spect.cutoffFreq = cutoffFreq;
spect.filterOrd = filterOrd;
spect.hardwareSampFreq = sampRate;
spect.gain = amp_gain;
spect.hardFiltOrd = hardFiltOrd;
spect.hardFiltFreq = hardFiltFreq;
spect.spectpar = spect_par;

% problems if samples needed < one cut


% derivation of normalisaion:
% rms = sqrt(BW) * sqrt(S) 
% on other hand  <f f> = Int dw/(2 pi T) |f(w)|^2 = BW * S 
% for S = 2 * |f(w)|^2/T
% f(w) = int dt e^(iw t) f(t) = FFT * T / N
% thus S = 2 * |FFT|^2 * T /N^2 = 2 |FFT/N|^2 / F_min as done above.

% spectrum of white noise:
% correlator = delta(t - t') = <x^2> * delta => RMS = 1/sqrt(delta)
% FT = 1/ (2 pi sqrt(delta))

% < f(x) f(x- x') > = int 0 .. (delta - x')/delta dx a^2  = a^2 max (0, (1 - |x'|/delta)) 
% int = delta a^2 = prefactor of delta function.    
% thus limiting form for corellator is delta(t) * a^2 * delta
% for a = sqrt(1/delta) = sqrt(samplerate), get a FT of 1. Spectrum of 2.
% Thus need a = sqrt(samplerate/2).
% does not fit exactly pi/2 works better. (pretty well)


% improvement possibilities:
% downsample data before plotting
% reduce plotting width
% reasonable default: 1500 samples - more than on screen
% 15 s of data
