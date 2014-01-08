function varargout = smpSnglShot(varargin)
% varargout = smpSnglShot(varargin)
% varargin: new raw data (one for each channel), previously computed data, analysis definition.
% varargout: newly computed data added to older.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


params = varargin{end};

nread = size(params.readout, 2);
nchan = length(varargin)-1-length(params.datadef); 

if nchan ~= length(nargout) % could allow pass-through of new data.
    error('Number of inputs and outputs inconsistent.');
end

for i = 1:nchan
    varargin{i} = reshape(varargin{i}, params.period, length(varargin{i})/params.period); % could trim timetrace
end

data = [varargin{1:nchan}]; % catenate channels


flags = [any(strcmp('mean', {params.datadef.type})), any(strcmp('hist', {params.datadef.type})), ...
    any(strcmp('bin', {params.datadef.type})), any(strcmp('mom', {params.datadef.type})), ...
    any(strcmp('corr', {params.datadef.type})), any(strcmp('bmom', {params.datadef.type})), 
    any(strcmp('bcorr', {params.datadef.type}))];
    
if any(flags)    
    if nread > 0
        mdata = zeros(size(data, 1), nchan * nread);
        for j = 0:nchan-1
            for i = 1:nread
                mdata(:, i+j*nread) = ...
                    mean(data(:, params.readout(1, i) + (0:params.readout(2, i)-1) + j * params.period), 2);
            end
        end
    else
        mdata = data;
    end
end

if any(flags([3, 6, 7])) 
    bdata = false(size(mdata));
    if ndims(thresh) == 1 % same threshold for each channel
        for j = 0:nchan-1
            bdata(j * nread + (1:nread)) = mdata(j * nread + (1:nread)) > params.thresh(j+1);
        end
    else %different thresholds given for each readout.
        bdata  = mdata > repmat(params.thresh(:)', size(data, 1), 1);
    end
end

for i = 1:length(params.datadef)    
    varargout{i} = varargin{i + nchan};
    ind = sum(isfinite(varargin{1+nchan}(:, 1, 1)))+1;
    
    switch params.dataadef(i).type
        case 'ave' % average of all traces (not possible/meaningful with driver level averaging)
             varargout{i}(ind, :) = mean(data);
             
        case 'mean' % per pulse average 
             varargout{i}(ind, :, :) = mdata;
             
        case 'hist' % histogram also good for probabilities
            hl = datadef(i).par{1};

            if size(hl, 2) == nchan
                for j = 1:nchan
                    ind2 = (1:nread) * (j-1) *nread;
                    varargout{i}(ind, :, ind2) = histc(mdata(:, ind2), hl(:, j));
                end
            else
                for j = 1:nchan*nread
                    varargout{i}(ind, :, j) = histc(mdata(:, j), hl(:, j));
                end
            end
            
        case 'bin' % binary - saves memory, useful for correlations.
                %derived from ave, param(1) = threshold
             varargout{i}(ind, :, :) = bdata;
    
        case 'mom' % moments of distribution
            mdef = datadef(i).par{1};
            if length(datadef(i).par) == 1
                lag = zeros(size(mdef));
            else                
                lag = repmat(datadef(i).par{2}, 1, size(mdef, 2)./size(lag, 2));
                %use same lags for all moments computed if not specified individually                
            end
            
            maxlag = max(lag(:));
            tmp = ones(size(data, 1)-maxlag, size(mdef, 1));
               
            for k = 0:maxlag             
                for j = 1:size(mdef, 1)
                    mask = mdef(j, :) > 0 & lag(j, :) == k;
                    tmp(mask , :) = tmp(mask , :) .* mdata(1+k:end-maxlag+k, mdef(j, mask));
                end
            end

            varargout{i}(ind, :) = mean(tmp);
            
%             if length(datadef(i).par) == 1 || size(datadef(i).par{2}, 2) == 1
%                 % same or no offset
%                 for j = 1:size(mdef, 1)
%                     mask = mdef(j, :) > 0;
%                     tmp(mask , :) = tmp(mask , :) .* mdata(:, mdef(j, mask));
%                 end
%             else
%                 for j = 1:size(mdef, 1)
%                     for
%                     mask = mdef(j, :) > 0;
%                     tmp(mask , :) = tmp(mask , :) .* mdata(:, mdef(j, mask));
%                 end
%             end

                            
        case 'corr'  %correlation functions, based on ave or bin, specify channels, pulses, offest

        case 'bmom' % moments of binary distribution
                
        case 'bcorr'  %binary correlation functions, based on ave or bin, specify channels, pulses, offest

    end
    
end
