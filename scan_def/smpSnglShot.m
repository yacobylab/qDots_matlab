function varargout = smpSnglShot(varargin)
% varargout = smpSnglShot(varargin)
% varargin: new raw data (one for each channel) + analysis definition.
% varargout: analysis results.
% The last element of varargin is a struct defining the processing to be done.
% Its fields are: 
%    datadef: struct array specifying data to be computed, has fields
%       type: (g)ave, mean, (g)hist, bin, (b)mom, corr, (b)corr
%               prefix 'g' averages/sums over all pulses, 'b' uses binary data
%       par: Cell array with parameters
%           mean: mean value of each measurement.
%                 par{1} is an optional list of measurements to be recorded (default all)
%           bin: binary value of each measurement.
%                 par{1} is an optional list of measurements to be recorded (default all)
%           hist: par{1} = histogram limits, either one per channel or one per measurement
%                 par{2} is an optional list of measurements to be analyzed (default all)
%           mom: par{1} = index matrix to readouts, each column specifying a moment
%                    An empty vector specifies the mean of all measurements and channels.                        
%                par{2} = time lag, matrix or vector for all moments computed (defaults to 0)
%           corr: par{1} = 2 x n matrix specifing pairs to be correlated
%                    an empty vector specifies autocorrelations of all measuremens and channels.
%                 par{2} = max lag.
%    readout: 2 x n array specyfing first sample and # of samples to be averaged
%             for each measurement. [] = no averaging, e.g. for using driver level averaging.
%
%    n: [#samples per measurement group, #reps per group, #groups, [#outer reps]]
%    
%    thresh: discrimination threshold, either one per channel or one per measurement.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


params = varargin{end};
n = params.n;

nread = size(params.readout, 2);
nchan = length(varargin)-1; 

if length(params.datadef) ~= nargout 
    error('Number of inputs and outputs inconsistent.');
end

% several pulses in a row: 
for i = 1:nchan
    if length(n) == 3
        varargin{i} = reshape(permute(reshape(varargin{i}, n(1), n(2), n(3)), [2 1 3]), [n(2), n(1)*n(3)]);
        n(4) = 1;
    else % whole sequence repeated
        varargin{i} = reshape(permute(reshape(varargin{i}, n), [2 4 1 3]), [n(2)*n(4), n(1)*n(3)]);
    end
    % could trim timetrace
end
% data order: inner repetitions * outer repetitions, channels * pulses * groups

% assume same measurements for each group if ?.
% if 0
%     datadef = repmat(params.datadef, 1, params.n(3));
% elseif  length(params.pardef) == n(2)* n(3);
%     datadef = params.datadef;
% end

data = [varargin{1:nchan}]; % catenate channels
% reorganize sequential pulse switching.
%reshape(permute(reshape(data, [params.nrep, nperiods/params.nrep, params.period]), [1 3 2]),...
%    [params.nrep, nperiods*params.period/params.nrep]);

varargin = []; 

flags = [any(strcmp('ave', {params.datadef.type})) || any(strcmp('gave', {params.datadef.type})), ...
    any(strcmp('mean', {params.datadef.type})), ...
    any(strcmp('hist', {params.datadef.type})) || any(strcmp('ghist', {params.datadef.type})), ...
    any(strcmp('bin', {params.datadef.type})), ...
    any(strcmp('mom', {params.datadef.type})), any(strcmp('corr', {params.datadef.type})), ...
    any(strcmp('bmom', {params.datadef.type})), any(strcmp('bcorr', {params.datadef.type}))];
    
if any(flags(2:8))    
    if nread > 0 
        mdata = zeros(n(2)*n(4), nchan * nread);
        for j = 0:nchan-1
            for i = 1:nread
                mdata(:, i+j*nread) = ...
                    mean(data(:, params.readout(1, i) + (1:params.readout(2, i)) ...
                    + j*params.n(1)*params.n(3)), 2);
            end
        end
    else % driver level averaging.
        mdata = data;
        nread = n(1) * n(3);
    end
end

if ~flags(1) 
    data = []; % free memory
end

if any(flags([4, 7, 8])) 
    bdata = false(size(mdata));
    if numel(params.thresh) == nchan % single threshold per channel
        for j = 0:nchan-1
            bdata(:, j * nread + (1:nread)) = mdata(:, j * nread + (1:nread)) > params.thresh(j+1);
        end
    else %different thresholds given for each readout.
        bdata  = mdata > repmat(params.thresh(:)', n(2)*n(4), 1);
    end
end

for i = 1:length(params.datadef)    
    
    switch params.datadef(i).type
        case 'ave'  % average of all traces (not possible/meaningful with driver level averaging)
            if isempty(params.datadef(i).par) 
                varargout{i} = reshape(mean(data), [n(1), n(3) * nchan])';
            else
                data = reshape(data, [n(2) * n(4), n(1), n(3) * nchan]);
                varargout{i} = squeeze(mean(data(:, params.datadef(i).par{1}, :)))';
            end
            data = [];            
                
        case 'gave' % also average over measurements.
            if isempty(params.datadef(i).par)
                varargout{i} = squeeze(mean(reshape(mean(data), [n(1), n(3), nchan]), 2))';
            else
                data = reshape(data, [n(2) * n(4), n(1), n(3), nchan]);
                varargout{i} = squeeze(mean(mean(data(:, params.datadef(i).par{1}, :), 1), 3))';                
            end
            data = [];

        case 'mean' % per pulse average 
            if length(params.datadef(i).par) >= 1
                varargout{i}(:, :) = mdata(:, params.datadef(i).par{1})';
            else
                varargout{i}(:, :) = mdata';
            end
             
        case 'hist' % histogram also good for probabilities
            hl = params.datadef(i).par{1};

            if length(params.datadef(i).par) >= 2
                mask = params.datadef(i).par{2}';
            else
                mask = 1:nchan*nread;
            end
            varargout{i} = zeros(nread*nchan, size(hl, 2));

            if size(hl, 1) == nchan % one set of limits per channel
                for j = 1:nchan
                    ind2 = mask(mask > (j-1)*nread & mask <= j*nread);
                    varargout{i}(ind2, :) = histc(mdata(:, ind2), hl(j, :))';
                end
            else % different limits for each measurement.
                for j = mask
                    varargout{i}(j, :) = histc(mdata(:, j), hl(j, :))';
                end
            end
            varargout{i}(:, end) = [];
            
            % alternative to separate case below
            %if params.datadef(i).type(1) == 'g'
            %    varargout{i} = sum(varargout{i}, 1);
            %end

        case 'ghist'
            hl = params.datadef(i).par{1};

            if length(params.datadef(i).par) >= 2
                mask = params.datadef(i).par{2}';
            else
                mask = 1:nchan*nread; % not particularly useful unless channels use same hist lims.
            end
            varargout{i} = zeros(nchan, size(hl, 2));

            %ind2 = mask(mask > (j-1)*nread & mask <= j*nread);
            varargout{i} = histc(reshape(mdata(:, mask), n(2)*n(4)*length(mask), 1), hl(j, :))';
            varargout{i}(:, end) = [];

        case 'bin' % binary - saves memory, useful for correlations.
                %derived from ave, param(1) = threshold
                if length(params.datadef(i).par) >= 1
                    varargout{i}(:, :) = bdata(:, params.datadef(i).par{1})';
                else
                    varargout{i}(:, :) = bdata';
                end
                
        case {'mom', 'bmom'} % moments of distribution
            mdef = params.datadef(i).par{1};
            if isempty(mdef)
                mdef  = 1:nread*nchan;
            end
            %  Not implemented (yet): A column with negative integer specifies all measurements for a single channel.
            %if mdef < 0 % all measurements of a channel. Could extend interpretation
            %    mdef = (1:nread) + nread * (abs(mdef)-1);
            %end
            if length(params.datadef(i).par) == 1
                lag = zeros(size(mdef));
            else                
                lag = repmat(params.datadef(i).par{2}, 1, size(mdef, 2)./size(params.datadef(i).par{2}, 2));
                %use same lags for all moments computed if not specified individually                
            end
            
            maxlag = max(lag(:));
            tmp = ones(n(2)*n(4)-maxlag, size(mdef, 2));
               
            for k = 0:maxlag             
                for j = 1:size(mdef, 1)
                    mask = mdef(j, :) > 0 & lag(j, :) == k;
                    if params.datadef(i).type(1) == 'b'
                        tmp(:, mask) = tmp(:, mask) & bdata(1+k:end-maxlag+k, mdef(j, mask));
                    else
                        tmp(:, mask) = tmp(:, mask) .* mdata(1+k:end-maxlag+k, mdef(j, mask));
                    end
                end
            end

            varargout{i} = mean(tmp);
            
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
                            
        case {'corr', 'bcorr'}  %correlation functions
            mdef = params.datadef(i).par{1};            
            if isempty(mdef)
                mdef  = repmat(1:nread*nchan, 2, 1);
            end

            varargout{i} = zeros(size(mdef, 2), 2*params.datadef(i).par{2}+1);
            
            for j = 1:size(mdef, 2);
                if params.datadef(i).type(1) == 'b'
                    varargout{i}(j, :) = xcorr(bdata(:, mdef(1, j) ), bdata(:, mdef(2, j)), ...
                        params.datadef(i).par{2}, 'unbiased')';
                else
                    varargout{i}(j, :) = xcorr(mdata(:, mdef(1, j) ), mdata(:, mdef(2, j)), ...
                        params.datadef(i).par{2}, 'unbiased')';
                end
            end
    end
end
