function [staticConnMeasures, staticDelaysMeasures] = computeStaticConn_from_SpikeTrains_v1(spikeTrains,params)

staticConnMeasures = [];
staticDelaysMeasures = [];

asdf = SparseToASDF(spikeTrains, 1);

if params.doTE
    disp('Computing Transfer Entropy')
    [staticConnMeasures.peakTE, ~, ~, staticDelaysMeasures.TEdelays] = ASDFTE(asdf, 1:params.max_delay, 1, 1); % using all the pairs in asdf, with delay of 1 to 3 bins
    % Remove diagonal elements
    staticConnMeasures.peakTE = staticConnMeasures.peakTE - diag(diag(staticConnMeasures.peakTE));
end
if params.doHO
    disp('Computing Higher Order Transfer Entropy')
    [staticConnMeasures.peakHOTE, ~, ~, staticDelaysMeasures.HOTEdelays] = ASDFTE(asdf, 1:params.max_delay, params.HOTE_l, params.HOTE_k); % using all the pairs in asdf, with delay of 1 to 3 bins
    staticConnMeasures.peakHOTE = staticConnMeasures.peakHOTE - diag(diag(staticConnMeasures.peakHOTE));
end

% Compute cross covariance
if params.doXCov || params.doXCorr
    disp('Computing cross correlation')
    % xCorr = xcorr(spikeTrains',params.max_delay);
    % xCorr = reshape(xCorr,2*params.max_delay+1,N,N); % elems ordered as [11,21,31...;12,22,32,...;13,23,33,...]
    % % Neurons can be either correlated or anticorrelated --> take abs value
    % [~,xCorrDelay]=max(abs(xCorr),[],1);
    % xCorrDelay=squeeze(xCorrDelay);

    staticConnMeasures.peakXCov = zeros(params.N,params.N);
    staticConnMeasures.peakXCorr = zeros(params.N,params.N);
    staticDelaysMeasures.XCovDelays = zeros(params.N,params.N);
    staticDelaysMeasures.XCorrDelays = zeros(params.N,params.N);
    if params.doSLjitter
        staticDelaysMeasures.XCorrDelaysSh = zeros(params.nShuff,params.N,params.N);
        staticDelaysMeasures.XCovDelaysSh = zeros(params.nShuff,params.N,params.N);
        staticConnMeasures.peakXCorrSh = zeros(params.nShuff,params.N,params.N);
        staticConnMeasures.peakXCovSh = zeros(params.nShuff,params.N,params.N);
    end
    
    for i = 1:params.N
        disp(['Computing corr for neuron ',num2str(i)])
        for j = 1:params.N
            if i~=j % Actually Izhikevich allows this type of connections!
                if params.doXCov
                    %tmpSpikeTrains = [spikeTrains(i,:);spikeTrains(j,:)]';
                    [xCov,lags] = xcov(spikeTrains(i,:),spikeTrains(j,:),params.max_delay,params.norm_method);
                    xCov(ceil(numel(xCov)/2))=[]; % eliminating 0 delays
                    lags(ceil(numel(xCov)/2))=[]; % eliminating 0 delays
                    negIdxs = (lags<0);  % computing xcov(a,b) shifts b by lag, so if we find maximum at lag=(-del) it means that the link is a->b with lag = del
                    %xCov = reshape(xCov,2*params.max_delay+1,2,2); % elems ordered as [11,21;12,22]
                    % Neurons can be either correlated or anticorrelated --> take abs value
                    [~,tmpxCovDelay]=max(abs(xCov(negIdxs)));
                    staticDelaysMeasures.XCovDelays(i,j) = -lags(tmpxCovDelay);
                    staticConnMeasures.peakXCov(i,j)=xCov(tmpxCovDelay);
                end
                if params.doXCorr
                    %tmpSpikeTrains = [spikeTrains(i,:);spikeTrains(j,:)]';
                    [xCorr,lags] = xcorr(spikeTrains(i,:),spikeTrains(j,:),params.max_delay,params.norm_method);
                    xCorr(ceil(numel(xCorr)/2))=[]; % eliminating 0 delays
                    lags(ceil(numel(xCorr)/2))=[]; % eliminating 0 delays
                    negIdxs = (lags<0);  % computing xcov(a,b) shifts b by lag, so if we find maximum at lag=(-del) it means that the link is a->b with lag = del
                    %xCov = reshape(xCov,2*params.max_delay+1,2,2); % elems ordered as [11,21;12,22]
                    % Neurons can be either correlated or anticorrelated --> take abs value
                    [~,tmpxCorrDelay]=max(abs(xCorr(negIdxs)));
                    staticDelaysMeasures.XCorrDelays(i,j) = -lags(tmpxCorrDelay);
                    staticConnMeasures.peakXCorr(i,j)=xCorr(tmpxCorrDelay);
                end

                % Compute Jitter null for XCorr and XCov
                if params.doSLjitter
                    spikeTrain_i = spikeTrains(i,:);
                    spikeTrain_j = spikeTrains(j,:);
                    % Jitter the receiver neuron
                    parfor shuffle = 1:params.nShuff
                        tmp_j_shuffle = zeros(length(spikeTrain_j),1);
                        tmp_j_spikes = find(spikeTrain_j)+randi([-params.max_delay params.max_delay],1,length(find(spikeTrain_j)));
                        tmp_j_spikes = tmp_j_spikes(tmp_j_spikes<length(spikeTrain_j) & tmp_j_spikes>0);
                        tmp_j_shuffle(tmp_j_spikes) = 1;
                        
                        % Compute XCorr for jittered receiver
                        if params.doXCorr
                            [tmp_CC,lags] = xcorr(spikeTrain_i,tmp_j_shuffle,params.max_delay,params.norm_method);
                            tmp_CC(ceil(numel(tmp_CC)/2))=[]; % eliminating 0 delays
                            lags(ceil(numel(tmp_CC)/2))=[]; % eliminating 0 delays
                            negIdxs = (lags<0);  % computing xcov(a,b) shifts b by lag, so if we find maximum at lag=(-del) it means that the link is a->b with lag = del
                            %xCov = reshape(xCov,2*params.max_delay+1,2,2); % elems ordered as [11,21;12,22]
                            % Neurons can be either correlated or anticorrelated --> take abs value
                            [~,tmpxCorrDelay]=max(abs(tmp_CC(negIdxs)));
                            parXCorrDelaysSh(shuffle) = -lags(tmpxCorrDelay);
                            parXCorrPeakSh(shuffle) =tmp_CC(tmpxCorrDelay);
                        end
                        % Compute XCov for jittered receiver
                        if params.doXCov
                            [tmp_CC,lags] = xcov(spikeTrain_i,tmp_j_shuffle,params.max_delay,params.norm_method);
                            tmp_CC(ceil(numel(tmp_CC)/2))=[]; % eliminating 0 delays
                            lags(ceil(numel(tmp_CC)/2))=[]; % eliminating 0 delays
                            negIdxs = (lags<0);  % computing xcov(a,b) shifts b by lag, so if we find maximum at lag=(-del) it means that the link is a->b with lag = del
                            %xCov = reshape(xCov,2*params.max_delay+1,2,2); % elems ordered as [11,21;12,22]
                            % Neurons can be either correlated or anticorrelated --> take abs value
                            [~,tmpxCovDelay]=max(abs(tmp_CC(negIdxs)));
                            parXCovDelaysSh(shuffle) = -lags(tmpxCovDelay);
                            parXCovPeakSh(shuffle)=tmp_CC(tmpxCovDelay);
                        end
                    end
                    staticConnMeasures.peakXCorrSh(:,i,j) = parXCorrPeakSh;
                    staticConnMeasures.peakXCovSh(:,i,j) = parXCovPeakSh;
                    staticDelaysMeasures.XCorrDelaysSh(:,i,j) = parXCorrDelaysSh;
                    staticDelaysMeasures.XCovDelaysSh(:,i,j) = parXCovDelaysSh;
                end
            end
        end
    end

    % Remove diagonal elements
    if params.doSLjitter
        staticConnMeasures.peakXCorrSh(shuffle,:,:) = staticConnMeasures.peakXCorrSh(shuffle,:,:) - diag(diag(staticConnMeasures.peakXCorrSh(shuffle,:,:)));
        staticConnMeasures.peakXCovSh(shuffle,:,:) = staticConnMeasures.peakXCovSh(shuffle,:,:) - diag(diag(staticConnMeasures.peakXCovSh(shuffle,:,:)));
    end
    staticConnMeasures.peakXCov = staticConnMeasures.peakXCov - diag(diag(staticConnMeasures.peakXCov));
    staticConnMeasures.peakXCorr = staticConnMeasures.peakXCorr - diag(diag(staticConnMeasures.peakXCorr));
end

end