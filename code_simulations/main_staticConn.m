% Izikevich model with dynamic inference of connectivity
% Used for Fig.3 (panels C to F) and Fig.4 

clear all; clc;

% PATH setting
PATH = ['path_to_main_directory\results\static_connectivity\']; % path to directory where results have to be saved

% Parameters setting
% Attention: to reproduce the paper results (in particular Fig.5C) it is
% essential to run the script several times changing the randSeed variable from 1 to 5 (5 repetitions of the simulation)
randSeed = 1; % change this parameter to run different repetitions of the simulated network
rand('seed',randSeed);

% Network parameters
M=10;                 % number of synapses per neuron
D=20;                  % maximal conduction delay 
% excitatory neurons   % inhibitory neurons      % total number 
Ne=80;                Ni=20;                   N=Ne+Ni;

% Analysis parameters
doHO = 1; % compute higher-order quantities
doXCov = 1;
doXCorr = 1;
norm_method = 'normalized'; % normalization method for XCorr and XCov either 'none' or 'normalized'
minutes_range = [5 10 20 30 50 70 90 120 180]; % 5:15:95

max_delay = 50; % maximum delay to compute metrics

null_thresh = (99:-1:0);
% Initialize ratios
TEtp= zeros(numel(minutes_range),numel(null_thresh)); TEfp = TEtp; TEprec = TEtp; TErec = TEtp;
HOTEtp= TEtp; HOTEfp = TEtp; HOTEprec = TEtp; HOTErec = TEtp;
XCovtp= TEtp; XCovfp = TEtp; XCovprec = TEtp; XCovrec = TEtp;
XCorrtp= TEtp; XCorrfp = TEtp; XCorrprec = TEtp; XCorrrec = TEtp;

% Initialize ratios
TEtp= zeros(numel(minutes_range),numel(null_thresh)); TEfp = TEtp; TEprec = TEtp; TErec = TEtp;
HOTEtp= TEtp; HOTEfp = TEtp; HOTEprec = TEtp; HOTErec = TEtp;
XCovtp= TEtp; XCovfp = TEtp; XCovprec = TEtp; XCovrec = TEtp;
XCorrtp= TEtp; XCorrfp = TEtp; XCorrprec = TEtp; XCorrrec = TEtp;

for mIdx = 1:numel(minutes_range)
    minutes = minutes_range(mIdx);
    disp(['Analysis for ',num2str(minutes),' minutes'])

    % Simulate Izhikevich network for the selected parameters
    [spikeTrains,~,post,delays]=Izhikevich_network_v1(minutes,M,D,Ne,Ni);

    %% FC connectivity analysis

    % Ground truth connectivity matrix
    disp('Computing ground truth connectivity matrix')
    gtConn = zeros(N,N);
    gtDelay = nan(N,N);
    for sIdx = 1:N %sender
        for rIdx = 1:N %receiver
            tmpConn = sum(post(sIdx,:) == rIdx);
            if tmpConn
                gtConn(sIdx, rIdx) = tmpConn;
            end
        end
        for dIdx = 1:D
            tmpDel = delays(sIdx,dIdx);
            if ~isempty(tmpDel{1})
                for nIdx = 1:numel(tmpDel{1})
                    gtDelay(sIdx, post(sIdx,tmpDel{1}(nIdx))) = dIdx;
                end
            end
        end
    end        

    %timeWin = 1000*minutes*60;
    disp('Computing Transfer Entropy')
    asdf = SparseToASDF(spikeTrains, 1);

    [peakTE, CI, allTE, TEdelays] = ASDFTE(asdf, 1:max_delay, 1, 1); % using all the pairs in asdf, with delay of 1 to 3 bins
    peakTE = peakTE - diag(diag(peakTE));
    if doHO
        [peakHOTE, HOCI, allTE, HOTEdelays] = ASDFTE(asdf, 1:max_delay, 5, 5); % using all the pairs in asdf, with delay of 1 to 3 bins
        peakHOTE = peakHOTE - diag(diag(peakHOTE));
    end    

    % Compute cross covariance
    if doXCov || doXCorr
        disp('Computing cross correlation')

        peakXCov = zeros(N,N);
        peakXCorr = zeros(N,N);
        XCovDelays = zeros(N,N);
        XCorrDelays = zeros(N,N);

        for i = 1:N
            disp(['Computing correlation for emitter neuron ',num2str(i)])
            for j = 1:N
                if i~=j
                    if doXCov
                        [xCov,lags] = xcov(spikeTrains(i,:),spikeTrains(j,:),max_delay,norm_method);
                        xCov(ceil(numel(xCov)/2))=[]; % eliminating 0 delays
                        lags(ceil(numel(xCov)/2))=[]; % eliminating 0 delays
                        negIdxs = (lags<0);  % computing xcov(a,b) shifts b by lag, so if we find maximum at lag=(-del) it means that the link is a->b with lag = del
                        % Neurons can be either correlated or anticorrelated --> take abs value
                        [~,tmpxCovDelay]=max(abs(xCov(negIdxs)));
                        XCovDelays(i,j) = -lags(tmpxCovDelay);
                        peakXCov(i,j)=xCov(tmpxCovDelay);
                    end
                    if doXCorr
                        [xCorr,lags] = xcorr(spikeTrains(i,:),spikeTrains(j,:),max_delay,norm_method);
                        xCorr(ceil(numel(xCorr)/2))=[]; % eliminating 0 delays
                        lags(ceil(numel(xCorr)/2))=[]; % eliminating 0 delays
                        negIdxs = (lags<0);  % computing xcov(a,b) shifts b by lag, so if we find maximum at lag=(-del) it means that the link is a->b with lag = del
                        % Neurons can be either correlated or anticorrelated --> take abs value
                        [~,tmpxCorrDelay]=max(abs(xCorr(negIdxs)));
                        XCorrDelays(i,j) = -lags(tmpxCorrDelay);
                        peakXCorr(i,j)=xCorr(tmpxCorrDelay);
                    end
                end
            end
        end
        peakXCov = peakXCov - diag(diag(peakXCov));
        peakXCorr = peakXCorr - diag(diag(peakXCorr));
    end
    
    %% Compute precision and recall

    % Compute TP and FP
    null_thresh = (99:-1:0);
    clear TEtp CItp TEfp CIfp
    idxT = (gtConn == 1); %idxs of true connections
    idxF = (gtConn == 0); %idxs of false connections

    for idx = 1:numel(null_thresh)
        TEthr = prctile(peakTE(:),null_thresh(idx));
        tmpTEmap = (peakTE > TEthr);
        CIthr = prctile(CI(:),null_thresh(idx));
        tmpCImap = (CI > CIthr);
        if doXCov
            XCovthr = prctile(abs(peakXCov(:)),null_thresh(idx)); % XC = cross correlation
            tmpXCovmap = (abs(peakXCov) > XCovthr);
        end
        if doXCorr
            XCorrthr = prctile(abs(peakXCorr(:)),null_thresh(idx)); % XC = cross correlation
            tmpXCorrmap = (abs(peakXCorr) > XCorrthr);
        end
        if doHO
            HOTEthr = prctile(peakHOTE(:),null_thresh(idx));
            tmpHOTEmap = (peakHOTE > HOTEthr);
            HOCIthr = prctile(HOCI(:),null_thresh(idx));
            tmpHOCImap = (HOCI > HOCIthr);
        end

        % True positives computed as the average of inferred connectivity = 1
        % over gd present and absent connections respectively
        TEtp(idx) = mean(tmpTEmap(idxT)); % TPR
        TEfp(idx) = mean(tmpTEmap(idxF)); % FPR
        TEprec(idx) = sum(sum(tmpTEmap(idxT)))/sum(sum(tmpTEmap)); % precision = TP/(TP+FP) = TP/(inferred P)
        TErec(idx) = sum(sum(tmpTEmap(idxT)))/sum(sum(gtConn)); % recall = TP/(TP+FN) = TP/(Actual P)

        if doXCov
            XCovtp(idx) = mean(tmpXCovmap(idxT));
            XCovfp(idx) = mean(tmpXCovmap(idxF));
            XCovprec(idx) = sum(sum(tmpXCovmap(idxT)))/sum(sum(tmpXCovmap)); % precision = TP/(TP+FP) = TP/(inferred P)
            XCovrec(idx) = sum(sum(tmpXCovmap(idxT)))/sum(sum(gtConn)); % recall = TP/(TP+FN) = TP/(Actual P)
        end
        if doXCorr
            XCorrtp(idx) = mean(tmpXCorrmap(idxT));
            XCorrfp(idx) = mean(tmpXCorrmap(idxF));
            XCorrprec(idx) = sum(sum(tmpXCorrmap(idxT)))/sum(sum(tmpXCorrmap)); % precision = TP/(TP+FP) = TP/(inferred P)
            XCorrrec(idx) = sum(sum(tmpXCorrmap(idxT)))/sum(sum(gtConn)); % recall = TP/(TP+FN) = TP/(Actual P)
        end
        if doHO
            HOTEtp(idx) = mean(tmpHOTEmap(idxT));
            HOTEfp(idx) = mean(tmpHOTEmap(idxF));
            HOTEprec(idx) = sum(sum(tmpHOTEmap(idxT)))/sum(sum(tmpHOTEmap));
            HOTErec(idx) = sum(sum(tmpHOTEmap(idxT)))/sum(sum(gtConn));
        end
    end

    %% Compute correlation between real and inferred delays on pairs having a GT link
    
    % TE delays
    tmpCorr = corrcoef(TEdelays(~isnan(gtDelay(:))),gtDelay(~isnan(gtDelay(:))));
    delayCorr.TE = tmpCorr(1,2);
    delayErr.TE = mean(TEdelays(~isnan(gtDelay(:))) - gtDelay(~isnan(gtDelay(:))));
    delayMatch.TE = mean(TEdelays(~isnan(gtDelay(:))) == gtDelay(~isnan(gtDelay(:))));

    if doHO
        % HOTE delays
        tmpCorr = corrcoef(HOTEdelays(~isnan(gtDelay(:))),gtDelay(~isnan(gtDelay(:))));
        delayCorr.HOTE = tmpCorr(1,2);
        delayErr.HOTE = mean(HOTEdelays(~isnan(gtDelay(:))) - gtDelay(~isnan(gtDelay(:))));
        delayMatch.HOTE = mean(HOTEdelays(~isnan(gtDelay(:))) == gtDelay(~isnan(gtDelay(:))));
    end
    if doXCov
        % XCov delays
        tmpCorr = corrcoef(XCovDelays(~isnan(gtDelay(:))),gtDelay(~isnan(gtDelay(:))));
        delayCorr.XCov = tmpCorr(1,2);
        delayErr.XCov = mean(XCovDelays(~isnan(gtDelay(:))) - gtDelay(~isnan(gtDelay(:))));
        delayMatch.XCov = mean(XCovDelays(~isnan(gtDelay(:))) == gtDelay(~isnan(gtDelay(:))));
    end
    if doXCorr
        % XCorr delays
        tmpCorr = corrcoef(XCorrDelays(~isnan(gtDelay(:))),gtDelay(~isnan(gtDelay(:))));
        delayCorr.XCorr = tmpCorr(1,2);
        delayErr.XCorr = mean(XCorrDelays(~isnan(gtDelay(:))) - gtDelay(~isnan(gtDelay(:))));
        delayMatch.XCorr = mean(XCorrDelays(~isnan(gtDelay(:))) == gtDelay(~isnan(gtDelay(:))));
    end
    fname = ['staticConnInference_rep',num2str(randSeed),'_min',num2str(minutes),'_',date];
    save([PATH,fname],'-v7.3')
end

