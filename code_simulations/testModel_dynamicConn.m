% Izikevich model with dynamic inference of connectivity
%
% v15 (for BI resubmission): v14 + heterogeneity of firing rates

clear all; clc

randSeed = 1;
delay_consistency = 1; % parameters controlling whether we want to estimate DFC at the delay inferred from static FC or not
rand('seed',randSeed);

% Parameters setting

% Simulation parameters
M=10;                 % number of synapses per neuron
D=20;                  % maximal conduction delay 
% excitatory neurons   % inhibitory neurons      % total number 
Ne=80;                Ni=20;                   N=Ne+Ni;

% Connectivity analysis parameters
params.doTE = 1;
params.doHO = 1; % compute higher-order Transfer Entropy
params.HOTE_l = 5; % dimensionality of the past of the emitter for HOTE
params.HOTE_k = 5; % dimensionality of the past of the receiver for HOTE
params.doXCov = 1; % compute cross-covariance
params.doXCorr = 0; % compute cross-correlation
params.doSLjitter = 0; % do Jitter null hypothesis
params.nShuff = 100; % number of null hypotheses surrogates
params.N = N;

minutes = 180; 
% Parameters for the moving window TE analysis
new_bin_minutes_range = [1,2,5,10,20,30]; % [1,2,5,10,20,30];
time_jump_minutes = 1;

params.max_delay = 50; % maximum delay to compute metrics
delay_err = 0;
sel_thresh = 95; % quantile threshold to build the connected subnetwork for DFC
params.norm_method = 'normalized'; % either 'none' or 'normalized'

[spikeTrains,s_time,post,delays]=Izhikevich_network_v1(minutes,M,D,Ne,Ni);

figure()
imagesc(spikeTrains);
xlabel('time [ms]')
ylabel('neurons')
title(['rasterplot of 3s of activity'])
axis([50000 50300 0 N]);

figure()
hold on
plot(squeeze(s_time(1,1,:)),'linewidth',2)
plot(squeeze(s_time(10,2,:)),'linewidth',2)
plot(squeeze(s_time(Ne+1,1,:)),'linewidth',2)
xlabel('time [s]')
ylabel('syn. weight')
legend('exc 1','exc 2','inh 1')
title('Example synaptic weights over 5 minutes')

figure()
hold on
histogram(mean(spikeTrains(1:Ne,:),2)*1000,10);
histogram(mean(spikeTrains(Ne+1:Ne+Ni,:),2)*1000,10);
legend('exc.','inh.')
xlabel('time [ms]')
ylabel('neurons')
title(['Firing rate distributions'])

%% Compute and plot average autocorrelation of synaptic weights
autCorr = zeros(size(s_time,1)*size(s_time,2),size(s_time,3)*2-1);
for i = 1:N
    for j = 1:M
        idx = (i-1)*M+j;
        autCorr(idx,:) = xcorr(squeeze(s_time(i,j,:)),'unbiased');
    end
end
figure()
x_vals = -(size(s_time,3)-1):(size(s_time,3)-1);
maxIdx = find(mean(autCorr,1)==max(mean(autCorr,1)));
% [minValue,closestIndex] = min(abs(mean(autCorr,1)-max(mean(autCorr,1))/2));
plot(x_vals,mean(autCorr,1),'linewidth',2)
% xline(closestIndex-maxIdx,'linewidth',2)
% text((3)*(closestIndex-maxIdx),3000,[num2str(closestIndex-maxIdx),'s'])
xlim([-(size(s_time,3)-1)/2,(size(s_time,3)-1)/2])
title('mean synaptic weight autocorrelation')

%% Info analysis

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

figure()
imagesc(gtDelay)
colorbar()
%gdConnRed = gdConn(neuIdx,neuIdx); % reduced conn. matrix
%
[staticConnMeasures, staticDelaysMeasures] = computeStaticConn_from_SpikeTrains_v1(spikeTrains,params);
%timeWin = 1000*minutes*60;

%% Evaluate TE precision

% Compute TP and FP
thresh = (99:-1:0);
clear TEtp CItp TEfp CIfp
idxT = (gtConn == 1); %idxs of true connections
idxF = (gtConn == 0); %idxs of false connections

for idx = 1:numel(thresh)
    if params.doTE
        TEthr = prctile(staticConnMeasures.peakTE(:),thresh(idx));
        tmpTEmap = (staticConnMeasures.peakTE > TEthr);
    end
    if params.doXCov
        XCovthr = prctile(abs(staticConnMeasures.peakXCov(:)),thresh(idx)); % XC = cross correlation
        tmpXCovmap = (abs(staticConnMeasures.peakXCov) > XCovthr);
    end
    if params.doXCorr
        XCorrthr = prctile(abs(staticConnMeasures.peakXCorr(:)),thresh(idx)); % XC = cross correlation
        tmpXCorrmap = (abs(staticConnMeasures.peakXCorr) > XCorrthr);
    end
    if params.doHO
        HOTEthr = prctile(staticConnMeasures.peakHOTE(:),thresh(idx));
        tmpHOTEmap = (staticConnMeasures.peakHOTE > HOTEthr);
    end

    % True positives computed as the average of inferred connectivity = 1
    % over gd present and absent connections respectively
    if params.doTE
        TEtp(idx) = mean(tmpTEmap(idxT)); % TPR
        TEfp(idx) = mean(tmpTEmap(idxF)); % FPR
        TEprec(idx) = sum(sum(tmpTEmap(idxT)))/sum(sum(tmpTEmap)); % precision = TP/(TP+FP) = TP/(inferred P)
        TErec(idx) = sum(sum(tmpTEmap(idxT)))/sum(sum(gtConn)); % recall = TP/(TP+FN) = TP/(Actual P)
    end
    
    if params.doXCov
        XCovtp(idx) = mean(tmpXCovmap(idxT));
        XCovfp(idx) = mean(tmpXCovmap(idxF));
        XCovprec(idx) = sum(sum(tmpXCovmap(idxT)))/sum(sum(tmpXCovmap)); % precision = TP/(TP+FP) = TP/(inferred P)
        XCovrec(idx) = sum(sum(tmpXCovmap(idxT)))/sum(sum(gtConn)); % recall = TP/(TP+FN) = TP/(Actual P)
    end
    if params.doXCorr
        XCorrtp(idx) = mean(tmpXCorrmap(idxT));
        XCorrfp(idx) = mean(tmpXCorrmap(idxF));
        XCorrprec(idx) = sum(sum(tmpXCorrmap(idxT)))/sum(sum(tmpXCorrmap)); % precision = TP/(TP+FP) = TP/(inferred P)
        XCorrrec(idx) = sum(sum(tmpXCorrmap(idxT)))/sum(sum(gtConn)); % recall = TP/(TP+FN) = TP/(Actual P)
    end
    if params.doHO
        HOTEtp(idx) = mean(tmpHOTEmap(idxT));
        HOTEfp(idx) = mean(tmpHOTEmap(idxF));
        HOTEprec(idx) = sum(sum(tmpHOTEmap(idxT)))/sum(sum(tmpHOTEmap));
        HOTErec(idx) = sum(sum(tmpHOTEmap(idxT)))/sum(sum(gtConn));
    end
end

figure()
if params.doTE
    scatter(TEfp, TEtp)
end
hold on
if params.doXCov
    scatter(XCovfp, XCovtp)
end
if params.doXCorr
    scatter(XCorrfp, XCorrtp)
end
if params.doHO
    scatter(HOTEfp, HOTEtp)
end
plot(0:0.1:1,0:0.1:1)
legend('TE','XCov','XCorr','HOTE')
xlabel('FPR')
ylabel('TPR')
title(['TE ROC, ',num2str(minutes),' minutes'])

figure()
scatter(TErec, TEprec)
hold on
if params.doXCov
    scatter(XCovrec, XCovprec)
end
if params.doXCorr
    scatter(XCorrrec, XCorrprec)
end
if params.doHO
    scatter(HOTErec, HOTEprec)
end
yline(0.1)
ylim([0,1])
legend('TE','XCov','XCorr','HOTE')
xlabel('recall')
ylabel('precision')
title(['TE PR, ',num2str(minutes),' minutes'])


%% Correlation between real and inferred delays on pairs having a GT link
% TE delays
tmpCorr = corrcoef(staticDelaysMeasures.TEdelays(~isnan(gtDelay(:))),gtDelay(~isnan(gtDelay(:))));
delayCorr = tmpCorr(1,2);
delayErr = mean(staticDelaysMeasures.TEdelays(~isnan(gtDelay(:))) - gtDelay(~isnan(gtDelay(:))));
delayMatch = mean(staticDelaysMeasures.TEdelays(~isnan(gtDelay(:))) == gtDelay(~isnan(gtDelay(:))));

figure()
histogram(gtDelay(~isnan(gtDelay(:))),'normalization','probability')
hold on
histogram(staticDelaysMeasures.TEdelays(~isnan(gtDelay(:))),'normalization','probability')
legend('Ground truth', 'Inferred')
text(30, 0.25, ['delays correlation = ',num2str(delayCorr,2)],'FontSize',14);
text(30, 0.23, ['avg delay err = ',num2str(delayErr,2),'ms'],'FontSize',14);
text(30, 0.21, ['delays match = ',num2str(delayMatch*100),'%'],'FontSize',14);
xlim([0,50])
xlabel('Delay [ms]')
ylabel('Prob')

% HOTE delays
if params.doHO
    tmpCorr = corrcoef(staticDelaysMeasures.HOTEdelays(~isnan(gtDelay(:))),gtDelay(~isnan(gtDelay(:))));
    delayCorr = tmpCorr(1,2);
    delayErr = mean(staticDelaysMeasures.HOTEdelays(~isnan(gtDelay(:))) - gtDelay(~isnan(gtDelay(:))));
    delayMatch = mean(staticDelaysMeasures.HOTEdelays(~isnan(gtDelay(:))) == gtDelay(~isnan(gtDelay(:))));

    figure()
    histogram(gtDelay(~isnan(gtDelay(:))),'normalization','probability')
    hold on
    histogram(staticDelaysMeasures.HOTEdelays(~isnan(gtDelay(:))),'normalization','probability')
    legend('Ground truth', 'Inferred')
    text(30, 0.25, ['delays correlation = ',num2str(delayCorr,2)],'FontSize',14);
    text(30, 0.23, ['avg delay err = ',num2str(delayErr,2),'ms'],'FontSize',14);
    text(30, 0.21, ['delays match = ',num2str(delayMatch*100),'%'],'FontSize',14);
    xlim([0,50])
    xlabel('Delay [ms]')
    ylabel('Prob')
    title('Corr GT and HOTE inferred delays')
end

% XCov delays
if params.doXCov
    tmpCorr = corrcoef(staticDelaysMeasures.XCovDelays(~isnan(gtDelay(:))),gtDelay(~isnan(gtDelay(:))));
    delayCorr = tmpCorr(1,2);
    delayErr = mean(staticDelaysMeasures.XCovDelays(~isnan(gtDelay(:))) - gtDelay(~isnan(gtDelay(:))));
    delayMatch = mean(staticDelaysMeasures.XCovDelays(~isnan(gtDelay(:))) == gtDelay(~isnan(gtDelay(:))));

    figure()
    histogram(gtDelay(~isnan(gtDelay(:))),'normalization','probability')
    hold on
    histogram(staticDelaysMeasures.XCovDelays(~isnan(gtDelay(:))),'normalization','probability')
    legend('Ground truth', 'Inferred')
    text(30, 0.25, ['delays correlation = ',num2str(delayCorr,2)],'FontSize',14);
    text(30, 0.23, ['avg delay err = ',num2str(delayErr,2),'ms'],'FontSize',14);
    text(30, 0.21, ['delays match = ',num2str(delayMatch*100),'%'],'FontSize',14);
    xlim([0,50])
    xlabel('Delay [ms]')
    ylabel('Prob')
    title('Corr GT and XCov inferred delays')
end

% XCorr delays
if params.doXCorr
    tmpCorr = corrcoef(staticDelaysMeasures.XCorrDelays(~isnan(gtDelay(:))),gtDelay(~isnan(gtDelay(:))));
    delayCorr = tmpCorr(1,2);
    delayErr = mean(staticDelaysMeasures.XCorrDelays(~isnan(gtDelay(:))) - gtDelay(~isnan(gtDelay(:))));
    delayMatch = mean(staticDelaysMeasures.XCorrDelays(~isnan(gtDelay(:))) == gtDelay(~isnan(gtDelay(:))));

    figure()
    histogram(gtDelay(~isnan(gtDelay(:))),'normalization','probability')
    hold on
    histogram(staticDelaysMeasures.XCorrDelays(~isnan(gtDelay(:))),'normalization','probability')
    legend('Ground truth', 'Inferred')
    text(30, 0.25, ['delays correlation = ',num2str(delayCorr,2)],'FontSize',14);
    text(30, 0.23, ['avg delay err = ',num2str(delayErr,2),'ms'],'FontSize',14);
    text(30, 0.21, ['delays match = ',num2str(delayMatch*100),'%'],'FontSize',14);
    xlim([0,50])
    xlabel('Delay [ms]')
    ylabel('Prob')
    title('Corr GT and XCorr inferred delays')
end
save(['full_seed1_IzhiModel_v15_',num2str(minutes),'min_static_',date],'-v7.3')

% Compute dynamical functional connectivity
meanCorrTE = zeros(1,numel(new_bin_minutes_range));
meanCorrHOTE = zeros(1,numel(new_bin_minutes_range));
meanCorrXCov = zeros(1,numel(new_bin_minutes_range));
meanCorrXCorr = zeros(1,numel(new_bin_minutes_range));

% Inferred and GT synaptic weights for each measure
match_TE_time = cell(1,numel(new_bin_minutes_range));
match_HOTE_time = cell(1,numel(new_bin_minutes_range));
match_XCov_time = cell(1,numel(new_bin_minutes_range));
match_XCorr_time = cell(1,numel(new_bin_minutes_range));
match_TE_gtS_time = cell(1,numel(new_bin_minutes_range)); % Groud truth synaptic weights
match_HOTE_gtS_time = cell(1,numel(new_bin_minutes_range));
match_XCov_gtS_time = cell(1,numel(new_bin_minutes_range));
match_XCorr_gtS_time = cell(1,numel(new_bin_minutes_range));

% Correlations between gt and inferred pairs for each measure
pairCorrTE = cell(1,numel(new_bin_minutes_range));
pairCorrHOTE = cell(1,numel(new_bin_minutes_range));
pairCorrXCov = cell(1,numel(new_bin_minutes_range));
pairCorrXCorr = cell(1,numel(new_bin_minutes_range));
    
for bIdx = 1:numel(new_bin_minutes_range)
    new_bin_minutes = new_bin_minutes_range(bIdx);
    disp(['Time window ',num2str(new_bin_minutes),' minutes'])
    
    new_bin = new_bin_minutes*60*1000;
    time_jump = time_jump_minutes*60*1000;

    % TE selected network
    TEthr = prctile(staticConnMeasures.peakTE(:),sel_thresh);
    TEmap = (staticConnMeasures.peakTE > TEthr);
    G = digraph(TEmap);
    TE_connPairs = G.Edges.EndNodes; % Pairs that we inferred as connected and for which we will compute TE
    % Idx 1 is row, idx 2 is column of TEmap
    TE_connDelay = zeros(1,size(TE_connPairs,1)); % Inferred delays for the putatively connected pairs
    for pIdx = 1:size(TE_connPairs,1)
        TE_connDelay(pIdx) = staticDelaysMeasures.TEdelays(TE_connPairs(pIdx,1),TE_connPairs(pIdx,2));
    end
    
    % HOTE selected network
    HOTE_connPairs = zeros(size(TE_connPairs));
    HOTE_connDelay = zeros(size(TE_connDelay));
    if params.doHO
        HOthr = prctile(staticConnMeasures.peakHOTE(:),sel_thresh);
        HOTEmap = (staticConnMeasures.peakHOTE > HOthr);
        G = digraph(HOTEmap);
        HOTE_connPairs = G.Edges.EndNodes; % Pairs that we inferred as connected and for which we will compute TE
        % Idx 1 is row, idx 2 is column of TEmap
        HOTE_connDelay = zeros(1,size(HOTE_connPairs,1)); % Inferred delays for the putatively connected pairs
        for pIdx = 1:size(HOTE_connPairs,1)
            HOTE_connDelay(pIdx) = staticDelaysMeasures.HOTEdelays(HOTE_connPairs(pIdx,1),HOTE_connPairs(pIdx,2));
        end
    end
    
    % XCov selected network
    XCov_connPairs = zeros(size(TE_connPairs));
    XCov_connDelay = zeros(size(TE_connDelay));
    if params.doXCov
        XCovthr = prctile(staticConnMeasures.peakXCov(:),sel_thresh);
        XCmap = (staticConnMeasures.peakXCov > XCovthr);
        G = digraph(XCmap);
        XCov_connPairs = G.Edges.EndNodes; % Pairs that we inferred as connected and for which we will compute TE
        % Idx 1 is row, idx 2 is column of TEmap
        XCov_connDelay = zeros(1,size(XCov_connPairs,1)); % Inferred delays for the putatively connected pairs
        for pIdx = 1:size(XCov_connPairs,1)
            XCov_connDelay(pIdx) = staticDelaysMeasures.XCovDelays(XCov_connPairs(pIdx,1),XCov_connPairs(pIdx,2));
        end
    end
    % XCorr selected network
    XCorr_connPairs = zeros(size(TE_connPairs));
    XCorr_connDelay = zeros(size(TE_connDelay));
    if params.doXCorr
        XCorrthr = prctile(staticConnMeasures.peakXCorr(:),sel_thresh);
        XCmap = (staticConnMeasures.peakXCorr > XCorrthr);
        G = digraph(XCmap);
        XCorr_connPairs = G.Edges.EndNodes; % Pairs that we inferred as connected and for which we will compute TE
        % Idx 1 is row, idx 2 is column of TEmap
        XCorr_connDelay = zeros(1,size(XCorr_connPairs,1)); % Inferred delays for the putatively connected pairs
        for pIdx = 1:size(XCorr_connPairs,1)
            XCorr_connDelay(pIdx) = staticDelaysMeasures.XCorrDelays(XCorr_connPairs(pIdx,1),XCorr_connPairs(pIdx,2));
        end
    end

    %%Compute dynamical functional connectivity
    % Consider eliminating from firings / TE_time all neurons that aren't in
    % any significant pair
    [TE_time,HOTE_time,XCov_time,XCorr_time] = movWin_connectivity_v2(spikeTrains,TE_connPairs,TE_connDelay,...
                                        params.doHO,HOTE_connPairs,HOTE_connDelay,...
                                        params.doXCov,XCov_connPairs,XCov_connDelay,params.max_delay,...
                                        params.doXCorr,XCorr_connPairs,XCorr_connDelay,...
                                        delay_err,new_bin,time_jump,delay_consistency,params.norm_method);

    %
    % Average the synaptic weight over windows of new_bin length with time_jump
    s_time_perm = permute(s_time,[1,3,2]);
    s_time_downsampled = temporal_rebinning(s_time_perm,new_bin/1000,'movmean',time_jump/1000);
    s_time_downsampled = permute(s_time_downsampled,[1,3,2]);

    % Find matching pairs 
    match_pairs.TE = [];
    match_pairs.HOTE = [];
    match_pairs.XCov = [];
    match_pairs.XCorr = [];
    
    for pairIdx = 1:size(TE_connPairs) % diff measures have same size
        idxs = [TE_connPairs(pairIdx,1),TE_connPairs(pairIdx,2)];
        if gtConn(idxs(1),idxs(2)) == 1
            match_pairs.TE = [match_pairs.TE, pairIdx];
        end
        if params.doHO
            idxs = [HOTE_connPairs(pairIdx,1),HOTE_connPairs(pairIdx,2)];
            if gtConn(idxs(1),idxs(2)) == 1
                match_pairs.HOTE = [match_pairs.HOTE, pairIdx];
            end
        end
        if params.doXCov
            idxs = [XCov_connPairs(pairIdx,1),XCov_connPairs(pairIdx,2)];
            if gtConn(idxs(1),idxs(2)) == 1
                match_pairs.XCov = [match_pairs.XCov, pairIdx];
            end
        end
        if params.doXCorr
            idxs = [XCorr_connPairs(pairIdx,1),XCorr_connPairs(pairIdx,2)];
            if gtConn(idxs(1),idxs(2)) == 1
                match_pairs.XCorr = [match_pairs.XCorr, pairIdx];
            end
        end
    end
    fraction_match.TE = numel(match_pairs.TE)/size(TE_time,1); % this should match the TPR at the chosen threshold
    fraction_match.HOTE = numel(match_pairs.HOTE)/size(HOTE_time,1);
    fraction_match.XCov = numel(match_pairs.XCov)/size(XCov_time,1);
    fraction_match.XCorr = numel(match_pairs.XCorr)/size(XCov_time,1);

    % Rearrange ground truth changes in connectivity for the matched (TP) pairs
    match_TE_time{bIdx} = TE_time(match_pairs.TE,:);
    match_HOTE_time{bIdx} = HOTE_time(match_pairs.HOTE,:);
    match_XCov_time{bIdx} = XCov_time(match_pairs.XCov,:);
    match_XCorr_time{bIdx} = XCorr_time(match_pairs.XCorr,:);

    % GT synaptic weights matching the ones found in different measures
    match_TE_gtS_time{bIdx} = zeros(numel(match_pairs.TE),size(TE_time,2)); % Groud truth synaptic weights
    match_HOTE_gtS_time{bIdx} = zeros(numel(match_pairs.HOTE),size(HOTE_time,2));
    match_XCov_gtS_time{bIdx} = zeros(numel(match_pairs.XCov),size(XCov_time,2));
    match_XCorr_gtS_time{bIdx} = zeros(numel(match_pairs.XCorr),size(XCov_time,2));
    
    % Correlations between FC and SC for each pair
    pairCorrTE{bIdx} = zeros(1,numel(match_pairs.TE));
    pairCorrHOTE{bIdx} = zeros(1,numel(match_pairs.HOTE));
    pairCorrXCov{bIdx} = zeros(1,numel(match_pairs.XCov));
    pairCorrXCorr{bIdx} = zeros(1,numel(match_pairs.XCorr));

    for pairIdx = 1:numel(match_pairs.TE)
        % Selecting idxs for the matched pair 
        idxs = [TE_connPairs(match_pairs.TE(pairIdx),1),TE_connPairs(match_pairs.TE(pairIdx),2)];
        idx2 = find(post(idxs(1),:) == idxs(2));
        tmp = s_time_downsampled(idxs(1),idx2,:);
        match_TE_gtS_time{bIdx}(pairIdx,:) = tmp;
    end
    for pairIdx = 1:numel(match_pairs.HOTE)
        % Selecting idxs for the matched pair 
        idxs = [HOTE_connPairs(match_pairs.HOTE(pairIdx),1),HOTE_connPairs(match_pairs.HOTE(pairIdx),2)];
        idx2 = find(post(idxs(1),:) == idxs(2));
        tmp = s_time_downsampled(idxs(1),idx2,:);
        match_HOTE_gtS_time{bIdx}(pairIdx,:) = tmp;
    end
    for pairIdx = 1:numel(match_pairs.XCov)
        % Selecting idxs for the matched pair 
        idxs = [XCov_connPairs(match_pairs.XCov(pairIdx),1),XCov_connPairs(match_pairs.XCov(pairIdx),2)];
        idx2 = find(post(idxs(1),:) == idxs(2));
        tmp = s_time_downsampled(idxs(1),idx2,:);
        match_XCov_gtS_time{bIdx}(pairIdx,:) = tmp;
    end
    for pairIdx = 1:numel(match_pairs.XCorr)
        % Selecting idxs for the matched pair 
        idxs = [XCorr_connPairs(match_pairs.XCorr(pairIdx),1),XCorr_connPairs(match_pairs.XCorr(pairIdx),2)];
        idx2 = find(post(idxs(1),:) == idxs(2));
        tmp = s_time_downsampled(idxs(1),idx2,:);
        match_XCorr_gtS_time{bIdx}(pairIdx,:) = tmp;
    end

    if new_bin_minutes == new_bin_minutes_range(end)
        % Plot and compute accuracy metrics
        x_axis=1:(new_bin_minutes/time_jump_minutes):size(match_TE_time{bIdx},2);
        x_labels=time_jump_minutes:new_bin_minutes:minutes-time_jump_minutes;
        % Plot TE DFC
        figure()
        subplot(1,2,1)
        imagesc(match_TE_gtS_time{bIdx}(1:20,:))
        colorbar()
        set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
        xlabel('minutes')
        ylabel('pairs')
        title('Ground truth')

        subplot(1,2,2)
        imagesc(match_TE_time{bIdx}(1:20,:))
        colorbar()
        set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
        xlabel('minutes')
        ylabel('pairs')
        title('Inferred (TE)')
        sgtitle(['TE synaptic weights over time, ',num2str(new_bin_minutes),'min window, ',num2str(time_jump_minutes),'min time jump']);

        %Plot HOTE DFC
        if params.doHO
            subplot(1,2,1)
            imagesc(match_HOTE_gtS_time{bIdx}(1:20,:))
            colorbar()
            set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
            xlabel('minutes')
            ylabel('pairs')
            title('Ground truth')

            subplot(1,2,2)
            imagesc(match_HOTE_time{bIdx}(1:20,:))
            colorbar()
            set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
            xlabel('minutes')
            ylabel('pairs')
            title('Inferred (HOTE)')
            sgtitle(['HOTE synaptic weights over time, ',num2str(new_bin_minutes),'min window, ',num2str(time_jump_minutes),'min time jump']);
        end
        
        % Plot XCov GT
        if params.doXCov
            subplot(1,2,1)
            imagesc(match_XCov_gtS_time{bIdx}(1:20,:))
            colorbar()
            set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
            xlabel('minutes')
            ylabel('pairs')
            title('Ground truth')

            subplot(1,2,2)
            imagesc(match_XCov_time{bIdx}(1:20,:))
            colorbar()
            set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
            xlabel('minutes')
            ylabel('pairs')
            title('Inferred (HOTE)')
            sgtitle(['XC synaptic weights over time, ',num2str(new_bin_minutes),'min window, ',num2str(time_jump_minutes),'min time jump']);
        end
        
        % Plot XCorr GT
        if params.doXCorr
            subplot(1,2,1)
            imagesc(match_XCorr_gtS_time{bIdx}(1:20,:))
            colorbar()
            set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
            xlabel('minutes')
            ylabel('pairs')
            title('Ground truth')

            subplot(1,2,2)
            imagesc(match_XCorr_time{bIdx}(1:20,:))
            colorbar()
            set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
            xlabel('minutes')
            ylabel('pairs')
            title('Inferred (HOTE)')
            sgtitle(['XC synaptic weights over time, ',num2str(new_bin_minutes),'min window, ',num2str(time_jump_minutes),'min time jump']);
        end
    end

    % Compute correlation between DFC and DSC
    for pairIdx = 1:numel(match_pairs.TE)  % average correlation on true positives
        tmp = corrcoef(match_TE_gtS_time{bIdx}(pairIdx,:),match_TE_time{bIdx}(pairIdx,:));
        pairCorrTE{bIdx}(pairIdx) = tmp(1,2);
        meanCorrTE(bIdx) = nanmean(pairCorrTE{bIdx});
    end
    if params.doHO
        for pairIdx = 1:numel(match_pairs.HOTE)
            tmp = corrcoef(match_HOTE_gtS_time{bIdx}(pairIdx,:),match_HOTE_time{bIdx}(pairIdx,:));
            pairCorrHOTE{bIdx}(pairIdx) = tmp(1,2);
            meanCorrHOTE(bIdx) = nanmean(pairCorrHOTE{bIdx});
        end
    end
    if params.doXCov
        for pairIdx = 1:numel(match_pairs.XCov)
            tmp = corrcoef(match_XCov_gtS_time{bIdx}(pairIdx,:),match_XCov_time{bIdx}(pairIdx,:));
            pairCorrXCov{bIdx}(pairIdx) = tmp(1,2);
            meanCorrXCov(bIdx) = nanmean(pairCorrXCov{bIdx});
        end
    end
    if params.doXCorr
        for pairIdx = 1:numel(match_pairs.XCorr)
            tmp = corrcoef(match_XCorr_gtS_time{bIdx}(pairIdx,:),match_XCorr_time{bIdx}(pairIdx,:));
            pairCorrXCorr{bIdx}(pairIdx) = tmp(1,2);
            meanCorrXCorr(bIdx) = nanmean(pairCorrXCorr{bIdx});
        end
    end

    if new_bin_minutes == new_bin_minutes_range(end)
        figure();
        histogram(pairCorrTE{bIdx},20,'normalization','probability')
        xlabel('correlation')
        title(['Histogram of corr(real, TE) synaptic weight over time, time window ',num2str(new_bin_minutes),'min'])
        text(meanCorrTE(end),0.1,['mean ',num2str(meanCorrTE(end))]);

        if params.doHO
            figure();
            histogram(pairCorrHOTE{bIdx},20,'normalization','probability')
            xlabel('correlation')
            title(['Histogram of corr(real, HOTE) synaptic weight over time, time window ',num2str(new_bin_minutes),'min'])
            text(meanCorrHOTE(end),0.1,['mean ',num2str(meanCorrHOTE(end))]);
        end
        if params.doXCov
            figure();
            histogram(pairCorrXCov{bIdx},20,'normalization','probability')
            xlabel('correlation')
            title(['Histogram of corr(real, XCov) synaptic weight over time, time window ',num2str(new_bin_minutes),'min'])
            text(meanCorrXCov(end),0.1,['mean ',num2str(meanCorrXCov(end))]);
        end
        if params.doXCorr
            figure();
            histogram(pairCorrXCorr{bIdx},20,'normalization','probability')
            xlabel('correlation')
            title(['Histogram of corr(real, XC) synaptic weight over time, time window ',num2str(new_bin_minutes),'min'])
            text(meanCorrXCorr(end),0.1,['mean ',num2str(meanCorrXCorr(end))]);
        end
    end
end

fname = ['dynamicConnInference_rep',num2str(randSeed),'_',num2str(minutes),'min_varyWindow_',num2str(delay_consistency),'_DelConsist_',date];
save(fname,'-v7.3')
