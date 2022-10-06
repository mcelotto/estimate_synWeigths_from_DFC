%% Script to reproduce plots in Figure 5

clear all; clc;
PATH = ['path_to_main_directory\results\dynamic_connectivity\']; % path to directory where the dynamic connectivity results have been saved

norm_method='normalized'; % 'none' or 'normalized'
nReps = 1;
new_bin_minutes_range = [1 2 5 10 20 30]; % sizes of sliding windows
time_jumps = new_bin_minutes_range;
doHO = 1; doXCov = 1; doXCorr = 1;
% Load simulation results
for repIdx = 1:nReps
    % Load delay-consistent results
    files = dir([PATH,'dynamicConnInference_rep',num2str(repIdx),'_180min_varyWindow_1_DelConsist*']);
    fname = files.name;
    dc{1} = load([PATH,fname],'pairCorrTE','match_TE_gtS_time','match_TE_time','pairCorrHOTE','match_HOTE_gtS_time','match_HOTE_time','pairCorrXCov','match_XCov_gtS_time','match_XCov_time','doXCov','doXCorr','meanCorrXCov','XCov_connPairs','x_axis','x_labels','time_jump_minutes','meanCorrTE','match_pairs','TE_connPairs','meanCorrHOTE','HOTE_connPairs','post','s_time','gtConn');
    % Load non delay-consistent results
    files = dir([PATH,'dynamicConnInference_rep',num2str(repIdx),'_180min_varyWindow_0_DelConsist*']);
    fname = files.name;
    dc{2} = load([PATH,fname],'pairCorrTE','match_TE_gtS_time','match_TE_time','pairCorrHOTE','match_HOTE_gtS_time','match_HOTE_time','pairCorrXCov','match_XCov_gtS_time','match_XCov_time','doXCov','doXCorr','meanCorrXCov','XCov_connPairs','x_axis','x_labels','time_jump_minutes','meanCorrTE','match_pairs','TE_connPairs','meanCorrHOTE','HOTE_connPairs');

    % store GT connectivity, post=synaptic neurons identity and GT synaptic
    % weight evolution for this repetition of the simulation
    gtConn{repIdx} = dc{1}.gtConn;
    post{repIdx} = dc{1}.post;
    s_time{repIdx} = dc{1}.s_time;
    connPairs{repIdx}.TE = dc{1}.TE_connPairs;
    connPairs{repIdx}.HOTE = dc{1}.HOTE_connPairs;
    connPairs{repIdx}.XCov = dc{1}.XCov_connPairs;
    
    % XCov selected network:
    % In the main script I accidentally forgot to take the absolute value
    % of XCov before selecting the synapses on which I computed DFC. Therefore
    % I selected only excitatory synapses for XCov, while I neglected inhibitory
    % synapses (N=200). Here, I matched the number of XCov synapses to the excitatory 
    % ones computed using TE and HOTE
  
    rng('default')
    selected_pairs = min([sum(~isnan(dc{1}.pairCorrHOTE{1})),sum(~isnan(dc{1}.pairCorrTE{1}))]);
    rIdx = randperm(selected_pairs);
    rIdx = rIdx(1:selected_pairs); 
    dc{1}.match_pairs.XCov = dc{1}.match_pairs.XCov(rIdx);
    dc{2}.match_pairs.XCov = dc{2}.match_pairs.XCov(rIdx);
    connPairs{repIdx}.XCov = connPairs{repIdx}.XCov(dc{1}.match_pairs.XCov,:);

    % Store delay-consistent and non delay-consistent results in common
    % data structures
    dc_label = {'dc1','dc0'}; % delay-consistency structure fieldnames
    for dcIdx=1:2
        dcLab=dc_label{dcIdx};
        
        tmppairCorrTE.(dcLab) = dc{dcIdx}.pairCorrTE;
        tmppairCorrHOTE.(dcLab) = dc{dcIdx}.pairCorrHOTE;
        tmppairCorrXCov.(dcLab) = dc{dcIdx}.pairCorrXCov;
        tmpMeanCorrTE.(dcLab) = dc{dcIdx}.meanCorrTE;
        tmpMeanCorrHOTE.(dcLab) = dc{dcIdx}.meanCorrHOTE;
        tmpMeanCorrXCov.(dcLab) = dc{dcIdx}.meanCorrXCov;

        for bIdx = 1:numel(new_bin_minutes_range)
            pairCorrTE{repIdx}.(dcLab).minTJ{bIdx}=tmppairCorrTE.(dcLab){bIdx};
            pairCorrHOTE{repIdx}.(dcLab).minTJ{bIdx}=tmppairCorrHOTE.(dcLab){bIdx};
            pairCorrXCov{repIdx}.(dcLab).minTJ{bIdx}=tmppairCorrXCov.(dcLab){bIdx}(rIdx);

        end
        meanCorrTE{repIdx}.(dcLab).minTJ=tmpMeanCorrTE.(dcLab);
        meanCorrHOTE{repIdx}.(dcLab).minTJ=tmpMeanCorrHOTE.(dcLab);
        meanCorrXCov{repIdx}.(dcLab).minTJ=tmpMeanCorrXCov.(dcLab);

        % Initialize matched_connections
        for bIdx = 1:numel(new_bin_minutes_range)
            match_XCov{repIdx}.(dcLab){bIdx} = dc{dcIdx}.match_XCov_time{bIdx}(rIdx,:);
            match_XCov_gtS.(dcLab){bIdx} = dc{dcIdx}.match_XCov_gtS_time{bIdx}(rIdx,:);
        end
        match_TE{repIdx}.(dcLab) = dc{dcIdx}.match_TE_time;
        match_TE_gtS{repIdx}.(dcLab) = dc{dcIdx}.match_TE_gtS_time;
        match_HOTE{repIdx}.(dcLab) = dc{dcIdx}.match_HOTE_time;
        match_HOTE_gtS{repIdx}.(dcLab) = dc{dcIdx}.match_HOTE_gtS_time;
        
        matched_pairs{repIdx}.(dcLab) = dc{dcIdx}.match_pairs;
    end

    % Common pairs to all measures
    tmpComm_pairs = intersect(connPairs{repIdx}.TE,connPairs{repIdx}.HOTE,'rows');
    Comm_pairs{repIdx} = intersect(tmpComm_pairs,connPairs{repIdx}.XCov,'rows');
    % discard non GT common pairs
    for idx = 1:size(Comm_pairs{repIdx},1)
        if ~gtConn{repIdx}(Comm_pairs{repIdx}(idx,1),Comm_pairs{repIdx}(idx,2))
            Comm_pairs{repIdx}(idx,:)=[];
        end
    end

% Comm_pairs{repIdx} = intersect(TE_connPairs,XCov_connPairs,'rows');
    %% Compute non-autocorr weights
    doHO = 1; doXCov = 1; doXCorr = 0;
    for dcIdx=1:2
        dcLab=dc_label{dcIdx};
        for bIdx = 1:numel(new_bin_minutes_range)
            new_bin_minutes = new_bin_minutes_range(bIdx);

            for TJ_idx = 1:numel(time_jumps)
                TJ = time_jumps(TJ_idx);
                t_idxs = 1:(TJ)/(time_jumps(1)):size(match_TE_gtS{repIdx}.(dcLab){bIdx},2);

                for pairIdx = 1:numel(pairCorrTE{repIdx}.(dcLab).minTJ{bIdx})  % average correlation on true positives
                    tmp = corrcoef(match_TE_gtS{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs),match_TE{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs));
                    pairCorrTE{repIdx}.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    tmp = corrcoef(diff(match_TE_gtS{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs)),diff(match_TE{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs)));
                    pairCorrTEdiff{repIdx}.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                end
                meanCorrTE{repIdx}.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrTE{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                meanCorrTEdiff{repIdx}.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrTEdiff{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                semCorrTE{repIdx}.(dcLab).varyTJ(bIdx,TJ_idx) = std(pairCorrTE{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrTE{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx))));
                %semCorrTE{repIdx}.(dcLab).minTJ(bIdx) = std(pairCorrTE{repIdx}.(dcLab).minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrTE{repIdx}.(dcLab).minTJ{bIdx})));

                if doHO
                    for pairIdx = 1:numel(pairCorrHOTE{repIdx}.(dcLab).minTJ{bIdx})
                        tmp = corrcoef(match_HOTE_gtS{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs),match_HOTE{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs));
                        pairCorrHOTE{repIdx}.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                        tmp = corrcoef(diff(match_HOTE_gtS{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs)),diff(match_HOTE{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs)));
                        pairCorrHOTEdiff{repIdx}.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    end
                    meanCorrHOTE{repIdx}.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrHOTE{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                    meanCorrHOTEdiff{repIdx}.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrHOTEdiff{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                    semCorrHOTE{repIdx}.(dcLab).varyTJ(bIdx,TJ_idx) = std(pairCorrHOTE{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrHOTE{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx))));
                    %semCorrHOTE{repIdx}.(dcLab).minTJ(bIdx) = std(pairCorrHOTE{repIdx}.(dcLab).minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrHOTE{repIdx}.(dcLab).minTJ{bIdx})));
                end
                if doXCov
                    for pairIdx = 1:numel(pairCorrXCov{repIdx}.(dcLab).minTJ{bIdx})
                        tmp = corrcoef(match_XCov_gtS.(dcLab){bIdx}(pairIdx,t_idxs),match_XCov{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs));
                        pairCorrXCov{repIdx}.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                        tmp = corrcoef(diff(match_XCov_gtS.(dcLab){bIdx}(pairIdx,t_idxs)),diff(match_XCov{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs)));
                        pairCorrXCovdiff{repIdx}.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    end
                    meanCorrXCov{repIdx}.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrXCov{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                    meanCorrXCovdiff{repIdx}.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrXCovdiff{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                    semCorrXCov{repIdx}.(dcLab).varyTJ(bIdx,TJ_idx) = std(pairCorrXCov{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrXCov{repIdx}.(dcLab).varyTJ{bIdx}(:,TJ_idx))));
                    %semCorrXCov{repIdx}.(dcLab).minTJ(bIdx) = std(pairCorrXCov{repIdx}.(dcLab).minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrXCov{repIdx}.(dcLab).minTJ{bIdx})));
                end
            end
        end
    end

    %% Compute sample-matched DFC and DSC correlations

    max_minutes = max(new_bin_minutes_range);
    sIdx = find(new_bin_minutes_range==max_minutes);
    new_bin_minutes = new_bin_minutes_range(sIdx);
    maxTJ = time_jumps(end);
    n_samples = numel(1:maxTJ/(time_jumps(1)):size(match_TE_gtS{repIdx}.(dcLab){sIdx},2)); % set size of 10 no autocorrelated

    for dcIdx=1:2
        dcLab=dc_label{dcIdx};
        for bIdx = 1:sIdx
            new_bin_minutes = new_bin_minutes_range(bIdx);
            % Take time points so that they're uncorrelated
            for TJ_idx = 1:numel(time_jumps)
                TJ = time_jumps(TJ_idx);
                t_idxs = 1:(TJ)/(time_jumps(1)):size(match_TE_gtS{repIdx}.(dcLab){bIdx},2);

                nDiff = numel(t_idxs)-n_samples;
                t_idxs(1:floor(nDiff/2))=[];
                t_idxs(end-ceil(nDiff/2)+1:end)=[];

                %for idx = 1:size(TE_connPairs,1)
                for idx = matched_pairs{repIdx}.(dcLab).TE
                    ridx = find(dc{1}.post(connPairs{repIdx}.TE(idx,1),:)==connPairs{repIdx}.TE(idx,2));
                    s_time_TE{repIdx}(idx,:) = dc{1}.s_time(connPairs{repIdx}.TE(idx,1),ridx,:);
                end
                s_time_dawnsamp = temporal_rebinning(s_time_TE{repIdx},TJ*60,'movmean',TJ*60);
                nDiff = size(s_time_dawnsamp,2)-n_samples;
                t_idxs_gt = 1:size(s_time_dawnsamp,2);
                t_idxs_gt(1:floor(nDiff/2))=[];
                t_idxs_gt(end-ceil(nDiff/2)+1:end)=[];

                for pairIdx = 1:numel(pairCorrTE{repIdx}.(dcLab).minTJ{bIdx})  % average correlation on true positives
                    %tmp = corrcoef(match_TE_gtS{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs),match_TE{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs));
                    tmp = corrcoef(s_time_dawnsamp(pairIdx,t_idxs_gt),match_TE{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs));
                    pairCorrTE{repIdx}.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    meanCorrTE{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrTE{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                    semCorrTE{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrTE{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrTE{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                    tmp = corrcoef(diff(s_time_dawnsamp(pairIdx,t_idxs_gt)),diff(match_TE{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs)));
                    pairCorrTEdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    meanCorrTEdiff{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrTEdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                    semCorrTEdiff{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrTEdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrTEdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                end
                if doHO
                    %for idx = 1:size(HOTE_connPairs,1)
                    for idx = matched_pairs{repIdx}.(dcLab).HOTE
                        ridx = find(dc{1}.post(connPairs{repIdx}.HOTE(idx,1),:)==connPairs{repIdx}.HOTE(idx,2));
                        s_time_HOTE{repIdx}(idx,:) = dc{1}.s_time(connPairs{repIdx}.HOTE(idx,1),ridx,:);
                    end
                    s_time_dawnsamp = temporal_rebinning(s_time_HOTE{repIdx},TJ*60,'movmean',TJ*60);
                    nDiff = size(s_time_dawnsamp,2)-n_samples;
                    t_idxs_gt = 1:size(s_time_dawnsamp,2);
                    t_idxs_gt(1:floor(nDiff/2))=[];
                    t_idxs_gt(end-ceil(nDiff/2)+1:end)=[];

                    for pairIdx = 1:numel(pairCorrHOTE{repIdx}.(dcLab).minTJ{bIdx})
                        tmp = corrcoef(s_time_dawnsamp(pairIdx,t_idxs_gt),match_HOTE{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs));
                        pairCorrHOTE{repIdx}.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                        meanCorrHOTE{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrHOTE{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                        semCorrHOTE{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrHOTE{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrHOTE{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                        tmp = corrcoef(diff(s_time_dawnsamp(pairIdx,t_idxs_gt)),diff(match_HOTE{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs)));
                        pairCorrHOTEdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                        meanCorrHOTEdiff{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrHOTEdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                        semCorrHOTEdiff{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrHOTEdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrHOTEdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                    end
                end
                if doXCov
                    %for idx = 1:size(XCov_connPairs,1)
                    for idx = 1:numel(matched_pairs{repIdx}.(dcLab).XCov)
                        ridx = find(dc{1}.post(connPairs{repIdx}.XCov(idx,1),:)==connPairs{repIdx}.XCov(idx,2));
                        s_time_XCov{repIdx}(idx,:) = dc{1}.s_time(connPairs{repIdx}.XCov(idx,1),ridx,:);
                    end
                    s_time_dawnsamp = temporal_rebinning(s_time_XCov{repIdx},TJ*60,'movmean',TJ*60);
                    nDiff = size(s_time_dawnsamp,2)-n_samples;
                    t_idxs_gt = 1:size(s_time_dawnsamp,2);
                    t_idxs_gt(1:floor(nDiff/2))=[];
                    t_idxs_gt(end-ceil(nDiff/2)+1:end)=[];
                    for pairIdx = 1:numel(pairCorrXCov{repIdx}.(dcLab).minTJ{bIdx})
                        tmp = corrcoef(s_time_dawnsamp(pairIdx,t_idxs_gt),match_XCov{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs));
                        pairCorrXCov{repIdx}.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                        meanCorrXCov{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrXCov{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                        semCorrXCov{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrXCov{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrXCov{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                        tmp = corrcoef(diff(s_time_dawnsamp(pairIdx,t_idxs_gt)),diff(match_XCov{repIdx}.(dcLab){bIdx}(pairIdx,t_idxs)));
                        pairCorrXCovdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                        meanCorrXCovdiff{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrXCovdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                        semCorrXCovdiff{repIdx}.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrXCovdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrXCovdiff{repIdx}.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                    end
                end
            end
        end
    end

    % Compute SEM for difference of dc1 and dc0 (the pairs are matched since
    % the delay consistency only affects the DFC and not the DSC)
    for bIdx = 1:numel(new_bin_minutes_range)
        semCorrTE{repIdx}.diff.minTJ(bIdx) = std(pairCorrTE{repIdx}.dc1.minTJ{bIdx}-pairCorrTE{repIdx}.dc0.minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrTE{repIdx}.(dcLab).minTJ{bIdx})));
        semCorrTE{repIdx}.diff.varyTJ(bIdx) = std(pairCorrTE{repIdx}.dc1.varyTJ{bIdx}-pairCorrTE{repIdx}.dc0.varyTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrTE{repIdx}.(dcLab).varyTJ{bIdx})));
        semCorrHOTE{repIdx}.diff.minTJ(bIdx) = std(pairCorrHOTE{repIdx}.dc1.minTJ{bIdx}-pairCorrHOTE{repIdx}.dc0.minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrHOTE{repIdx}.(dcLab).minTJ{bIdx})));
        semCorrHOTE{repIdx}.diff.varyTJ(bIdx) = std(pairCorrHOTE{repIdx}.dc1.varyTJ{bIdx}-pairCorrHOTE{repIdx}.dc0.varyTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrHOTE{repIdx}.(dcLab).varyTJ{bIdx})));
        semCorrXCov{repIdx}.diff.minTJ(bIdx) = std(pairCorrXCov{repIdx}.dc1.minTJ{bIdx}-pairCorrXCov{repIdx}.dc0.minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrXCov{repIdx}.(dcLab).minTJ{bIdx})));
        semCorrXCov{repIdx}.diff.varyTJ(bIdx) = std(pairCorrXCov{repIdx}.dc1.varyTJ{bIdx}-pairCorrXCov{repIdx}.dc0.varyTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrXCov{repIdx}.(dcLab).varyTJ{bIdx})));
    end

    for bIdx = 1:sIdx
        semCorrTE{repIdx}.diff.NsampleMatch(bIdx) = std(pairCorrTE{repIdx}.dc1.NsampleMatch{bIdx}-pairCorrTE{repIdx}.dc0.NsampleMatch{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrTE{repIdx}.(dcLab).NsampleMatch{bIdx})));
        semCorrHOTE{repIdx}.diff.NsampleMatch(bIdx) = std(pairCorrHOTE{repIdx}.dc1.NsampleMatch{bIdx}-pairCorrHOTE{repIdx}.dc0.NsampleMatch{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrHOTE{repIdx}.(dcLab).NsampleMatch{bIdx})));
        semCorrXCov{repIdx}.diff.NsampleMatch(bIdx) = std(pairCorrXCov{repIdx}.dc1.NsampleMatch{bIdx}-pairCorrXCov{repIdx}.dc0.NsampleMatch{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrXCov{repIdx}.(dcLab).NsampleMatch{bIdx})));
    end
end

%% pool correlation values across repetitions
DFC_types = fieldnames(pairCorrTE{repIdx}.dc1);

for DFCidx = [2 3] % minTJ is not used and has different size
    DFCLab = DFC_types{DFCidx};
    for dcIdx = 1:2
    	dcLab = dc_label{dcIdx};
        
        pooledMeanTE.(dcLab).(DFCLab) = zeros(numel(new_bin_minutes_range),numel(new_bin_minutes_range)); % pairwise correlations 
        pooledMeanHOTE.(dcLab).(DFCLab) = zeros(numel(new_bin_minutes_range),numel(new_bin_minutes_range));
        pooledMeanXCov.(dcLab).(DFCLab) = zeros(numel(new_bin_minutes_range),numel(new_bin_minutes_range));
        
        pooledSemTE.(dcLab).(DFCLab) = zeros(numel(new_bin_minutes_range),numel(new_bin_minutes_range)); % pairwise correlations 
        pooledSemHOTE.(dcLab).(DFCLab) = zeros(numel(new_bin_minutes_range),numel(new_bin_minutes_range));
        pooledSemXCov.(dcLab).(DFCLab) = zeros(numel(new_bin_minutes_range),numel(new_bin_minutes_range));
        
        for bIdx = 1:numel(new_bin_minutes_range)
            pooledCorrTE{bIdx}.(dcLab).(DFCLab) = []; % pairwise correlations
            pooledCorrHOTE{bIdx}.(dcLab).(DFCLab) = [];
            pooledCorrXCov{bIdx}.(dcLab).(DFCLab) = [];
            
            for repIdx = 1:numel(nReps) 
                pooledCorrTE{bIdx}.(dcLab).(DFCLab) = [pooledCorrTE{bIdx}.(dcLab).(DFCLab); pairCorrTE{repIdx}.(dcLab).(DFCLab){bIdx}];
                pooledCorrHOTE{bIdx}.(dcLab).(DFCLab) = [pooledCorrHOTE{bIdx}.(dcLab).(DFCLab); pairCorrHOTE{repIdx}.(dcLab).(DFCLab){bIdx}];
                pooledCorrXCov{bIdx}.(dcLab).(DFCLab) = [pooledCorrXCov{bIdx}.(dcLab).(DFCLab); pairCorrXCov{repIdx}.(dcLab).(DFCLab){bIdx}];
            end
        
            % pooled mean
            pooledMeanTE.(dcLab).(DFCLab)(bIdx,:) = nanmean(pooledCorrTE{bIdx}.(dcLab).(DFCLab),1);
            pooledMeanHOTE.(dcLab).(DFCLab)(bIdx,:) = nanmean(pooledCorrHOTE{bIdx}.(dcLab).(DFCLab),1);
            pooledMeanXCov.(dcLab).(DFCLab)(bIdx,:) = nanmean(pooledCorrXCov{bIdx}.(dcLab).(DFCLab),1);
            % pooled SEM
            pooledSemTE.(dcLab).(DFCLab)(bIdx,:) = std(pooledCorrTE{bIdx}.(dcLab).(DFCLab),[],1,'omitnan')./sqrt(sum(~isnan(pooledCorrTE{bIdx}.(dcLab).(DFCLab)),1));
            pooledSemHOTE.(dcLab).(DFCLab)(bIdx,:) = std(pooledCorrHOTE{bIdx}.(dcLab).(DFCLab),[],1,'omitnan')./sqrt(sum(~isnan(pooledCorrHOTE{bIdx}.(dcLab).(DFCLab)),1));
            pooledSemXCov.(dcLab).(DFCLab)(bIdx,:) = std(pooledCorrXCov{bIdx}.(dcLab).(DFCLab),[],1,'omitnan')./sqrt(sum(~isnan(pooledCorrXCov{bIdx}.(dcLab).(DFCLab)),1));

        end
    end
    for bIdx = 1:numel(new_bin_minutes_range)
        pooledSemTE.diff.(DFCLab)(bIdx,:) = std(pooledCorrTE{bIdx}.dc1.(DFCLab)-pooledCorrTE{bIdx}.dc0.(DFCLab),'omitnan')/sqrt(sum(~isnan(pooledCorrTE{bIdx}.(dcLab).(DFCLab))));
        pooledSemHOTE.diff.(DFCLab)(bIdx,:) = std(pooledCorrHOTE{bIdx}.dc1.(DFCLab)-pooledCorrHOTE{bIdx}.dc0.(DFCLab),'omitnan')/sqrt(sum(~isnan(pooledCorrHOTE{bIdx}.(dcLab).(DFCLab))));
        pooledSemXCov.diff.(DFCLab)(bIdx,:) = std(pooledCorrXCov{bIdx}.dc1.(DFCLab)-pooledCorrXCov{bIdx}.dc0.(DFCLab),'omitnan')/sqrt(sum(~isnan(pooledCorrXCov{bIdx}.(dcLab).(DFCLab))));
    end
end
%% Lineplot DSC and DFC (for delay consistency = 1) (Fig.5A)
rng(6)
new_bin_minutes = 10; %
repIdx = 1;
bIdx = find(new_bin_minutes_range==new_bin_minutes);
% Take uncorrelated points
t_idxs = 1:new_bin_minutes:size(match_TE_gtS{repIdx}.dc1{bIdx},2);
%t_idxs = 1:1:size(match_TE_gtS{repIdx}.dc1{bIdx},2);

while 1
    rIdx = randperm(size(Comm_pairs{repIdx},1));
    rIdx = rIdx(1:3);
    commIdxs = Comm_pairs{repIdx}(rIdx,:);
    if sum(commIdxs(:,1)>81)==0 % avoid picking inhibitory synapses
        break;
    end
end

for pIdx = 1:numel(rIdx)
    tmpi = find(connPairs{repIdx}.TE(:,1)==commIdxs(pIdx,1));
    tmpj = find(connPairs{repIdx}.TE(:,2)==commIdxs(pIdx,2));
    TEidxs(pIdx) = intersect(tmpi,tmpj);
    
    tmpi = find(connPairs{repIdx}.HOTE(:,1)==commIdxs(pIdx,1));
    tmpj = find(connPairs{repIdx}.HOTE(:,2)==commIdxs(pIdx,2));
    HOTEidxs(pIdx) = intersect(tmpi,tmpj);
    
    tmpi = find(connPairs{repIdx}.XCov(:,1)==commIdxs(pIdx,1));
    tmpj = find(connPairs{repIdx}.XCov(:,2)==commIdxs(pIdx,2));
    XCovidxs(pIdx) = intersect(tmpi,tmpj);
end

% Convert pairs global indeces to indeces of the matched GT pairs
TEidxs=find(ismember(matched_pairs{repIdx}.dc1.TE,TEidxs)==1);
HOTEidxs=find(ismember(matched_pairs{repIdx}.dc1.HOTE,HOTEidxs)==1);
XCovidxs=find(ismember(matched_pairs{repIdx}.dc1.XCov,XCovidxs)==1);

xaxis = t_idxs+new_bin_minutes;
figure()
subplot(2,2,1)
h=plot(xaxis,match_TE_gtS{repIdx}.dc1{bIdx}(TEidxs,t_idxs)','linewidth',2);
legend('Syn. 1','Syn. 2','Syn. 3')
xlim([xaxis(1),xaxis(end)])
xlabel('Time [min]')
ylabel('a.u.')
title('Ground truth')
subplot(2,2,2)
plot(xaxis,match_TE{repIdx}.dc1{bIdx}(TEidxs,t_idxs)','linewidth',2)
xlim([xaxis(1),xaxis(end)])
xlabel('Time [min]')
ylabel('bits')
title('TE')
subplot(2,2,3)
plot(xaxis,match_HOTE{repIdx}.dc1{bIdx}(HOTEidxs,t_idxs)','linewidth',2)
xlim([xaxis(1),xaxis(end)])
xlabel('Time [min]')
ylabel('bits')
title('HOTE')
subplot(2,2,4)
plot(xaxis,match_XCov{repIdx}.dc1{bIdx}(XCovidxs,t_idxs)','linewidth',2)
xlim([xaxis(1),xaxis(end)])
xlabel('Time [min]')
ylabel('corr.')
title('XCov')
sgtitle(['Dynamic connectivity for example synapses'])

cols=[h(1).Color;h(2).Color;h(3).Color];
%% Plot correlation histograms (Fig.5B)
% Compute correlation between DFC and DSC
show_hist_minutes = [10]; % 
bIdx = find(new_bin_minutes_range==show_hist_minutes);
TJ_idx = find(time_jumps==show_hist_minutes);
nPlots = 1+doHO+doXCov+doXCorr;

exampleSyn.HOTE = pairCorrHOTE{repIdx}.dc1.varyTJ{bIdx}(HOTEidxs);
exampleSyn.TE = pairCorrTE{repIdx}.dc1.varyTJ{bIdx}(TEidxs,TJ_idx);
exampleSyn.XCov = pairCorrXCov{repIdx}.dc1.varyTJ{bIdx}(XCovidxs,TJ_idx);

minX = min([pooledCorrTE{bIdx}.dc1.varyTJ(:,TJ_idx);pooledCorrTE{bIdx}.dc1.varyTJ(:,TJ_idx)]);
maxX = max([pooledCorrTE{bIdx}.dc1.varyTJ(:,TJ_idx);pooledCorrTE{bIdx}.dc1.varyTJ(:,TJ_idx)]);
edges = linspace(minX,maxX,20);

for bIdx = 1:numel(show_hist_minutes)
    new_bin_minutes = show_hist_minutes(bIdx);
    bIdx = find(new_bin_minutes_range == new_bin_minutes);

    figure()
    subplot(1,nPlots,1)
    hold on
    h=histogram(pooledCorrTE{bIdx}.dc1.varyTJ(:,TJ_idx),edges,'normalization','probability');
    xlabel('correlation')
    %title(['corr(DSC, DFC_{TE}) over time, time window ',num2str(new_bin_minutes),'min'])
    title('TE')
    xlim([minX,maxX])
    ylim([0,0.15])
    text(h.BinEdges(2),max(h.Values)-0.03,['mean=',num2str(nanmean(pooledCorrTE{bIdx}.dc1.varyTJ(:,TJ_idx)))],'FontSize',18);
    bin=[];
    for idx=1:numel(exampleSyn.TE)
        tmp=find(exampleSyn.TE(idx)>h.BinEdges);
        bin(idx)=tmp(end);
        xval=0.5*(h.BinEdges(bin(idx))+h.BinEdges(bin(idx)+1));
        yval = (h.Values(bin(idx)))-(sum(bin==tmp(end))-1)*0.01;
        scatter(xval,yval,[],cols(idx,:),'filled');
    end
    
    if doHO
        subplot(1,nPlots,2)
        hold on
        h1=histogram(pooledCorrHOTE{bIdx}.dc1.varyTJ(:,TJ_idx),edges,'normalization','probability');
        xlabel('correlation')
        title('HOTE')
        xlim([minX,maxX])
        ylim([0,0.15])
        %title(['corr(DSC, DFC_{HOTE}) over time, time window ',num2str(new_bin_minutes),'min'])
        text(h1.BinEdges(2),max(h1.Values)-0.03,['mean=',num2str(nanmean(pooledCorrHOTE{bIdx}.dc1.varyTJ(:,TJ_idx)))],'FontSize',18);
        bin=[];
        for idx=1:numel(exampleSyn.TE)
            tmp=find(exampleSyn.HOTE(idx)>h1.BinEdges);
            bin(idx)=tmp(end);
            xval=0.5*(h1.BinEdges(bin(idx))+h1.BinEdges(bin(idx)+1));
            yval = (h1.Values(bin(idx)))-(sum(bin==tmp(end))-1)*0.01;
            scatter(xval,yval,[],cols(idx,:),'filled');
        end
    end
    if doXCov
        subplot(1,nPlots,3)
        hold on
        h2=histogram(pooledCorrXCov{bIdx}.dc1.varyTJ(:,TJ_idx),edges,'normalization','probability');
        xlabel('correlation')
        %title(['corr(DSC, DFC_{XCov}) over time, time window ',num2str(new_bin_minutes),'min'])
        title('XCov')
        xlim([minX,maxX])
        text(h2.BinEdges(2),max(h2.Values)-0.04,['mean=',num2str(nanmean(pooledCorrXCov{bIdx}.dc1.varyTJ(:,TJ_idx)))],'FontSize',18);
        bin=[];
        for idx=1:numel(exampleSyn.TE)
            tmp=find(exampleSyn.XCov(idx)>h2.BinEdges);
            bin(idx)=tmp(end);
            xval=0.5*(h2.BinEdges(bin(idx))+h2.BinEdges(bin(idx)+1));
            yval = (h2.Values(bin(idx)))-(sum(bin==tmp(end))-1)*0.01;
            scatter(xval,yval,[],cols(idx,:),'filled');
        end
    end
end
%% Plot DSC and DFC correlation with sliding window size (Fig.5C)
cols = [0.4660 0.6740 0.1880;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0 0.4470 0.7410;0.8500 0.3250 0.0980];
%corrType = {'min1','varyTJ','NsampleMatch'};
corrType = {'NsampleMatch'};

%cols = [1 0 0; 0 1 0; 0 0 1];

for dcIdx=1:2
    dcLab=dc_label{dcIdx};
    for idx=1:numel(corrType)
        figure()
        % on the diagonal we have T = TJ
        l1=plot(new_bin_minutes_range(1:size(pooledMeanTE.(dcLab).(corrType{idx}),1)),diag(pooledMeanTE.(dcLab).(corrType{idx})),'color',cols(1,:),'linewidth',2);
        hold on
        l2=plot(new_bin_minutes_range(1:size(pooledMeanXCov.(dcLab).(corrType{idx}),1)),diag(pooledMeanXCov.(dcLab).(corrType{idx})),'color',cols(3,:),'linewidth',2);
        l3=plot(new_bin_minutes_range(1:size(pooledMeanHOTE.(dcLab).(corrType{idx}),1)),diag(pooledMeanHOTE.(dcLab).(corrType{idx})),'color',cols(2,:),'linewidth',2);
        legend('TE','XCov','HOTE','Autoupdate','off')
        shadedErrorBar(new_bin_minutes_range(1:size(pooledMeanTE.(dcLab).(corrType{idx}),1)),diag(pooledMeanTE.(dcLab).(corrType{idx})),diag(pooledSemTE.(dcLab).(corrType{idx})),'lineProps',{'color',cols(1,:)})
        shadedErrorBar(new_bin_minutes_range(1:size(pooledMeanXCov.(dcLab).(corrType{idx}),1)),diag(pooledMeanXCov.(dcLab).(corrType{idx})),diag(pooledSemXCov.(dcLab).(corrType{idx})),'lineProps',{'color',cols(3,:)})
        shadedErrorBar(new_bin_minutes_range(1:size(pooledMeanHOTE.(dcLab).(corrType{idx}),1)),diag(pooledMeanHOTE.(dcLab).(corrType{idx})),diag(pooledSemHOTE.(dcLab).(corrType{idx})),'lineProps',{'color',cols(2,:)})
        xlabel('window size [min]')
        ylabel('correlation')
        xlim([min(new_bin_minutes_range(1:size(pooledMeanTE.(dcLab).(corrType{idx}),1))) max(new_bin_minutes_range(1:size(pooledMeanTE.(dcLab).(corrType{idx}),1)))])
        ylim([0.05 0.85])
        title(['Mean correlation between DFC and DSC, ',corrType{idx},'; delay cons.=',dcLab(3)])
    end

end
%

for idx=1:numel(corrType)
    figure()
    % on the diagonal we have T = TJ
    pooledMeanTE.diff.(corrType{idx})=diag(pooledMeanTE.dc1.(corrType{idx}))-diag(pooledMeanTE.dc0.(corrType{idx}));
    l1=plot(new_bin_minutes_range(1:size(pooledMeanTE.(dcLab).(corrType{idx}),1)),pooledMeanTE.diff.(corrType{idx}),'color',cols(1,:),'linewidth',2);
    hold on
    pooledMeanXCov.diff.(corrType{idx})=diag(pooledMeanXCov.dc1.(corrType{idx}))-diag(pooledMeanXCov.dc0.(corrType{idx}));
    l2=plot(new_bin_minutes_range(1:size(pooledMeanXCov.(dcLab).(corrType{idx}),1)),pooledMeanXCov.diff.(corrType{idx}),'color',cols(3,:),'linewidth',2);
    pooledMeanHOTE.diff.(corrType{idx})=diag(pooledMeanHOTE.dc1.(corrType{idx}))-diag(pooledMeanHOTE.dc0.(corrType{idx}));
    l3=plot(new_bin_minutes_range(1:size(pooledMeanHOTE.(dcLab).(corrType{idx}),1)),pooledMeanHOTE.diff.(corrType{idx}),'color',cols(2,:),'linewidth',2);
    legend('TE','XCov','HOTE','Autoupdate','off')
    shadedErrorBar(new_bin_minutes_range(1:size(pooledMeanTE.diff.(corrType{idx}),1)),pooledMeanTE.diff.(corrType{idx}),pooledSemTE.diff.(corrType{idx}),'lineProps',{'color',cols(1,:)})
    shadedErrorBar(new_bin_minutes_range(1:size(pooledMeanXCov.diff.(corrType{idx}),1)),pooledMeanXCov.diff.(corrType{idx}),pooledSemTE.diff.(corrType{idx}),'lineProps',{'color',cols(3,:)})
    shadedErrorBar(new_bin_minutes_range(1:size(pooledMeanHOTE.diff.(corrType{idx}),1)),pooledMeanHOTE.diff.(corrType{idx}),pooledSemTE.diff.(corrType{idx}),'lineProps',{'color',cols(2,:)})
    xlabel('window size [min]')
    ylabel('correlation')
    xlim([min(new_bin_minutes_range(1:size(pooledMeanTE.(dcLab).(corrType{idx}),1))) max(new_bin_minutes_range(1:size(pooledMeanTE.(dcLab).(corrType{idx}),1)))])
    title(['Mean correlation gain from delay consistency, ',corrType{idx}])
end
