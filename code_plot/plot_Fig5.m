%% Static connectivity plot
% For BI resubmision, adjust XCov pairs by taking abs value before
% selecting them

clear all; clc;
%rng(4)
norm_method='normalized'; % 'none' or 'normalized'
new_bin_minutes_range = [1 2 5 10 20 30];
time_jumps = [1, 2, 5, 10, 20, 30];
doHO = 1; doXCov = 1; doXCorr = 1;
% Load simulation results
if strcmp(norm_method,'none')
    % delat consistency
    load('D:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\full_seed1_IzhiModel_v10_90min_noXCorr_varyWindow_23-Feb-2022','pairCorrTE','pairCorrHOTE','match_TE_gtS_time','match_TE_time','match_HOTE_gtS_time','match_HOTE_time','doHO','x_axis','x_labels','time_jump_minutes','meanCorrTE','meanCorrHOTE','match_pairs','TE_connPairs','HOTE_connPairs')
    load('D:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\full_seed1_IzhiModel_v10_90min_XCov_varyWindow_24-Feb-2022','pairCorrXCov','match_XCov_gtS_time','match_XCov_time','doXCov','doXCorr','meanCorrXCov','XCov_connPairs')
    % no Delay consistency
    ndc = load('D:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\full_seed1_IzhiModel_v12_90min_varyWindow_noNorm_0_DelConsist_04-Mar-2022','pairCorrTE','pairCorrHOTE','new_bin_minutes_range','match_TE_gtS_time','match_TE_time','match_HOTE_gtS_time','match_HOTE_time','doHO','x_axis','x_labels','time_jump_minutes','meanCorrTE','meanCorrHOTE','pairCorrXCov','new_bin_minutes_range','match_XCov_gtS_time','match_XCov_time','doXCov','doXCorr','meanCorrXCov');
elseif strcmp(norm_method,'normalized')
    % Load consistent delay XCov and TE
    load('E:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\BI_resubmission\full_seed1_IzhiModel_v15_180min_varyWindow_1_DelConsist_15-Sep-2022','pairCorrTE','match_TE_gtS_time','match_TE_time','pairCorrHOTE','match_HOTE_gtS_time','match_HOTE_time','pairCorrXCov','match_XCov_gtS_time','match_XCov_time','doXCov','doXCorr','meanCorrXCov','XCov_connPairs','x_axis','x_labels','time_jump_minutes','meanCorrTE','match_pairs','TE_connPairs','meanCorrHOTE','HOTE_connPairs');
%     % Load consistent delay XCov
%     load('D:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\full_seed1_IzhiModel_v12_normalized_90min_varyWindow_1_DelConsist_08-Mar-2022.mat','pairCorrXCov','match_XCov_gtS_time','match_XCov_time','doXCov','doXCorr','meanCorrXCov','XCov_connPairs')
%     Load Non-consistent delay TE and HOTE
%     ndc = load('D:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\full_seed1_IzhiModel_v12_90min_varyWindow_noNorm_0_DelConsist_04-Mar-2022','pairCorrTE','pairCorrHOTE','new_bin_minutes_range','match_TE_gtS_time','match_TE_time','match_HOTE_gtS_time','match_HOTE_time','doHO','x_axis','x_labels','time_jump_minutes','meanCorrTE','meanCorrHOTE','new_bin_minutes_range','match_XCov_time','doXCorr','meanCorrXCov');
%     Load Non-consistent delay XCov
    ndc = load('E:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\BI_resubmission\full_seed1_IzhiModel_v15_180min_varyWindow_0_DelConsist_17-Sep-2022','pairCorrTE','match_TE_gtS_time','match_TE_time','pairCorrHOTE','match_HOTE_gtS_time','match_HOTE_time','pairCorrXCov','match_XCov_gtS_time','match_XCov_time','doXCov','doXCorr','meanCorrXCov','XCov_connPairs','x_axis','x_labels','time_jump_minutes','meanCorrTE','match_pairs','TE_connPairs','meanCorrHOTE','HOTE_connPairs');
%     newfields = fieldnames(ndc2);
%     for fIdx=1:numel(newfields)
%         ndc.(newfields{fIdx})=ndc2.(newfields{fIdx});
%     end
end

% XCov selected network
% In te script testModel I accidentally forgot to take the absolute value
% of XCov before selecting the synapses on which I computed DFC. Therefore
% I selected all excitatory synapses for XCov, while I neglected inhibitory
% synapses (N=200). Therefore, to keep 300 excitatory synapses as for the
% other measues I delete here 200 random excitatory synapses from XCov.
rng(1)
rIdx = randperm(size(XCov_connPairs,1));
rIdx = rIdx(1:300);
XCov_connPairs = XCov_connPairs(rIdx,:);
ndc.XCov_connPairs = ndc.XCov_connPairs(rIdx,:);

tmpPairCorrTE.dc1 = pairCorrTE;
tmpPairCorrHOTE.dc1 = pairCorrHOTE;
tmpPairCorrXCov.dc1 = pairCorrXCov;
tmpMeanCorrTE.dc1 = meanCorrTE;
tmpMeanCorrHOTE.dc1 = meanCorrHOTE;
tmpMeanCorrXCov.dc1 = meanCorrXCov;

tmpPairCorrTE.dc0 = ndc.pairCorrTE;
tmpPairCorrHOTE.dc0 = ndc.pairCorrHOTE;
tmpPairCorrXCov.dc0 = ndc.pairCorrXCov;
tmpMeanCorrTE.dc0 = ndc.meanCorrTE;
tmpMeanCorrHOTE.dc0 = ndc.meanCorrHOTE;
tmpMeanCorrXCov.dc0 = ndc.meanCorrXCov;
clear pairCorrTE pairCorrHOTE pairCorrXCov meanCorrTE meanCorrHOTE meanCorrXCov

for bIdx = 1:numel(new_bin_minutes_range)
    pairCorrTE.dc0.minTJ{bIdx}=tmpPairCorrTE.dc0{bIdx};
    pairCorrHOTE.dc0.minTJ{bIdx}=tmpPairCorrHOTE.dc0{bIdx};
    pairCorrXCov.dc0.minTJ{bIdx}=tmpPairCorrXCov.dc0{bIdx}(rIdx);
    
    pairCorrTE.dc1.minTJ{bIdx}=tmpPairCorrTE.dc1{bIdx};
    pairCorrHOTE.dc1.minTJ{bIdx}=tmpPairCorrHOTE.dc1{bIdx};
    pairCorrXCov.dc1.minTJ{bIdx}=tmpPairCorrXCov.dc1{bIdx}(rIdx);
end
meanCorrTE.dc0.minTJ=tmpMeanCorrTE.dc0;
meanCorrHOTE.dc0.minTJ=tmpMeanCorrHOTE.dc0;
meanCorrXCov.dc0.minTJ=tmpMeanCorrXCov.dc0;

meanCorrTE.dc1.minTJ=tmpMeanCorrTE.dc1(1:numel(new_bin_minutes_range));
meanCorrHOTE.dc1.minTJ=tmpMeanCorrHOTE.dc1(1:numel(new_bin_minutes_range));
meanCorrXCov.dc1.minTJ=tmpMeanCorrXCov.dc1(1:numel(new_bin_minutes_range));

% Initialize meatched_connections
for bIdx = 1:numel(new_bin_minutes_range)
    match_XCov.dc1{bIdx} = match_XCov_time{bIdx}(rIdx,:);
    match_XCov_gtS.dc1{bIdx} = match_XCov_gtS_time{bIdx}(rIdx,:);
    match_XCov.dc0{bIdx} = ndc.match_XCov_time{bIdx}(rIdx,:);
    match_XCov_gtS.dc0{bIdx} = ndc.match_XCov_gtS_time{bIdx}(rIdx,:);
end
match_TE.dc1 = match_TE_time;
match_TE_gtS.dc1 = match_TE_gtS_time;
match_HOTE.dc1 = match_HOTE_time;
match_HOTE_gtS.dc1 = match_HOTE_gtS_time;
match_TE.dc0 = ndc.match_TE_time;
match_TE_gtS.dc0 = ndc.match_TE_gtS_time;
match_HOTE.dc0 = ndc.match_HOTE_time;
match_HOTE_gtS.dc0 = ndc.match_HOTE_gtS_time;

dc_label = {'dc0','dc1'};

clear ndc

% Common pairs to all measures
tmpComm_pairs = intersect(TE_connPairs,XCov_connPairs,'rows');
Comm_pairs = intersect(tmpComm_pairs,HOTE_connPairs,'rows');
% Comm_pairs = intersect(TE_connPairs,XCov_connPairs,'rows');
%% Compute non-autocorr weights
doHO = 1; doXCov = 1; doXCorr = 0;
for dcIdx=1:2
    dcLab=dc_label{dcIdx};
    for bIdx = 1:numel(new_bin_minutes_range)
        new_bin_minutes = new_bin_minutes_range(bIdx);
        
        for TJ_idx = 1:numel(time_jumps)
            TJ = time_jumps(TJ_idx);
            t_idxs = 1:(TJ)/(time_jumps(1)):size(match_TE_gtS.(dcLab){bIdx},2);

            for pairIdx = 1:numel(pairCorrTE.(dcLab).minTJ{bIdx})  % average correlation on true positives
                tmp = corrcoef(match_TE_gtS.(dcLab){bIdx}(pairIdx,t_idxs),match_TE.(dcLab){bIdx}(pairIdx,t_idxs));
                pairCorrTE.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                tmp = corrcoef(diff(match_TE_gtS.(dcLab){bIdx}(pairIdx,t_idxs)),diff(match_TE.(dcLab){bIdx}(pairIdx,t_idxs)));
                pairCorrTEdiff.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
            end
            meanCorrTE.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrTE.(dcLab).varyTJ{bIdx}(:,TJ_idx));
            meanCorrTEdiff.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrTEdiff.(dcLab).varyTJ{bIdx}(:,TJ_idx));
            semCorrTE.(dcLab).varyTJ(bIdx,TJ_idx) = std(pairCorrTE.(dcLab).varyTJ{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrTE.(dcLab).varyTJ{bIdx}(:,TJ_idx))));
            %semCorrTE.(dcLab).minTJ(bIdx) = std(pairCorrTE.(dcLab).minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrTE.(dcLab).minTJ{bIdx})));

            if doHO
                for pairIdx = 1:numel(pairCorrHOTE.(dcLab).minTJ{bIdx})
                    tmp = corrcoef(match_HOTE_gtS.(dcLab){bIdx}(pairIdx,t_idxs),match_HOTE.(dcLab){bIdx}(pairIdx,t_idxs));
                    pairCorrHOTE.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    tmp = corrcoef(diff(match_HOTE_gtS.(dcLab){bIdx}(pairIdx,t_idxs)),diff(match_HOTE.(dcLab){bIdx}(pairIdx,t_idxs)));
                    pairCorrHOTEdiff.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                end
                meanCorrHOTE.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrHOTE.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                meanCorrHOTEdiff.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrHOTEdiff.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                semCorrHOTE.(dcLab).varyTJ(bIdx,TJ_idx) = std(pairCorrHOTE.(dcLab).varyTJ{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrHOTE.(dcLab).varyTJ{bIdx}(:,TJ_idx))));
                %semCorrHOTE.(dcLab).minTJ(bIdx) = std(pairCorrHOTE.(dcLab).minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrHOTE.(dcLab).minTJ{bIdx})));
            end
            if doXCov
                for pairIdx = 1:numel(pairCorrXCov.(dcLab).minTJ{bIdx})
                    tmp = corrcoef(match_XCov_gtS.(dcLab){bIdx}(pairIdx,t_idxs),match_XCov.(dcLab){bIdx}(pairIdx,t_idxs));
                    pairCorrXCov.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    tmp = corrcoef(diff(match_XCov_gtS.(dcLab){bIdx}(pairIdx,t_idxs)),diff(match_XCov.(dcLab){bIdx}(pairIdx,t_idxs)));
                    pairCorrXCovdiff.(dcLab).varyTJ{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                end
                meanCorrXCov.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrXCov.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                meanCorrXCovdiff.(dcLab).varyTJ(bIdx,TJ_idx) = nanmean(pairCorrXCovdiff.(dcLab).varyTJ{bIdx}(:,TJ_idx));
                semCorrXCov.(dcLab).varyTJ(bIdx,TJ_idx) = std(pairCorrXCov.(dcLab).varyTJ{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrXCov.(dcLab).varyTJ{bIdx}(:,TJ_idx))));
                %semCorrXCov.(dcLab).minTJ(bIdx) = std(pairCorrXCov.(dcLab).minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrXCov.(dcLab).minTJ{bIdx})));
            end
        end
    end
end

%% Compute same n samples weights
%load('E:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\full_seed1_IzhiModel_v15_90min_3sTJ_varyWindow_1_DelConsist_06-Sep-2022','post','s_time')
load('E:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\BI_resubmission\full_seed1_IzhiModel_v15_180min_varyWindow_1_DelConsist_15-Sep-2022','post','s_time')

max_minutes = 30;
sIdx = find(new_bin_minutes_range==max_minutes);
new_bin_minutes = new_bin_minutes_range(sIdx);
maxTJ = time_jumps(end);
n_samples = numel(1:maxTJ/(time_jumps(1)):size(match_TE_gtS.(dcLab){sIdx},2)); % set size of 10 no autocorrelated

for dcIdx=1:2
    dcLab=dc_label{dcIdx};
    for bIdx = 1:sIdx
        new_bin_minutes = new_bin_minutes_range(bIdx);
        % Take time points so that they're uncorrelated
        for TJ_idx = 1:numel(time_jumps)
            TJ = time_jumps(TJ_idx);
            t_idxs = 1:(TJ)/(time_jumps(1)):size(match_TE_gtS.(dcLab){bIdx},2);

            nDiff = numel(t_idxs)-n_samples;
            t_idxs(1:floor(nDiff/2))=[];
            t_idxs(end-ceil(nDiff/2)+1:end)=[];

            for idx = 1:size(TE_connPairs,1)
                ridx = find(post(TE_connPairs(idx,1),:)==TE_connPairs(idx,2));
                s_time_TE(idx,:) = s_time(TE_connPairs(idx,1),ridx,:);
            end
            s_time_dawnsamp = temporal_rebinning(s_time_TE,TJ*60,'movmean',TJ*60);
            nDiff = size(s_time_dawnsamp,2)-n_samples;
            t_idxs_gt = 1:size(s_time_dawnsamp,2);
            t_idxs_gt(1:floor(nDiff/2))=[];
            t_idxs_gt(end-ceil(nDiff/2)+1:end)=[];
            
            for pairIdx = 1:numel(pairCorrTE.(dcLab).minTJ{bIdx})  % average correlation on true positives
                %tmp = corrcoef(match_TE_gtS.(dcLab){bIdx}(pairIdx,t_idxs),match_TE.(dcLab){bIdx}(pairIdx,t_idxs));
                tmp = corrcoef(s_time_dawnsamp(pairIdx,t_idxs_gt),match_TE.(dcLab){bIdx}(pairIdx,t_idxs));
                pairCorrTE.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                meanCorrTE.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrTE.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                semCorrTE.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrTE.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrTE.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                tmp = corrcoef(diff(s_time_dawnsamp(pairIdx,t_idxs_gt)),diff(match_TE.(dcLab){bIdx}(pairIdx,t_idxs)));
                pairCorrTEdiff.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                meanCorrTEdiff.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrTEdiff.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                semCorrTEdiff.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrTEdiff.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrTEdiff.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
            end
            if doHO
                for idx = 1:size(HOTE_connPairs,1)
                    ridx = find(post(HOTE_connPairs(idx,1),:)==HOTE_connPairs(idx,2));
                    s_time_HOTE(idx,:) = s_time(HOTE_connPairs(idx,1),ridx,:);
                end
                s_time_dawnsamp = temporal_rebinning(s_time_HOTE,TJ*60,'movmean',TJ*60);
                nDiff = size(s_time_dawnsamp,2)-n_samples;
                t_idxs_gt = 1:size(s_time_dawnsamp,2);
                t_idxs_gt(1:floor(nDiff/2))=[];
                t_idxs_gt(end-ceil(nDiff/2)+1:end)=[];
            
                for pairIdx = 1:numel(pairCorrHOTE.(dcLab).minTJ{bIdx})
                    tmp = corrcoef(s_time_dawnsamp(pairIdx,t_idxs_gt),match_HOTE.(dcLab){bIdx}(pairIdx,t_idxs));
                    pairCorrHOTE.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    meanCorrHOTE.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrHOTE.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                    semCorrHOTE.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrHOTE.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrHOTE.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                    tmp = corrcoef(diff(s_time_dawnsamp(pairIdx,t_idxs_gt)),diff(match_HOTE.(dcLab){bIdx}(pairIdx,t_idxs)));
                    pairCorrHOTEdiff.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    meanCorrHOTEdiff.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrHOTEdiff.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                    semCorrHOTEdiff.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrHOTEdiff.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrHOTEdiff.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                end
            end
            if doXCov
                for idx = 1:size(XCov_connPairs,1)
                    ridx = find(post(XCov_connPairs(idx,1),:)==XCov_connPairs(idx,2));
                    s_time_XCov(idx,:) = s_time(XCov_connPairs(idx,1),ridx,:);
                end
                s_time_dawnsamp = temporal_rebinning(s_time_XCov,TJ*60,'movmean',TJ*60);
                nDiff = size(s_time_dawnsamp,2)-n_samples;
                t_idxs_gt = 1:size(s_time_dawnsamp,2);
                t_idxs_gt(1:floor(nDiff/2))=[];
                t_idxs_gt(end-ceil(nDiff/2)+1:end)=[];
                for pairIdx = 1:numel(pairCorrXCov.(dcLab).minTJ{bIdx})
                    tmp = corrcoef(s_time_dawnsamp(pairIdx,t_idxs_gt),match_XCov.(dcLab){bIdx}(pairIdx,t_idxs));
                    pairCorrXCov.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    meanCorrXCov.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrXCov.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                    semCorrXCov.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrXCov.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrXCov.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                    tmp = corrcoef(diff(s_time_dawnsamp(pairIdx,t_idxs_gt)),diff(match_XCov.(dcLab){bIdx}(pairIdx,t_idxs)));
                    pairCorrXCovdiff.(dcLab).NsampleMatch{bIdx}(pairIdx,TJ_idx) = tmp(1,2);
                    meanCorrXCovdiff.(dcLab).NsampleMatch(bIdx,TJ_idx) = nanmean(pairCorrXCovdiff.(dcLab).NsampleMatch{bIdx}(:,TJ_idx));
                    semCorrXCovdiff.(dcLab).NsampleMatch(bIdx,TJ_idx) = std(pairCorrXCovdiff.(dcLab).NsampleMatch{bIdx}(:,TJ_idx),'omitnan')/sqrt(sum(~isnan(pairCorrXCovdiff.(dcLab).NsampleMatch{bIdx}(:,TJ_idx))));
                end
            end
        end
    end
end

% Compute SEM for difference of dc1 and dc0 (the pairs are matched since
% the delay consistency only affects the DFC and not the DSC)
for bIdx = 1:numel(new_bin_minutes_range)
    semCorrTE.diff.minTJ(bIdx) = std(pairCorrTE.dc1.minTJ{bIdx}-pairCorrTE.dc0.minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrTE.(dcLab).minTJ{bIdx})));
    semCorrTE.diff.varyTJ(bIdx) = std(pairCorrTE.dc1.varyTJ{bIdx}-pairCorrTE.dc0.varyTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrTE.(dcLab).varyTJ{bIdx})));
    semCorrHOTE.diff.minTJ(bIdx) = std(pairCorrHOTE.dc1.minTJ{bIdx}-pairCorrHOTE.dc0.minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrHOTE.(dcLab).minTJ{bIdx})));
    semCorrHOTE.diff.varyTJ(bIdx) = std(pairCorrHOTE.dc1.varyTJ{bIdx}-pairCorrHOTE.dc0.varyTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrHOTE.(dcLab).varyTJ{bIdx})));
    semCorrXCov.diff.minTJ(bIdx) = std(pairCorrXCov.dc1.minTJ{bIdx}-pairCorrXCov.dc0.minTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrXCov.(dcLab).minTJ{bIdx})));
    semCorrXCov.diff.varyTJ(bIdx) = std(pairCorrXCov.dc1.varyTJ{bIdx}-pairCorrXCov.dc0.varyTJ{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrXCov.(dcLab).varyTJ{bIdx})));
end

for bIdx = 1:sIdx
    semCorrTE.diff.NsampleMatch(bIdx) = std(pairCorrTE.dc1.NsampleMatch{bIdx}-pairCorrTE.dc0.NsampleMatch{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrTE.(dcLab).NsampleMatch{bIdx})));
    semCorrHOTE.diff.NsampleMatch(bIdx) = std(pairCorrHOTE.dc1.NsampleMatch{bIdx}-pairCorrHOTE.dc0.NsampleMatch{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrHOTE.(dcLab).NsampleMatch{bIdx})));
    semCorrXCov.diff.NsampleMatch(bIdx) = std(pairCorrXCov.dc1.NsampleMatch{bIdx}-pairCorrXCov.dc0.NsampleMatch{bIdx},'omitnan')/sqrt(sum(~isnan(pairCorrXCov.(dcLab).NsampleMatch{bIdx})));
end

%% Plot DSC-DFC corr trend with TJ < T
figure()
h = axes;
%plot(time_jumps*60,meanCorrXCovdiff.dc1.NsampleMatch')
for idx = 1:size(meanCorrXCovdiff.dc1.NsampleMatch,1)
    hold on
    errorbar(time_jumps*60,meanCorrXCovdiff.dc1.NsampleMatch(idx,:)',semCorrXCovdiff.dc1.NsampleMatch(idx,:)')
end
xlabel('log(time precision [s]')
ylabel('\DeltaW corr.')
title('\DeltaW correlation vs temporal resolution for different T')
legend('T = 1min','T = 2min', 'T = 5min','AutoUpdate','off')
xline(60,'k--')

h.XScale='log';


figure()
h = axes;
%plot(time_jumps*60,meanCorrXCovdiff.dc1.NsampleMatch')
for idx = 1:size(meanCorrXCov.dc1.NsampleMatch,1)
    hold on
    errorbar(time_jumps*60,meanCorrXCov.dc1.NsampleMatch(idx,:)',semCorrXCov.dc1.NsampleMatch(idx,:)')
end
xlabel('log(time precision [s]')
ylabel('W corr.')
title('Weights correlation vs temporal resolution for different T')
legend('T = 1min','T = 2min', 'T = 5min','AutoUpdate','off')
xline(60,'k--')

h.XScale='log';
%% Plot DSC vs DFC
rIdx = randperm(200);
rIdx = rIdx(1:20);
dcLab = 'dc1';

for bIdx = [1,3,numel(new_bin_minutes_range)]
    new_bin_minutes = new_bin_minutes_range(bIdx);
    % Take time points so that they're uncorrelated
    %t_idxs = 1:new_bin_minutes:size(match_TE_gtS{bIdx},2);
    % Take all points
    t_idxs = 1:size(match_TE_gtS.(dcLab){bIdx},2);
    % Plot TE DFC
    figure()
    subplot(1,2,1)
    imagesc(match_TE_gtS.(dcLab){bIdx}(1:20,t_idxs))
    colorbar()
    set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
    xlabel('minutes')
    ylabel('pairs')
    title('Ground truth')

    subplot(1,2,2)
    imagesc(match_TE.(dcLab){bIdx}(1:20,t_idxs))
    colorbar()
    set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
    xlabel('minutes')
    ylabel('pairs')
    title('Measured (TE)')
    if bIdx == 1
        tmp = caxis;
            cmax=tmp(2);
            caxis=[tmp(1),cmax*(1/2)];
    end
    %sgtitle(['TE synaptic weights over time, ',num2str(new_bin_minutes),'min window, ',num2str(time_jump_minutes),'min time jump']);
    sgtitle(['TE synaptic weights over time, ',num2str(new_bin_minutes),'min window']);

    %Plot HOTE DFC
    if doHO
        figure()
        subplot(1,2,1)
        imagesc(match_HOTE_gtS.(dcLab){bIdx}(1:20,t_idxs))
        colorbar()
        set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
        xlabel('minutes')
        ylabel('pairs')
        title('Ground truth')

        subplot(1,2,2)
        imagesc(match_HOTE.(dcLab){bIdx}(1:20,t_idxs))
        colorbar()
        set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
        xlabel('minutes')
        ylabel('pairs')
        title('Measured (HOTE)')
        if bIdx == 1
            tmp = caxis;
            cmax=tmp(2);
            caxis=[tmp(1),cmax*(1/2)];
        end
        %sgtitle(['HOTE synaptic weights over time, ',num2str(new_bin_minutes),'min window, ',num2str(time_jump_minutes),'min time jump']);
        sgtitle(['HOTE synaptic weights over time, ',num2str(new_bin_minutes),'min window']);
    end

    % Plot XCov GT
    if doXCov
        figure()
        subplot(1,2,1)
        imagesc(match_XCov_gtS.(dcLab){bIdx}(1:20,t_idxs))
        colorbar()
        set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
        xlabel('minutes')
        ylabel('pairs')
        title('Ground truth')

        subplot(1,2,2)
        imagesc(match_XCov.(dcLab){bIdx}(1:20,t_idxs))
        colorbar()
        set(gca,'XTick',x_axis,'XTickLabels',x_labels)   % X Label Ticks
        xlabel('minutes')
        ylabel('pairs')
        title('Measured (XCov)')
        if bIdx == 1
            tmp = caxis;
            cmax=tmp(2);
            caxis=[tmp(1),cmax*(1/2)];
        end
        %sgtitle(['XCov synaptic weights over time, ',num2str(new_bin_minutes),'min window, ',num2str(time_jump_minutes),'min time jump']);
        sgtitle(['XCov synaptic weights over time, ',num2str(new_bin_minutes),'min window']);
    end
end

%% Lineplot DSC and DFC (for delay consistency = 1)
rng(6)
new_bin_minutes = 10;
bIdx = find(new_bin_minutes_range==new_bin_minutes);
% Take uncorrelated points
t_idxs = 1:new_bin_minutes:size(match_TE_gtS.dc1{bIdx},2);
%t_idxs = 1:1:size(match_TE_gtS.dc1{bIdx},2);

while 1
    rIdx = randperm(size(Comm_pairs,1));
    rIdx = rIdx(1:3);
    commIdxs = Comm_pairs(rIdx,:);
    if sum(commIdxs(:,1)>81)==0
        break;
    end
end

for pIdx = 1:numel(rIdx)
    tmpi = find(TE_connPairs(:,1)==commIdxs(pIdx,1));
    tmpj = find(TE_connPairs(:,2)==commIdxs(pIdx,2));
    TEidxs(pIdx) = intersect(tmpi,tmpj);
    
    tmpi = find(HOTE_connPairs(:,1)==commIdxs(pIdx,1));
    tmpj = find(HOTE_connPairs(:,2)==commIdxs(pIdx,2));
    HOTEidxs(pIdx) = intersect(tmpi,tmpj);
    
    tmpi = find(XCov_connPairs(:,1)==commIdxs(pIdx,1));
    tmpj = find(XCov_connPairs(:,2)==commIdxs(pIdx,2));
    XCovidxs(pIdx) = intersect(tmpi,tmpj);
end

xaxis = t_idxs+new_bin_minutes;
figure()
subplot(2,2,1)
h=plot(xaxis,match_TE_gtS.dc1{bIdx}(TEidxs,t_idxs)','linewidth',2);
legend('Syn. 1','Syn. 2','Syn. 3')
xlim([xaxis(1),xaxis(end)])
xlabel('Time [min]')
ylabel('a.u.')
title('Ground truth')
subplot(2,2,2)
plot(xaxis,match_TE.dc1{bIdx}(TEidxs,t_idxs)','linewidth',2)
xlim([xaxis(1),xaxis(end)])
xlabel('Time [min]')
ylabel('bits')
title('TE')
subplot(2,2,3)
plot(xaxis,match_HOTE.dc1{bIdx}(HOTEidxs,t_idxs)','linewidth',2)
xlim([xaxis(1),xaxis(end)])
xlabel('Time [min]')
ylabel('bits')
title('HOTE')
subplot(2,2,4)
plot(xaxis,match_XCov.dc1{bIdx}(XCovidxs,t_idxs)','linewidth',2)
xlim([xaxis(1),xaxis(end)])
xlabel('Time [min]')
ylabel('corr.')
title('XCov')
sgtitle(['Dynamic connectivity for example synapses'])

cols=[h(1).Color;h(2).Color;h(3).Color];
%% Plot correlation histograms
% Compute correlation between DFC and DSC
show_hist_minutes = [10];
bIdx = find(new_bin_minutes_range==show_hist_minutes);
TJ_idx = find(time_jumps==show_hist_minutes);
nPlots = 1+doHO+doXCov+doXCorr;

exampleSyn.HOTE = pairCorrHOTE.dc1.varyTJ{bIdx}(HOTEidxs);
exampleSyn.TE = pairCorrTE.dc1.varyTJ{bIdx}(TEidxs,TJ_idx);
exampleSyn.XCov = pairCorrXCov.dc1.varyTJ{bIdx}(XCovidxs,TJ_idx);

% minX = min([pairCorrTE.dc1.varyTJ{bIdx},pairCorrHOTE.dc1.varyTJ{bIdx},pairCorrXCov.dc1.varyTJ{bIdx}]);
% maxX = max([pairCorrTE.dc1.varyTJ{bIdx},pairCorrHOTE.dc1.varyTJ{bIdx},pairCorrXCov.dc1.varyTJ{bIdx}]);
minX = min([pairCorrTE.dc1.varyTJ{bIdx}(:,TJ_idx);pairCorrXCov.dc1.varyTJ{bIdx}(:,TJ_idx)]);
maxX = max([pairCorrTE.dc1.varyTJ{bIdx}(:,TJ_idx);pairCorrXCov.dc1.varyTJ{bIdx}(:,TJ_idx)]);
edges = linspace(minX,maxX,20);

for bIdx = 1:numel(show_hist_minutes)
    new_bin_minutes = show_hist_minutes(bIdx);
    bIdx = find(new_bin_minutes_range == new_bin_minutes);

    figure()
    subplot(1,nPlots,1)
    hold on
    h=histogram(pairCorrTE.dc1.varyTJ{bIdx},edges,'normalization','probability');
    xlabel('correlation')
    %title(['corr(DSC, DFC_{TE}) over time, time window ',num2str(new_bin_minutes),'min'])
    title('TE')
    xlim([minX,maxX])
    ylim([0,0.15])
    text(h.BinEdges(2),max(h.Values)-0.03,['mean=',num2str(meanCorrTE.dc1.varyTJ(bIdx,TJ_idx),3)],'FontSize',18);
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
        h1=histogram(pairCorrHOTE.dc1.varyTJ{bIdx},edges,'normalization','probability');
        xlabel('correlation')
        title('HOTE')
        xlim([minX,maxX])
        ylim([0,0.15])
        %title(['corr(DSC, DFC_{HOTE}) over time, time window ',num2str(new_bin_minutes),'min'])
        text(h1.BinEdges(2),max(h1.Values)-0.03,['mean=',num2str(meanCorrHOTE.dc1.varyTJ(bIdx,TJ_idx),3)],'FontSize',18);
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
        h2=histogram(pairCorrXCov.dc1.varyTJ{bIdx},edges,'normalization','probability');
        xlabel('correlation')
        %title(['corr(DSC, DFC_{XCov}) over time, time window ',num2str(new_bin_minutes),'min'])
        title('XCov')
        xlim([minX,maxX])
        text(h2.BinEdges(2),max(h2.Values)-0.04,['mean=',num2str(meanCorrXCov.dc1.varyTJ(bIdx,TJ_idx),3)],'FontSize',18);
        bin=[];
        for idx=1:numel(exampleSyn.TE)
            tmp=find(exampleSyn.XCov(idx)>h2.BinEdges);
            bin(idx)=tmp(end);
            xval=0.5*(h2.BinEdges(bin(idx))+h2.BinEdges(bin(idx)+1));
            yval = (h2.Values(bin(idx)))-(sum(bin==tmp(end))-1)*0.01;
            scatter(xval,yval,[],cols(idx,:),'filled');
        end
    end
    if doXCorr
        subplot(1,nPlots,4)
        h3=histogram(pairCorrXCorr.dc1.varyTJ{bIdx},edges,'normalization','probability');
        xlabel('correlation')
        %title(['corr(DSC, DFC_{XCorr}) over time, time window ',num2str(new_bin_minutes),'min'])
        title('XCorr')
        xlim([minX,maxX])
        text(h2.BinEdges(2),max(h3.Values)-0.03,['mean=',num2str(meanCorrXCorr.dc1.varyTJ(bIdx))],'FontSize',18);
    end
end
%% Plot
cols = [0.4660 0.6740 0.1880;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0 0.4470 0.7410;0.8500 0.3250 0.0980];
%corrType = {'min1','varyTJ','NsampleMatch'};
corrType = {'NsampleMatch'};

%cols = [1 0 0; 0 1 0; 0 0 1];

for dcIdx=1:2
    dcLab=dc_label{dcIdx};
    for idx=1:numel(corrType)
        figure()
        % on the diagonal we have T = TJ
        l1=plot(new_bin_minutes_range(1:size(meanCorrTE.(dcLab).(corrType{idx}),1)),diag(meanCorrTE.(dcLab).(corrType{idx})),'color',cols(1,:),'linewidth',2);
        hold on
        l2=plot(new_bin_minutes_range(1:size(meanCorrXCov.(dcLab).(corrType{idx}),1)),diag(meanCorrXCov.(dcLab).(corrType{idx})),'color',cols(3,:),'linewidth',2);
        l3=plot(new_bin_minutes_range(1:size(meanCorrHOTE.(dcLab).(corrType{idx}),1)),diag(meanCorrHOTE.(dcLab).(corrType{idx})),'color',cols(2,:),'linewidth',2);
        legend('TE','XCov','HOTE','Autoupdate','off')
        shadedErrorBar(new_bin_minutes_range(1:size(meanCorrTE.(dcLab).(corrType{idx}),1)),diag(meanCorrTE.(dcLab).(corrType{idx})),diag(semCorrTE.(dcLab).(corrType{idx})),'lineProps',{'color',cols(1,:)})
        shadedErrorBar(new_bin_minutes_range(1:size(meanCorrXCov.(dcLab).(corrType{idx}),1)),diag(meanCorrXCov.(dcLab).(corrType{idx})),diag(semCorrXCov.(dcLab).(corrType{idx})),'lineProps',{'color',cols(3,:)})
        shadedErrorBar(new_bin_minutes_range(1:size(meanCorrHOTE.(dcLab).(corrType{idx}),1)),diag(meanCorrHOTE.(dcLab).(corrType{idx})),diag(semCorrHOTE.(dcLab).(corrType{idx})),'lineProps',{'color',cols(2,:)})
        xlabel('window size [min]')
        ylabel('correlation')
        xlim([1 max(new_bin_minutes_range(1:size(meanCorrTE.(dcLab).(corrType{idx}),1)))])
        ylim([0.05 0.85])
        title(['Mean correlation between DFC and DSC, ',corrType{idx},'; delay cons.=',dcLab(3)])
    end

end
%

for idx=1:numel(corrType)
    figure()
    % on the diagonal we have T = TJ
    meanCorrTE.diff.(corrType{idx})=diag(meanCorrTE.dc1.(corrType{idx}))-diag(meanCorrTE.dc0.(corrType{idx}));
    l1=plot(new_bin_minutes_range(1:size(meanCorrTE.(dcLab).(corrType{idx}),1)),meanCorrTE.diff.(corrType{idx}),'color',cols(1,:),'linewidth',2);
    hold on
    meanCorrXCov.diff.(corrType{idx})=diag(meanCorrXCov.dc1.(corrType{idx}))-diag(meanCorrXCov.dc0.(corrType{idx}));
    l2=plot(new_bin_minutes_range(1:size(meanCorrXCov.(dcLab).(corrType{idx}),1)),meanCorrXCov.diff.(corrType{idx}),'color',cols(3,:),'linewidth',2);
    meanCorrHOTE.diff.(corrType{idx})=diag(meanCorrHOTE.dc1.(corrType{idx}))-diag(meanCorrHOTE.dc0.(corrType{idx}));
    l3=plot(new_bin_minutes_range(1:size(meanCorrHOTE.(dcLab).(corrType{idx}),1)),meanCorrHOTE.diff.(corrType{idx}),'color',cols(2,:),'linewidth',2);
    legend('TE','XCov','HOTE','Autoupdate','off')
    shadedErrorBar(new_bin_minutes_range(1:size(meanCorrTE.diff.(corrType{idx}),1)),meanCorrTE.diff.(corrType{idx}),semCorrTE.diff.(corrType{idx}),'lineProps',{'color',cols(1,:)})
    shadedErrorBar(new_bin_minutes_range(1:size(meanCorrXCov.diff.(corrType{idx}),1)),meanCorrXCov.diff.(corrType{idx}),semCorrXCov.diff.(corrType{idx}),'lineProps',{'color',cols(3,:)})
    shadedErrorBar(new_bin_minutes_range(1:size(meanCorrHOTE.diff.(corrType{idx}),1)),meanCorrHOTE.diff.(corrType{idx}),semCorrHOTE.diff.(corrType{idx}),'lineProps',{'color',cols(2,:)})
    xlabel('window size [min]')
    ylabel('correlation')
    xlim([1 max(new_bin_minutes_range(1:size(meanCorrTE.(dcLab).(corrType{idx}),1)))])
    title(['Mean correlation gain from delay consistency, ',corrType{idx}])
end
%% Plot
figure()
plot(new_bin_minutes_range,diag(meanCorrTE.dc1.varyTJ)) % on the diagonal we have T = TJ
hold on
plot(new_bin_minutes_range,diag(meanCorrXCov.dc1.varyTJ))
plot(new_bin_minutes_range,diag(meanCorrHOTE.dc1.varyTJ))
legend('TE','XCov','HOTE')
xlabel('recording length')
ylabel('correlation')
title('Mean correlation between DFC and DSC')

figure()
plot(new_bin_minutes_range,diag(meanCorrTEdiff.dc1.varyTJ)) % on the diagonal we have T = TJ
hold on
plot(new_bin_minutes_range,diag(meanCorrXCovdiff.dc1.varyTJ))
plot(new_bin_minutes_range,diag(meanCorrHOTEdiff.dc1.varyTJ))
legend('TE','XCov','HOTE')
xlabel('recording length')
ylabel('\DeltaW correlation')
title('Mean correlation between DFC and DSC')

figure()
plot(new_bin_minutes_range,diag(meanCorrTEdiff.dc1.NsampleMatch)) % on the diagonal we have T = TJ
hold on
plot(new_bin_minutes_range,diag(meanCorrXCovdiff.dc1.NsampleMatch))
plot(new_bin_minutes_range,diag(meanCorrHOTEdiff.dc1.NsampleMatch))
legend('TE','XCov','HOTE')
xlabel('recording length')
ylabel('\DeltaW correlation')
title('Mean correlation between DFC and DSC')
%% Empirical STDP measure
% Load spikeTrains, s_time and post
load('E:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\full_seed1_IzhiModel_v15_90min_3sTJ_varyWindow_1_DelConsist_06-Sep-2022','spikeTrains','s_time','s_time_downsampled','post','gtDelay')
%%
% Params setting
prec_min = 1; % selected precision in minutes
bIdx = find(new_bin_minutes_range==prec_min);
sel_TJ = time_jumps(3);
sel_TJ_sec = sel_TJ*60; % s_time_downsampled is in seconds
delay_range = [-20:20];

%downsamp_prec = 5; %s

% Get spike times from trains
spikeTimes = SparseToASDF(spikeTrains, 1);

% Get GT connectivity donsampled to match TJ
s_time_perm = permute(s_time,[1,3,2]);
s_time_TJdawnsamp = temporal_rebinning(s_time_perm,60,'movmean',sel_TJ_sec);
s_time_TJdawnsamp = permute(s_time_TJdawnsamp,[1,3,2]);

t_idxs = 1:sel_TJ/(time_jumps(1)):size(match_XCov.dc1{bIdx},2);
XCov_TJ_downsampled = match_XCov.dc1{bIdx}(:,t_idxs);

% Initialize data structures
empirical_STDP_TJ = cell(1,numel(delay_range));
empirical_STDP_gt = cell(1,numel(delay_range));
empirical_STDP_gt_downsamp = cell(1,numel(delay_range));

for pairIdx = 1:size(XCov_connPairs,1)/10
    pairIdx
    emit = XCov_connPairs(pairIdx,1);
    rec = XCov_connPairs(pairIdx,2);
    gtSynDel = gtDelay(emit,rec);
    % Check if you are also considering synapses from inh neurons (not
    % plastic)
    spikeTimes_emit = spikeTimes{emit};
    spikeTimes_rec = spikeTimes{rec};
    for t = spikeTimes_emit
        stdp_rec = find((spikeTimes_rec-t >= delay_range(1)+gtSynDel)&(spikeTimes_rec-t <= delay_range(end)+gtSynDel));
        if ~isempty(stdp_rec)
            time_diff = spikeTimes_rec(stdp_rec)-t-gtSynDel; % w.r.t. gtDelay
            for td = 1:numel(time_diff)
                dIdx = find(delay_range==time_diff(td));
                rec_idx = find(post(emit,:)==rec); % relative idx of the receiver in the syn. weights evolution matrix
                recTime_TJ = floor(spikeTimes_rec(stdp_rec(td))/(sel_TJ*60*1000)); % get Time in units of the selected TJ
                recTime_gt = floor(spikeTimes_rec(stdp_rec(td))/(1000)); % get Time in units of the selected TJ
                if ((recTime_TJ>0)&(recTime_TJ<size(XCov_TJ_downsampled,2)))
                    dW_nextTime_TJ = XCov_TJ_downsampled(pairIdx,recTime_TJ+1)-XCov_TJ_downsampled(pairIdx,recTime_TJ); % get dW in the next timeJump
                    empirical_STDP_TJ{dIdx}=[empirical_STDP_TJ{dIdx},dW_nextTime_TJ];
                end
                if ((recTime_gt>0)&(recTime_gt<size(s_time,3)-60))
                    dW_nextTime_sec = s_time(emit,rec_idx,recTime_gt+60)-s_time(emit,rec_idx,recTime_gt); % get dW in the next second (GT)
                    empirical_STDP_gt{dIdx}=[empirical_STDP_gt{dIdx},dW_nextTime_sec];
                end
                if ((recTime_TJ>0)&(recTime_TJ<size(XCov_TJ_downsampled,2)))
                    dW_nextTime_gt_downsamp = s_time_TJdawnsamp(emit,rec_idx,recTime_TJ+1)-s_time_TJdawnsamp(emit,rec_idx,recTime_TJ); % get dW in the next timeJump
                    empirical_STDP_gt_downsamp{dIdx}=[empirical_STDP_gt_downsamp{dIdx},dW_nextTime_gt_downsamp];
                end
            end
        end
    end
end

for delIdx = 1:numel(delay_range)
    mean_empirical_STDP_TJ(delIdx) = mean(empirical_STDP_TJ{delIdx});
    sem_empirical_STDP_TJ(delIdx) = std(empirical_STDP_TJ{delIdx})/sqrt(numel(empirical_STDP_TJ{delIdx}));
    mean_empirical_STDP_gt(delIdx) = mean(empirical_STDP_gt{delIdx});
    sem_empirical_STDP_gt(delIdx) = std(empirical_STDP_gt{delIdx})/sqrt(numel(empirical_STDP_gt{delIdx}));
    mean_empirical_STDP_TJ_downsamp(delIdx) = mean(empirical_STDP_gt_downsamp{delIdx});
    sem_empirical_STDP_TJ_downsamp(delIdx) = std(empirical_STDP_gt_downsamp{delIdx})/sqrt(numel(empirical_STDP_gt_downsamp{delIdx}));
end

% save(['empirical_stdp_',date])
%%
figure()
subplot(1,2,2)
hold on
plot(delay_range,mean_empirical_STDP_TJ)
errorbar(delay_range,mean_empirical_STDP_TJ,sem_empirical_STDP_TJ)
title('Measrued XCov')
xlabel('Delay [ms]')
ylabel('\DeltaW [corr.]')
subplot(1,2,1)
hold on
plot(delay_range,mean_empirical_STDP_TJ)
errorbar(delay_range,mean_empirical_STDP_TJ_downsamp,sem_empirical_STDP_TJ_downsamp)
title('Downsampled GT')
xlabel('Delay [ms]')
ylabel('\DeltaW [a.u.]')
sgtitle('STDP rule inference','FontSize',14)
% subplot(1,3,3)
% hold on
% plot(delay_range,mean_empirical_STDP_gt)
% errorbar(delay_range,mean_empirical_STDP_gt,sem_empirical_STDP_gt)
% title('GT 1min lag')

%% Fit double exponential to empirical STDP
x_start=find(delay_range == 1);
x_plus = x_start:numel(delay_range);
f = fit(delay_range(x_plus),mean_empirical_STDP_TJ(x_plus),'exp1');