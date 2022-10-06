%% Script to reproduce plots in Figure 4C-E

doHO = 1; doXCov = 1;

clear all; close all;

% Load simulation results
PATH = ['path_to_main_directory\results\static_connectivity\']; % path to directory where the static connectivity results have been saved

params.doOverlap = 1;
params.nRep = 5;
minutes_range = [5 10 20 30 50 70 90 120 180];
for rIdx = 1:params.nRep
    for mIdx = 1:numel(minutes_range)
        filename = dir([PATH,'staticConnInference_rep',num2str(rIdx),'_min',num2str(minutes_range(mIdx)),'_*']);
        
        tmpFile = load([PATH,filename.name],...
            'TErec','TEprec','HOTErec','HOTEprec','XCovrec','XCovprec','XCorrrec','XCorrprec','TEdelays','XCovDelays','HOTEdelays',...
            'gtConn','peakTE','peakHOTE','peakXCov','peakXCorr','minutes_range','delayCorr','delayErr','gtDelay');
        
        tmpMinIdx = find(tmpFile.minutes_range==minutes_range(mIdx));

        gtConn{rIdx,mIdx} = tmpFile.gtConn; gtDelay{rIdx,mIdx} = tmpFile.gtDelay; peakTE{rIdx,mIdx} = tmpFile.peakTE; peakXCorr{rIdx,mIdx} = tmpFile.peakXCorr; 
        peakHOTE{rIdx,mIdx} = tmpFile.peakHOTE; peakXCov{rIdx,mIdx} = tmpFile.peakXCov; TEdelays{rIdx,mIdx} = tmpFile.TEdelays; HOTEdelays{rIdx,mIdx} = tmpFile.HOTEdelays; XCovDelays{rIdx,mIdx} = tmpFile.XCovDelays;
        measureNames = fieldnames(tmpFile.delayCorr); 
        for measIdx = 1:numel(measureNames)
            measLab = measureNames{measIdx};
            delayCorr.(measLab)(rIdx,mIdx) = tmpFile.delayCorr.(measLab);
            delayErr.(measLab)(rIdx,mIdx) = tmpFile.delayErr.(measLab);
        end
    end
end

legend_labels=cell(1,numel(minutes_range)+1);
for i = 1:numel(minutes_range)
    legend_labels{i}=[num2str(minutes_range(i)),'min'];
end
legend_labels{i+1}='null hyp.';

%% Delay correlation and delay avg error (panels 4C and 4E)

cols = [0.4660 0.6740 0.1880;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0 0.4470 0.7410;0.8500 0.3250 0.0980];

% Plot panel 4D
figure()
hold on
plot(minutes_range,mean(delayCorr.TE,1),'linewidth',1.5,'color',cols(1,:))
h1=errorbar(minutes_range,mean(delayCorr.TE,1),std(delayCorr.TE,[],1)/sqrt(size(delayCorr.TE,1)),'linewidth',1.5,'color',cols(1,:));
plot(minutes_range,mean(delayCorr.HOTE,1),'linewidth',1.5,'color',cols(2,:))
h2=errorbar(minutes_range,mean(delayCorr.HOTE,1),std(delayCorr.HOTE,[],1)/sqrt(size(delayCorr.HOTE,1)),'linewidth',1.5,'color',cols(2,:));
plot(minutes_range,mean(delayCorr.XCov,1),'linewidth',1.5,'color',cols(3,:))
h3=errorbar(minutes_range,mean(delayCorr.XCov,1),std(delayCorr.XCov,[],1)/sqrt(size(delayCorr.XCov,1)),'linewidth',1.5,'color',cols(3,:));
xlim([minutes_range(1),minutes_range(end)])
legend([h1,h2,h3],'TE','HOTE','XCov')
xlabel('recording length [minutes]')
ylabel('correlation')
title(['Correlation with GT delays'],'FontSize',16)

% Plot panel 4E
figure()
hold on
plot(minutes_range,mean(delayErr.TE,1),'linewidth',1.5,'color',cols(1,:))
h1=errorbar(minutes_range,mean(delayErr.TE,1),std(delayErr.TE,[],1)/sqrt(size(delayErr.TE,1)),'linewidth',1.5,'color',cols(1,:));
plot(minutes_range,mean(delayErr.HOTE,1),'linewidth',1.5,'color',cols(2,:))
h2=errorbar(minutes_range,mean(delayErr.HOTE,1),std(delayErr.HOTE,[],1)/sqrt(size(delayErr.HOTE,1)),'linewidth',1.5,'color',cols(2,:));
plot(minutes_range,mean(delayErr.XCov,1),'linewidth',1.5,'color',cols(3,:))
h3=errorbar(minutes_range,mean(delayErr.XCov,1),std(delayErr.XCov,[],1)/sqrt(size(delayErr.XCov,1)),'linewidth',1.5,'color',cols(3,:));
xlim([minutes_range(1),minutes_range(end)])
legend([h1,h2,h3],'TE','HOTE','XCov')
xlabel('recording length [minutes]')
ylabel('average error [ms]')
title(['Average delays error'],'FontSize',16)

%% Delays scatter plots (Fig.4C)

doHO = 1; doXCov = 1;

fsize = 16;
markers_scale = 1;
% TE
figure()
ax = gca;
hold on
for gtd = 1:20
    for estd = 1:50
        weight = 0;
        for repIdx = 1:params.nRep
            weight = weight + sum((gtDelay{repIdx,end}(~isnan(gtDelay{repIdx,end}(:)))==gtd) & (TEdelays{repIdx,end}(~isnan(gtDelay{repIdx,end}(:)))==estd));
        end
        if weight > 0
            scatter(gtd, estd,markers_scale*weight, cols(1,:),'filled','MarkerEdgeColor','k');
        end
    end
end
text(5, 32, ['delays correlation = ',num2str(delayCorr.TE(end),2)],'FontSize',fsize);
text(5, 28, ['avg delay err = ',num2str(delayErr.TE(end),2),'ms'],'FontSize',fsize);
xlim([0,21])
plot(0:50,0:50,'k--','linewidth',2)
xlabel('Ground truth delay [ms]')
ylabel('TE delay [ms]')
title(['True vs TE delays, ',num2str(minutes_range(end)),' minutes'])
ax.FontSize = fsize;

% HOTE
if doHO
    figure()
    ax = gca;
    hold on
    for gtd = 1:20
        for estd = 1:50
            weight = 0;
            for repIdx = 1:params.nRep
                weight = weight + sum((gtDelay{repIdx,end}(~isnan(gtDelay{repIdx,end}(:)))==gtd) & (HOTEdelays{repIdx,end}(~isnan(gtDelay{repIdx,end}(:)))==estd));
            end
            if weight > 0
                scatter(gtd, estd,markers_scale*weight, cols(2,:),'filled','MarkerEdgeColor','k');
            end
        end
    end
    text(5, 32, ['delays correlation = ',num2str(delayCorr.HOTE(end),2)],'FontSize',fsize);
    text(5, 28, ['avg delay err = ',num2str(delayErr.HOTE(end),2),'ms'],'FontSize',fsize);
    %text(25, 0.18, ['delays match = ',num2str(delayMatch.TE(mIdx)*100),'%'],'FontSize',14);
    xlim([0,21])
    ylim([0,51])
    plot(0:50,0:50,'k--','linewidth',2)
    xlabel('Ground truth delay [ms]')
    ylabel('HOTE delay [ms]')
    title(['True vs HOTE delays, ',num2str(minutes_range(end)),' minutes'])
    ax.FontSize = fsize;
end

if doXCov
    % XCov
    figure();
    ax = gca;
    hold on
    for gtd = 1:20
        for estd = 1:50
            weight = 0;
            for repIdx = 1:params.nRep
                weight = weight + sum((gtDelay{repIdx,end}(~isnan(gtDelay{repIdx,end}(:)))==gtd) & (XCovDelays{repIdx,end}(~isnan(gtDelay{repIdx,end}(:)))==estd));
            end
            if weight > 0
                scatter(gtd, estd,markers_scale*weight, cols(3,:),'filled','MarkerEdgeColor','k');
            end
        end
    end
    text(5, 32, ['delays correlation = ',num2str(mean(delayCorr.XCov(:,end)),2)],'FontSize',fsize);
    text(5, 28, ['avg delay err = ',num2str(mean(delayErr.XCov(:,end)),2),'ms'],'FontSize',fsize);
    %text(25, 0.18, ['delays match = ',num2str(delayMatch.TE(mIdx)*100),'%'],'FontSize',14);
    xlim([0,21])
    ylim([0,51])
    plot(0:50,0:50,'k--','linewidth',2)
    xlabel('Ground truth delay [ms]')
    ylabel('XCov delay [ms]')
    title(['True vs XCov delays, ',num2str(minutes_range(end)),' minutes'])
    ax.FontSize = fsize;
end
