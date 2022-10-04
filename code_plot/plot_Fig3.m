%% Static connectivity plot
% Working with single min files + overlap measure
clear all; close all;

% Load simulation results
PATH = 'E:\Lemke_mat_files\sleepAnalysis\static_connectivity\BI_resubmission\';
params.doOverlap = 1;
params.nRep = 5;
minutes_range = [5 10 20 30 50 70 90 120 180];
for rIdx = 1:params.nRep
    rIdx
    for mIdx = 1:numel(minutes_range)
        filename = dir([PATH,'staticConnInference_rep',num2str(rIdx),'_min',num2str(minutes_range(mIdx)),'_*']);

        tmpFile = load([PATH,filename.name],...
            'TErec','TEprec','HOTErec','HOTEprec','XCovrec','XCovprec','XCorrrec','XCorrprec',...
            'gtConn','peakTE','peakHOTE','peakXCov','peakXCorr','minutes_range');
        
        tmpMinIdx = find(tmpFile.minutes_range==minutes_range(mIdx));

        gtConn{rIdx,mIdx} = tmpFile.gtConn; peakTE{rIdx,mIdx} = tmpFile.peakTE; peakXCorr{rIdx,mIdx} = tmpFile.peakXCorr; 
        peakHOTE{rIdx,mIdx} = tmpFile.peakHOTE; peakXCov{rIdx,mIdx} = tmpFile.peakXCov; 
    end
end

legend_labels=cell(1,numel(minutes_range)+1);
for i = 1:numel(minutes_range)
    legend_labels{i}=[num2str(minutes_range(i)),'min'];
end
legend_labels{i+1}='null hyp.';

%% Compute measures overlap (XCov, HOTE)

PR_measures = compute_PR_v1(gtConn,peakTE,peakHOTE,peakXCov,peakXCorr,minutes_range,params);
save(['Overlap_180min_',date])

%% Distribution of synaptic weights
PATH = 'E:\Lemke_mat_files\sleepAnalysis\dynamic_connectivity\BI_resubmission\';
filename = ['dynamicConnInference_rep1_180min_varyWindow_1_DelConsist_15-Sep-2022'];
load([PATH,filename],'gtDelay','s_time','minutes_range','Ne','M','staticConnMeasures');

prctiles = [1, 5, 10];
weights = abs(staticConnMeasures.peakXCov(:));
connMap = abs(staticConnMeasures.peakXCov);
cols = distinguishable_colors(numel(prctiles));
cm = [1,1,1;cols];

colormap(cm);

tmpCorrPast = zeros(100,100);
totPercConn = zeros(100,100);
figure()
hold on
histogram(weights,'normalization','probability')
for prcIdx = 1:numel(prctiles)
    tmpPerc = prctiles(prcIdx);
    tmpThresh(prcIdx) = prctile(weights,100-tmpPerc);
    h{prcIdx}=xline(tmpThresh(prcIdx),'color',cols(prcIdx,:));
    tmpConn = (connMap>tmpThresh(prcIdx));
    newConn = tmpConn-tmpCorrPast;
    tmpCorrPast = tmpConn;
    newConn = logical(newConn);
    totPercConn(newConn) = prcIdx;
end
title('Measured synaptic weights distribution','FontSize',14)
xlabel('XCov','FontSize',14)
ylabel('Probability','FontSize',14)
legend([h{1},h{2},h{3}],'99th','95th','90th')

figure()
imagesc(totPercConn)
colormap(cm)
title('Top inferred synapses','FontSize',14)
xline(80.5,'k--')
xlabel('Receiver','FontSize',14)
yline(80.5,'k--')
ylabel('Emitter','FontSize',14)

%% Plot precision recall curve for different measures
% Average PR curves across repetitions
cols = [0.4660 0.6740 0.1880;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0 0.4470 0.7410;0.8500 0.3250 0.0980];

mIdx = 9;

figure()
hold on
plot(squeeze(mean(PR_measures.TE.rec(:,mIdx,:),1)), squeeze(mean(PR_measures.TE.prec(:,mIdx,:),1)),'x-','linewidth',2,'color',cols(1,:))

plot(squeeze(mean(PR_measures.HOTE.rec(:,mIdx,:),1)), squeeze(mean(PR_measures.HOTE.prec(:,mIdx,:),1)),'x-','linewidth',2,'color',cols(2,:))

plot(squeeze(mean(PR_measures.XCov.rec(:,mIdx,:),1)), squeeze(mean(PR_measures.XCov.prec(:,mIdx,:),1)),'x-','linewidth',2,'color',cols(3,:))

plot(squeeze(mean(PR_measures.XCorr.rec(:,mIdx,:),1)), squeeze(mean(PR_measures.XCorr.prec(:,mIdx,:),1)),'x-','linewidth',2,'color',cols(4,:))

if params.doOverlap
    plot(squeeze(mean(PR_measures.HOTE_XCov.rec(:,mIdx,:),1)), squeeze(mean(PR_measures.HOTE_XCov.prec(:,mIdx,:),1)),'x-','linewidth',2,'color',cols(5,:))
end

yline(0.1,'k--','LineWidth',1.5)
legend('TE','HOTE','XCov','XCorr','Overlap_{HOTE-XCov}','null hyp.')
ylim([0,1])
xlabel('recall')
ylabel('precision')
title(['All metrics PR, ',num2str(minutes_range(mIdx)),' minutes'],'FontSize',16)

%% Plot AUPR with time for different measures
cols = [0.4660 0.6740 0.1880;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0 0.4470 0.7410;0.8500 0.3250 0.0980];
for rIdx = 1:params.nRep
    for mIdx = 1:numel(minutes_range)
        AUPR.TE(rIdx,mIdx) = trapz(squeeze(PR_measures.TE.rec(rIdx,mIdx,end:-1:1)), squeeze(PR_measures.TE.prec(rIdx,mIdx,end:-1:1)));
        AUPR.HOTE(rIdx,mIdx) = trapz(squeeze(PR_measures.HOTE.rec(rIdx,mIdx,end:-1:1)), squeeze(PR_measures.HOTE.prec(rIdx,mIdx,end:-1:1))); 
        AUPR.XCorr(rIdx,mIdx) = trapz(squeeze(PR_measures.XCorr.rec(rIdx,mIdx,end:-1:1)), squeeze(PR_measures.XCorr.prec(rIdx,mIdx,end:-1:1))); 
        AUPR.XCov(rIdx,mIdx) = trapz(squeeze(PR_measures.XCov.rec(rIdx,mIdx,end:-1:1)), squeeze(PR_measures.XCov.prec(rIdx,mIdx,end:-1:1))); 
        if params.doOverlap
            AUPR.HOTE_XCov(rIdx,mIdx) = trapz(squeeze(PR_measures.HOTE_XCov.rec(rIdx,mIdx,end:-1:1)), squeeze(PR_measures.HOTE_XCov.prec(rIdx,mIdx,end:-1:1))); 
        end
    end
end

figure()
h1=plot(minutes_range, mean(AUPR.TE,1),'linewidth',2,'color',cols(1,:));
hold on
errorbar(minutes_range, mean(AUPR.TE,1), std(AUPR.TE,[],1),'color',h1.Color)
h2=plot(minutes_range, mean(AUPR.HOTE,1),'linewidth',2,'color',cols(2,:));
errorbar(minutes_range, mean(AUPR.HOTE,1), std(AUPR.HOTE,[],1),'color',h2.Color)
h3=plot(minutes_range, mean(AUPR.XCov,1),'linewidth',2,'color',cols(3,:));
errorbar(minutes_range, mean(AUPR.XCov,1), std(AUPR.XCov,[],1),'color',h3.Color)
h4=plot(minutes_range, mean(AUPR.XCorr,1),'linewidth',2,'color',cols(4,:));
errorbar(minutes_range, mean(AUPR.XCorr,1), std(AUPR.XCorr,[],1),'color',h4.Color)
if params.doOverlap
    h5=plot(minutes_range, mean(AUPR.HOTE_XCov,1),'linewidth',2,'color',cols(5,:));
    errorbar(minutes_range, mean(AUPR.HOTE_XCov,1), std(AUPR.HOTE_XCov,[],1),'color',h5.Color)
end
legend([h1,h2,h3,h4,h5],'TE','HOTE','XCov','XCorr','Overlap_{HOTE-XCov}')
ylim([0,1])
xlim([minutes_range(1),minutes_range(end)])
xlabel('simulation length')
ylabel('AUPR')
title(['AUPR, with simulation length'],'FontSize',16)

%% Plot 5th and 10th prctile precision
figure()
subplot(2,1,1)
plot(minutes_range, mean(PR_measures.TE.prec(:,:,90),1),'linewidth',2,'color',cols(1,:))
% set(gca,'YScale','log');
hold on
errorbar(minutes_range, mean(PR_measures.TE.prec(:,:,90),1),std(PR_measures.TE.prec(:,:,90),[],1),'linewidth',2,'color',cols(1,:))
plot(minutes_range, mean(PR_measures.HOTE.prec(:,:,90),1),'linewidth',2,'color',cols(2,:))
errorbar(minutes_range, mean(PR_measures.HOTE.prec(:,:,90),1),std(PR_measures.HOTE.prec(:,:,90),[],1),'linewidth',2,'color',cols(2,:))
plot(minutes_range, mean(PR_measures.XCov.prec(:,:,90),1),'linewidth',2,'color',cols(3,:))
errorbar(minutes_range, mean(PR_measures.XCov.prec(:,:,90),1),std(PR_measures.XCov.prec(:,:,90),[],1),'linewidth',2,'color',cols(3,:))
plot(minutes_range, mean(PR_measures.XCorr.prec(:,:,90),1),'linewidth',2,'color',cols(4,:))
errorbar(minutes_range, mean(PR_measures.XCorr.prec(:,:,90),1),std(PR_measures.XCorr.prec(:,:,90),[],1),'linewidth',2,'color',cols(4,:))
if params.doOverlap
    plot(minutes_range, mean(PR_measures.HOTE_XCov.prec(:,:,90),1),'linewidth',2,'color',cols(5,:))
    errorbar(minutes_range, mean(PR_measures.HOTE_XCov.prec(:,:,90),1),std(PR_measures.HOTE_XCov.prec(:,:,90),[],1),'linewidth',2,'color',cols(5,:))
end
ylim([0,1])
xlim([minutes_range(1),minutes_range(end)])
xlabel('simulation length [min]')
ylabel('precision')
title('90th perctile (1000 pairs)')

subplot(2,1,2)
h1=plot(minutes_range, mean(PR_measures.TE.prec(:,:,95),1),'linewidth',2,'color',cols(1,:));
hold on
errorbar(minutes_range, mean(PR_measures.TE.prec(:,:,95),1),std(PR_measures.TE.prec(:,:,95),[],1),'linewidth',2,'color',cols(1,:))
h2=plot(minutes_range, mean(PR_measures.HOTE.prec(:,:,95),1),'linewidth',2,'color',cols(2,:));
errorbar(minutes_range, mean(PR_measures.HOTE.prec(:,:,95),1),std(PR_measures.HOTE.prec(:,:,95),[],1),'linewidth',2,'color',cols(2,:))
h3=plot(minutes_range, mean(PR_measures.XCov.prec(:,:,95),1),'linewidth',2,'color',cols(3,:));
errorbar(minutes_range, mean(PR_measures.XCov.prec(:,:,95),1),std(PR_measures.XCov.prec(:,:,95),[],1),'linewidth',2,'color',cols(3,:))
h4=plot(minutes_range, mean(PR_measures.XCorr.prec(:,:,95),1),'linewidth',2,'color',cols(4,:));
errorbar(minutes_range, mean(PR_measures.XCorr.prec(:,:,95),1),std(PR_measures.XCorr.prec(:,:,95),[],1),'linewidth',2,'color',cols(4,:))
if params.doOverlap
    h5=plot(minutes_range, mean(PR_measures.HOTE_XCov.prec(:,:,95),1),'linewidth',2,'color',cols(5,:));
    errorbar(minutes_range, mean(PR_measures.HOTE_XCov.prec(:,:,95),1),std(PR_measures.HOTE_XCov.prec(:,:,95),[],1),'linewidth',2,'color',cols(5,:))
end
legend([h1,h2,h3,h4,h5],'TE','HOTE','XCov','XCorr','Overlap_{HOTE-XCov}')
ylim([0,1])
xlim([minutes_range(1),minutes_range(end)])
xlabel('simulation length [min]')
ylabel('precision')
title('95th perctile (500 pairs)')
sgtitle(['Precision for fixed percentile'])


% Plot Overlap - maximum
diffOverlap = PR_measures.HOTE_XCov.prec(:,:,90) - max([PR_measures.XCov.prec(:,:,90);PR_measures.HOTE.prec(:,:,90)]);
meanDiffOverlap = mean(diffOverlap,1);
stdDiffOverlap = std(diffOverlap,[],1);

figure()
plot(minutes_range, meanDiffOverlap,'linewidth',2,'color',cols(5,:))
errorbar(minutes_range, meanDiffOverlap,stdDiffOverlap,'color',cols(5,:))
%% Plot fraction of inferred connections in all subgroups

for rIdx = 1:params.nRep
    for mIdx = 1:numel(minutes_range)
        % CM = connectivity measures
        CM{rIdx,mIdx}.gtConn.all=gtConn{rIdx,mIdx};
        CM{rIdx,mIdx}.peakXCov.all=peakXCov{rIdx,mIdx};
        CM{rIdx,mIdx}.peakXCorr.all=peakXCorr{rIdx,mIdx};
        CM{rIdx,mIdx}.peakHOTE.all = peakHOTE{rIdx,mIdx};
        CM{rIdx,mIdx}.peakTE.all = peakTE{rIdx,mIdx};

        measures = {'gtConn','peakXCov','XCovmap','peakHOTE','peakTE','HOTEmap','TEmap','XCorrmap','XCovmap','overlap_HOTE_XCov'};
        exc_idx = 1:80;
        inh_idx = 81:100;
        sel_thresh = 90;

        HOTEthr = prctile(CM{rIdx,mIdx}.peakHOTE.all(:),sel_thresh);
        CM{rIdx,mIdx}.HOTEmap.all = (CM{rIdx,mIdx}.peakHOTE.all >= HOTEthr);
        TEthr = prctile(CM{rIdx,mIdx}.peakTE.all(:),sel_thresh);
        CM{rIdx,mIdx}.TEmap.all = (CM{rIdx,mIdx}.peakTE.all >= TEthr);
        XCovthr = prctile(abs(CM{rIdx,mIdx}.peakXCov.all(:)),sel_thresh);
        CM{rIdx,mIdx}.XCovmap.all = (abs(CM{rIdx,mIdx}.peakXCov.all) >= XCovthr);
        XCorrthr = prctile(abs(CM{rIdx,mIdx}.peakXCorr.all(:)),sel_thresh);
        CM{rIdx,mIdx}.XCorrmap.all = (abs(CM{rIdx,mIdx}.peakXCorr.all) >= XCorrthr);
        CM{rIdx,mIdx}.overlap_HOTE_XCov.all = measures_overlap_prctile(peakXCov{rIdx,mIdx},peakHOTE{rIdx,mIdx},sel_thresh);

        for meIdx = 1:numel(measures)
            mLab = measures{meIdx};
            CM{rIdx,mIdx}.(mLab).EE = CM{rIdx,mIdx}.(mLab).all(1:80,1:80);
            CM{rIdx,mIdx}.(mLab).EI = CM{rIdx,mIdx}.(mLab).all(1:80,81:100);
            CM{rIdx,mIdx}.(mLab).IE = CM{rIdx,mIdx}.(mLab).all(81:100,1:80);
            CM{rIdx,mIdx}.(mLab).II = CM{rIdx,mIdx}.(mLab).all(81:100,81:100);
            CM{rIdx,mIdx}.(mLab).allE = CM{rIdx,mIdx}.(mLab).all(1:80,:);
            CM{rIdx,mIdx}.(mLab).allI = CM{rIdx,mIdx}.(mLab).all(81:100,:);
        end
    end
end

%% Plot fraction of connections per group comparing measures
sel_idx = numel(minutes_range);
connGroups = {'EE','EI','IE','II'};
x_spread=[-0.25,-0.15,-0.05,0.05,0.15,0.25];
measures = {'gtConn','HOTEmap','TEmap','XCovmap','XCorrmap','overlap_HOTE_XCov'};
barwidth = 0.1;
figure()
hold on
% cols = distinguishable_colors(6);
cols = [0.4660 0.6740 0.1880;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0 0.4470 0.7410;0.8500 0.3250 0.0980];
cols = [0 0 0; cols];

for rIdx = 1:params.nRep
    for meIdx=1:numel(measures)
        inferredEI.(measures{meIdx}).EE(rIdx) = sum(sum(CM{rIdx,sel_idx}.(measures{meIdx}).EE))/sum(sum(CM{sel_idx}.(measures{meIdx}).all));
        inferredEI.(measures{meIdx}).EI(rIdx) = sum(sum(CM{rIdx,sel_idx}.(measures{meIdx}).EI))/sum(sum(CM{sel_idx}.(measures{meIdx}).all));
        inferredEI.(measures{meIdx}).IE(rIdx) = sum(sum(CM{rIdx,sel_idx}.(measures{meIdx}).IE))/sum(sum(CM{sel_idx}.(measures{meIdx}).all));
        inferredEI.(measures{meIdx}).II(rIdx) = sum(sum(CM{rIdx,sel_idx}.(measures{meIdx}).II))/sum(sum(CM{sel_idx}.(measures{meIdx}).all));
        inferredEI.(measures{meIdx}).overall_err(rIdx)=(abs(sum(sum(CM{rIdx,sel_idx}.(measures{meIdx}).EE))-sum(sum(CM{rIdx,sel_idx}.gtConn.EE)))+abs(sum(sum(CM{rIdx,sel_idx}.(measures{meIdx}).EI))-sum(sum(CM{rIdx,sel_idx}.gtConn.EI)))...
            +abs(sum(sum(CM{rIdx,sel_idx}.(measures{meIdx}).IE))-sum(sum(CM{rIdx,sel_idx}.gtConn.IE)))+abs(sum(sum(CM{rIdx,sel_idx}.(measures{meIdx}).II))-sum(sum(CM{rIdx,sel_idx}.gtConn.II))))/(4*sum(sum(CM{rIdx,sel_idx}.gtConn.all)));
    end
end

clear h
for meIdx=1:numel(measures)
    h(meIdx)=bar(1+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).EE),barwidth,'FaceColor',cols(meIdx,:));
    errorbar(1+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).EE),std(inferredEI.(measures{meIdx}).EE),'color','k');
    scatter(1+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).EE),18,cols(meIdx,:),'filled');
    bar(2+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).EI),barwidth,'FaceColor',cols(meIdx,:));
    errorbar(2+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).EI),std(inferredEI.(measures{meIdx}).EI),'color','k');
    scatter(2+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).EI),18,cols(meIdx,:),'filled');
    bar(3+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).IE),barwidth,'FaceColor',cols(meIdx,:));
    errorbar(3+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).IE),std(inferredEI.(measures{meIdx}).IE),'color','k');
    scatter(3+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).IE),18,cols(meIdx,:),'filled');
    bar(4+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).II),barwidth,'FaceColor',cols(meIdx,:));
    errorbar(4+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).II),std(inferredEI.(measures{meIdx}).II),'color','k');
    scatter(4+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).II),18,cols(meIdx,:),'filled');
    
    % Compute avg error by averaging across groups the differences between single measures and ground truth 
    
    bar(5+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).overall_err(rIdx)),barwidth,'FaceColor',cols(meIdx,:));
    errorbar(5+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).overall_err),std(inferredEI.(measures{meIdx}).overall_err),'color','k');
    scatter(5+x_spread(meIdx),mean(inferredEI.(measures{meIdx}).overall_err),18,cols(meIdx,:),'filled');
end
legend([h],'GT','HOTE','TE','XCov','XCorr','Overlap_{HOTE-XCov}')
ax = gca;
ax.FontSize = 14;
set(gca, 'XTick', [1. 2. 3. 4. 5.], 'XTickLabel', {'E-E', 'E-I', 'I-E', 'I-I', 'Avg. err'})
ylim([0 0.7])
title(['Fraction of synapses in each connectivity group'])

%% E/I estimation performance (Fig.4)
cols = [0.4660 0.6740 0.1880;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0 0.4470 0.7410;0.8500 0.3250 0.0980];

figure()
exc_connections = [];
inh_connections = [];

for rIdx = 1:params.nRep
    exc_connections = [exc_connections;[CM{rIdx,sel_idx}.peakXCov.EE(logical(CM{rIdx,sel_idx}.gtConn.EE));CM{rIdx,sel_idx}.peakXCov.EI(logical(CM{rIdx,sel_idx}.gtConn.EI))]];
    inh_connections = [inh_connections;[CM{rIdx,sel_idx}.peakXCov.IE(logical(CM{rIdx,sel_idx}.gtConn.IE));CM{rIdx,sel_idx}.peakXCov.II(logical(CM{rIdx,sel_idx}.gtConn.II))]];
end
hold on
histogram(exc_connections,10,'normalization','probability','FaceColor',cols(1,:))
histogram(inh_connections,3,'normalization','probability','FaceColor',cols(4,:))
lgd=legend('Excitatory','Inhibitory','Autoupdate','off');
lgd.FontSize=16;
xline(0,'k--','linewidth',1.5)
xlabel('XCov [corr.]')
ylabel('prob.')
title(['XCov values distributions, ',num2str(minutes_range(sel_idx)),' minutes'],'FontSize',16)

for rIdx = 1:params.nRep
    for mIdx = 1:numel(minutes_range)
        EIvals.XCov(1,1)=sum(CM{rIdx,mIdx}.peakXCov.EE(logical(CM{rIdx,mIdx}.gtConn.EE))>=0)+...
            sum(CM{rIdx,mIdx}.peakXCov.EI(logical(CM{rIdx,mIdx}.gtConn.EI))>=0);
        EIvals.XCov(1,2)=sum(CM{rIdx,mIdx}.peakXCov.EE(logical(CM{rIdx,mIdx}.gtConn.EE))<0)+...
            sum(CM{rIdx,mIdx}.peakXCov.EI(logical(CM{rIdx,mIdx}.gtConn.EI))<0);
        EIvals.XCov(2,1)=sum(CM{rIdx,mIdx}.peakXCov.IE(logical(CM{rIdx,mIdx}.gtConn.IE))>=0)+...
            sum(CM{rIdx,mIdx}.peakXCov.II(logical(CM{rIdx,mIdx}.gtConn.II))>=0);
        EIvals.XCov(2,2)=sum(CM{rIdx,mIdx}.peakXCov.IE(logical(CM{rIdx,mIdx}.gtConn.IE))<0)+...
            sum(CM{rIdx,mIdx}.peakXCov.II(logical(CM{rIdx,mIdx}.gtConn.II))<0);

        tp.Exc.XCov(rIdx,mIdx)=EIvals.XCov(1,1)/(EIvals.XCov(1,1)+EIvals.XCov(1,2));
        tp.Inh.XCov(rIdx,mIdx)=EIvals.XCov(2,2)/(EIvals.XCov(2,1)+EIvals.XCov(2,2));
        fp.Exc.XCov(rIdx,mIdx)=EIvals.XCov(1,2)/(EIvals.XCov(1,1)+EIvals.XCov(1,2));
        fp.Inh.XCov(rIdx,mIdx)=EIvals.XCov(2,1)/(EIvals.XCov(2,1)+EIvals.XCov(2,2));
    end
end

figure()
h1=plot(minutes_range,mean(tp.Exc.XCov,1),'linewidth',2,'color',cols(1,:));
hold on
errorbar(minutes_range,mean(tp.Exc.XCov,1),std(tp.Exc.XCov,[],1),'linewidth',2,'color',cols(1,:))
h2=plot(minutes_range,mean(tp.Inh.XCov,1),'linewidth',2,'color',cols(4,:));
errorbar(minutes_range,mean(tp.Inh.XCov,1),std(tp.Inh.XCov,[],1),'linewidth',2,'color',cols(4,:))
lgd=legend([h1 h2],'Excitatory', 'Inhibitory');
lgd.FontSize=16;
title('Classification accuracy','FontSize',16)
xlim([minutes_range(1) minutes_range(end)])
xlabel('Simulation time [min]')
ylim([0.6,1])
ylabel('accuracy')
%% Once we'll have maps at all minutes
for minIdx = 1:numel(minutes_range)
    for meIdx = 1:numel(measures)
        fract_error(minIdx,meIdx)=abs(sum(sum(CM.(measures{meIdx}).EE))-sum(sum(CM.gtConn.EE)));
    end
end


