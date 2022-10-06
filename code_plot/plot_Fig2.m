%% Script to reproduce plots in Figure 2

clear all; close all;

% Load simulation results
PATH = ['path_to_main_directory\results\dynamic_connectivity\']; % path to directory where the dynamic connectivity results have been saved

params.doOverlap = 1;
params.nRep = 5;
minutes_simulation = 180;
files = dir([PATH,'dynamicConnInference_rep',num2str(repIdx),'_180min_varyWindow_1_DelConsist*']);
fname = files.name;

load([PATH,fname],'gtConn','gtDelay','s_time','minutes_range','Ne','M');

legend_labels=cell(1,numel(minutes_simulation)+1);
for i = 1:numel(minutes_simulation)
    legend_labels{i}=[num2str(minutes_simulation(i)),'min'];
end
legend_labels{i+1}='null hyp.';

%% Plot conn matrix
%colormap([1,1,1;colormap])
cm = [0.8,0.8,0.8;colormap];
cm = cm(1:end-1,:);

gtDelay(isnan(gtDelay)) = 0;

figure()
imagesc(gtDelay)
colormap(cm)
title('Connectivity matrix with synaptic delays','FontSize',15)
c=colorbar();
c.Label.String='Delay [ms]';
xlabel('Receiver')
ylabel('Emitter')
xline(80.5,'r--','linewidth',2)
yline(80.5,'r--','linewidth',2)

%% Synaptic weights examples
time_min = [1:600];
time_min = time_min/60;

figure()
hold on
plot(time_min,squeeze(s_time(1,1,1:600)),'linewidth',2)
plot(time_min,squeeze(s_time(10,2,1:600)),'linewidth',2)
plot(time_min,squeeze(s_time(Ne+1,1,1:600)),'linewidth',2)
xlabel('time [min]')
ylabel('syn. weight [a.u.]')
legend('Excit. 1','Excit. 2','Inhib. 1')
title('Example synaptic weights temporal evolution')
ylim([-7,10])

%% Synaptic weights autocorrelation

autCorr = zeros(Ne*M,size(s_time,3)*2-1);

for i = 1:Ne
    for j = 1:M
        idx = (i-1)*M+j;
        autCorr(idx,:) = xcov(squeeze(s_time(i,j,:)),'normalized');
        %autCorr(idx,:) = autocorr(squeeze(s_time(i,j,:)));
    end
end

figure()
x_vals = -(size(s_time,3)-1):(size(s_time,3)-1);
x_vals = x_vals/60;
[~,maxIdx] = max(mean(autCorr,1));
x_pos = [ceil(size(autCorr,2)/2):size(autCorr,2)];
[minValue,closestIndex] = min(abs(mean(autCorr(:,x_pos),1)-max(mean(autCorr,1))/2));
plot(x_vals,mean(autCorr,1),'linewidth',2)
shadedErrorBar(x_vals,mean(autCorr,1),std(autCorr,[],1)/sqrt(size(autCorr,1)))
xline(closestIndex/60+ceil(size(autCorr,2)/(2*60))-maxIdx/(60),'linewidth',2)
text((5)*(closestIndex/60+ceil(size(autCorr,2)/(2*60))-maxIdx/(60)),0.4,[num2str(closestIndex+ceil(size(autCorr,2)/2)-maxIdx),'s'],'FontSize',14)
xlim([-(size(s_time,3)-1)/(10*60),(size(s_time,3)-1)/(10*60)])
ylim([-0.2,1.1])
ylabel('Corr.')
xlabel('Lag [min]')
title('Mean synaptic weights autocorrelation')