
%% extract the data- three conditions
clear all
clc
%--------------------------------------------------------------------Condition 1
data1_name = dir ('F:\Fusion\Face-sequence\EEG\Pre_processing\Conditions\Condition1_Neutral\*.set');
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\Fusion\Face-sequence\EEG\Pre_processing\Conditions\Condition1_Neutral');

%initialization 
Neutral_data = zeros (EEG.nbchan,EEG.pnts,length(data1_name)); %channels*data_points*sub_numbs
Neutral_trials = zeros (1,length(data1_name));
Happy_data = zeros (EEG.nbchan,EEG.pnts,length(data1_name));
Happy_trials = zeros (1,length(data1_name));
N2H_data = zeros (EEG.nbchan,EEG.pnts,length(data1_name));
N2H_trials = zeros (1,length(data1_name));
H2N_data = zeros (EEG.nbchan,EEG.pnts,length(data1_name));
H2N_trials = zeros (1,length(data1_name));

for ii = 1:length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\Fusion\Face-sequence\EEG\Pre_processing\Conditions\Condition1_Neutral\');
    Neutral_trials(1,ii) = size (EEG.CSD,3);
    Neutral_data(:,:,ii) = squeeze (nanmean (EEG.CSD,3));
end
%--------------------------------------------------------------------Condition 2
data2_name = dir ('F:\Fusion\Face-sequence\EEG\Pre_processing\Conditions\Condition2_Happy\*.set');
for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\Fusion\Face-sequence\EEG\Pre_processing\Conditions\Condition2_Happy\');
    Happy_trials(1,ii) = size (EEG.CSD,3);
    Happy_data(:,:,ii) = squeeze (nanmean (EEG.CSD,3));
end
%--------------------------------------------------------------------Condition 3
data3_name = dir ('F:\Fusion\Face-sequence\EEG\Pre_processing\Conditions\Condition3_N2H\*.set');
for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\Fusion\Face-sequence\EEG\Pre_processing\Conditions\Condition3_N2H\');
    N2H_trials(1,ii) = size (EEG.CSD,3);
    N2H_data(:,:,ii) = squeeze (nanmean (EEG.CSD,3));
end
%--------------------------------------------------------------------Condition 4
data4_name = dir ('F:\Fusion\Face-sequence\EEG\Pre_processing\Conditions\Condition4_H2N\*.set');
for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\Fusion\Face-sequence\EEG\Pre_processing\Conditions\Condition4_H2N\');
    H2N_trials(1,ii) = size (EEG.CSD,3);
    H2N_data(:,:,ii) = squeeze (nanmean (EEG.CSD,3));
end
time_range = EEG.times;

save ssVEP_ERP   Neutral_data Neutral_trials Happy_data Happy_trials N2H_data N2H_trials H2N_data H2N_trials time_range EEG

%% --------------------------plot an ERP of each channels of ssVEP
clear all
clc
close all
load ssVEP_ERP
figure (1), clf,  
set (gcf,'color','w')
for chani = 1:16;
    subplot (4, 4,chani)
    plot (time_range, squeeze (nanmean(Neutral_data(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (time_range, squeeze (nanmean(Happy_data(chani,:,:),3)),'--','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (time_range, squeeze (nanmean(N2H_data(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (time_range, squeeze (nanmean(H2N_data(chani,:,:),3)),'-r','linewidth',2.5);
    hold on
    line ([-500,3000],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-5,10],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-200 2300],'ylim',[-6, 6],'ydir','reverse','linewidth',2);
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

figure (2), clf
set (gcf,'color','w')
for chani = 17:32;
    subplot (4, 4,chani-16)
    plot (time_range, squeeze (nanmean(Neutral_data(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (time_range, squeeze (nanmean(Happy_data(chani,:,:),3)),'--','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (time_range, squeeze (nanmean(N2H_data(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (time_range, squeeze (nanmean(H2N_data(chani,:,:),3)),'-r','linewidth',2.5);
    hold on
    line ([-200,3000],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-200 2300],'ylim',[-6,6],'ydir','reverse','linewidth',2);
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

figure (3), clf
set (gcf,'color','w')
for chani = 33:48;
    subplot (4, 4,chani-32)
    plot (time_range, squeeze (nanmean(Neutral_data(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (time_range, squeeze (nanmean(Happy_data(chani,:,:),3)),'--','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (time_range, squeeze (nanmean(N2H_data(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (time_range, squeeze (nanmean(H2N_data(chani,:,:),3)),'-r','linewidth',2.5);
    hold on
    line ([-200,3000],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-200 2300],'ylim',[-6,6],'ydir','reverse','linewidth',2);
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

figure (4), clf
set (gcf,'color','w')
for chani = 49:64;
    subplot (4, 4,chani-48)
    plot (time_range, squeeze (nanmean(Neutral_data(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (time_range, squeeze (nanmean(Happy_data(chani,:,:),3)),'--','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (time_range, squeeze (nanmean(N2H_data(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (time_range, squeeze (nanmean(H2N_data(chani,:,:),3)),'-r','linewidth',2.5);
    hold on
    line ([-200,2000],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-200 2300],'ylim',[-6,6],'ydir','reverse','linewidth',2);
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

%% --------------------------Plot an averaged ssVEP ERP of Occipital electrodes for each condition
% 27:30 O1 IZ OZ POZ 64 O2
figure, clf
set (gcf,'color','w')
    plot (time_range, squeeze (nanmean(nanmean(Neutral_data([27:30,64],:,:),3),1)),'-k','linewidth',3);
    hold on
    plot (time_range, squeeze (nanmean(nanmean(N2H_data([27:30,64],:,:),3),1)),'-b','linewidth',3);
    hold on
    plot (time_range, squeeze (nanmean(nanmean(Happy_data([27:30,64],:,:),3),1)),'-','color',[.5,.5,.5],'linewidth',3);
    hold on
    plot (time_range, squeeze (nanmean(nanmean(H2N_data([27:30,64],:,:),3),1)),'-r','linewidth',3);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-15,12],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-200 2300],'ylim',[-12,12],'ytick',-10:5:15,'ydir','reverse','linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    title ('Occipital','FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
%     ylabel (['Amplitude(','\mu','V)'],'FontSize',16,'fontweight','bold','fontname','arial black')
% legend Neutral  N2H Happy H2N
%% --------------------Plot an averaged ssVEP ERP (across all conditions and all subs)of Occipital electrodes
% 27:30 O1 IZ OZ POZ 64 O2
figure, clf
 erp = zeros (4,750);
 erp(1,:) = squeeze (nanmean(nanmean(Neutral_data([27:30,64],:,:),3),1));
 erp(2,:) = squeeze (nanmean(nanmean(N2H_data([27:30,64],:,:),3),1));
 erp(3,:)= squeeze (nanmean(nanmean(Happy_data([27:30,64],:,:),3),1));
 erp(4,:) = squeeze (nanmean(nanmean(H2N_data([27:30,64],:,:),3),1));
set (gcf,'color','w')
    plot (time_range, nanmean (erp),'-k','linewidth',3);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-15,15],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-200 2300],'ylim',[-12,12],'ytick',-15:5:15,'ydir','reverse','linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    title ('Occipital','FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel (['Amplitude(','\mu','V)'],'FontSize',28,'fontweight','bold','fontname','arial black')
%% --------------------------Plot an averaged ssVEP ERP of Frontal electrodes for each condition

figure, clf
set (gcf,'color','w')
    plot (time_range, squeeze (nanmean(nanmean(Neutral_data([1 3 33 34 36 37],:,:),3),1)),'-k','linewidth',3);
    hold on
    plot (time_range, squeeze (nanmean(nanmean(Happy_data([1 3 33 34 36 37],:,:),3),1)),'-','color',[.5,.5,.5],'linewidth',3);
    hold on
    plot (time_range, squeeze (nanmean(nanmean(N2H_data([1 3 33 34 36 37],:,:),3),1)),'-b','linewidth',3);
    hold on
    plot (time_range, squeeze (nanmean(nanmean(H2N_data([1 3 33 34 36 37],:,:),3),1)),'-r','linewidth',3);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-15,15],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-200 2300],'ylim',[-12,12],'ytick',-10:5:10,'ydir','reverse','linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    title ('Frontal','FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
%     ylabel (['Amplitude(','\mu','V)'],'FontSize',16,'fontweight','bold','fontname','arial black')
% legend Neutral Happy N2H H2N

%% ---------------Plot an averaged ssVEP ERP of Frontal electrodes across all conditions and all subjects
% 1 FP1 33 FPZ 34 FP2
figure, clf
 erp = zeros (4,750);
 erp(1,:) = squeeze (nanmean(nanmean(Neutral_data([1 3 33 34 36 37],:,:),3),1));
 erp(2,:) = squeeze (nanmean(nanmean(N2H_data([1 3 33 34 36 37],:,:),3),1));
 erp(3,:)= squeeze (nanmean(nanmean(Happy_data([1 3 33 34 36 37],:,:),3),1));
 erp(4,:) = squeeze (nanmean(nanmean(H2N_data([1 3 33 34 36 37],:,:),3),1));

set (gcf,'color','w')
    plot (time_range, nanmean (erp),'-k','linewidth',3);
    hold on
    line ([-200,3500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-18,12],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-200 2300],'ylim',[-12,12],'ytick',-15:5:15,'ydir','reverse','linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    title ('Frontal','FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel (['Amplitude(','\mu','V)'],'FontSize',28,'fontweight','bold','fontname','arial black')
    
%% --------------------------Plot left,center,right Frontal;left,centre,right Occipatal

% figure (7), clf
% set (gcf,'color','w')
% subplot(2,3,1)
%     plot (time_range, squeeze (nanmean(nanmean(Neutral_data(1:7,:,:),3),1)),'--k','linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(Happy_data(1:7,:,:),3),1)),'--','color',[.5,.5,.5],'linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(N2H_data(1:7,:,:),3),1)),'-b','linewidth',2.5');
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(H2N_data(1:7,:,:),3),1)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle',':');
%     hold on
%     line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle',':');
%     set (gca,'xlim',[-200 2500],'ylim',[-5,8],'ytick',-5:5:10,'ydir','reverse','linewidth',2);
%     title ('Left Frontal')
% subplot(2,3,2)
%     plot (time_range, squeeze (nanmean(nanmean(Neutral_data([33,37,38],:,:),3),1)),'--k','linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(Happy_data([33,37,38],:,:),3),1)),'--','color',[.5,.5,.5],'linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(N2H_data([33,37,38],:,:),3),1)),'-b','linewidth',2.5');
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(H2N_data([33,37,38],:,:),3),1)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle',':');
%     hold on
%     line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle',':');
%     set (gca,'xlim',[-200 2500],'ylim',[-5,8],'ytick',-5:5:10,'ydir','reverse','linewidth',2);
% subplot(2,3,3)
%     plot (time_range, squeeze (nanmean(nanmean(Neutral_data([34:36,39:42],:,:),3),1)),'--k','linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(Happy_data([34:36,39:42],:,:),3),1)),'--','color',[.5,.5,.5],'linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(N2H_data([34:36,39:42],:,:),3),1)),'-b','linewidth',2.5');
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(H2N_data([34:36,39:42],:,:),3),1)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle',':');
%     hold on
%     line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle',':');
%     set (gca,'xlim',[-200 2500],'ylim',[-5,8],'ytick',-5:5:10,'ydir','reverse','linewidth',2);
% subplot(2,3,4)
%     plot (time_range, squeeze (nanmean(nanmean(Neutral_data(25:27,:,:),3),1)),'--k','linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(Happy_data(25:27,:,:),3),1)),'--','color',[.5,.5,.5],'linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(N2H_data(25:27,:,:),3),1)),'-b','linewidth',2.5');
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(H2N_data(25:27,:,:),3),1)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle',':');
%     hold on
%     line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle',':');
%     set (gca,'xlim',[-200 2500],'ylim',[-5,8],'ytick',-5:5:10,'ydir','reverse','linewidth',2);
% subplot(2,3,5)
%     plot (time_range, squeeze (nanmean(nanmean(Neutral_data(28:30,:,:),3),1)),'--k','linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(Happy_data(28:30,:,:),3),1)),'--','color',[.5,.5,.5],'linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(N2H_data(28:30,:,:),3),1)),'-b','linewidth',2.5');
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(H2N_data(28:30,:,:),3),1)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle',':');
%     hold on
%     line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle',':');
%     set (gca,'xlim',[-200 2500],'ylim',[-5,8],'ytick',-5:5:10,'ydir','reverse','linewidth',2);
% subplot(2,3,6)
%     plot (time_range, squeeze (nanmean(nanmean(Neutral_data(62:64,:,:),3),1)),'--k','linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(Happy_data(62:64,:,:),3),1)),'--','color',[.5,.5,.5],'linewidth',2.5);
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(N2H_data(62:64,:,:),3),1)),'-b','linewidth',2.5');
%     hold on
%     plot (time_range, squeeze (nanmean(nanmean(H2N_data(62:64,:,:),3),1)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle',':');
%     hold on
%     line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle',':');
%     set (gca,'xlim',[-200 2500],'ylim',[-5,8],'ytick',-5:5:10,'ydir','reverse','linewidth',2);  
% legend Neutral Happy N2H H2N

%% ------------------------------Plot Global field power 

figure (8), clf
set (gcf,'color','w')
plot (time_range, std((squeeze (nanmean(Neutral_data(:,:,:),3))+ squeeze (nanmean(Happy_data(:,:,:),3))+ squeeze (nanmean(N2H_data(:,:,:),3))+squeeze (nanmean(H2N_data(:,:,:),3)))./3),'-k','linewidth',2.5);
set (gca,'xlim',[-200 2500],'ytick',1:1:4,'linewidth',2,'Box','off');
title ('GFP (Global Field Power)')