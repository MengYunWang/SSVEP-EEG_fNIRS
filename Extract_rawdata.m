%% Extract the data from the preprocessed fnirs data
% Step 1: Use HomerII to preprocess the data
% Step 2: Use CalculatCV to find bad channels and trials
%         then exclude them
% Step 3: Extract the data from each subject procResult.dc according to the
%         markers 's'
% Step 4: Use plot to exam the data
% Created by M-Y Wang
% 02-01-2019
%% Exclude the bad trials and bad channels
    clear all
    clc
    cd  F:\Fusion\Face-sequence\fNIRS\Rawdata
    
    trial_range = 5.1; % according to your segment
    load ('105-sequence.nirs','-mat')
    [trials,~] = find (s);
    procResult.dc (trials(1)-50:trials (1)+trial_range*50,:,[12 13]) = NaN;
    procResult.dc (:,:,[16 17]) = NaN;
    save 105-sequence.nirs 
    
    clear all
    trial_range = 5.1;
    load ('106-sequence.nirs','-mat')
    [trials,~] = find (s);
    procResult.dc (trials(1)-50:trials (1)+trial_range*50,:,8) = NaN;
    save 106-sequence.nirs
    
    clear all
    trial_range = 5.1;
    load ('108-sequence.nirs','-mat')
    procResult.dc (:,:,[12 17]) = NaN;
    save 108-sequence.nirs
    
    clear all
    trial_range = 5.1;
    load ('109-sequence.nirs','-mat')
    procResult.dc (:,:,12) = NaN;    
    save 109-sequence.nirs
    
    clear all
    trial_range = 5.1;
    load ('110-sequence.nirs','-mat')
    [trials,~] = find (s);
    procResult.dc (trials(23)-50:trials (23)+trial_range*50,:,[12 13]) = NaN;
    save 110-sequence.nirs     
    
    clear all
    trial_range = 5.1;
    load ('113-sequence.nirs','-mat') % exclude this subject
    [trials,~] = find (s);
    procResult.dc (trials(91)-50:trials (91)+trial_range*50,:,[15 17]) = NaN;    
    procResult.dc (trials(8)-50:trials (8)+trial_range*50,:,16) = NaN;
    procResult.dc (trials(31)-50:trials (31)+trial_range*50,:,17) = NaN;
    procResult.dc (trials(61)-50:trials (61)+trial_range*50,:,17) = NaN;
    procResult.dc (trials(91)-50:trials (91)+trial_range*50,:,17) = NaN;
    procResult.dc (:,:,[8 12]) = NaN;
    save 113-sequence.nirs   
    
    clear all
    trial_range = 5.1;
    load ('116-sequence.nirs','-mat')
    [trials,~] = find (s);
    procResult.dc (trials(22)-50:trials (22)+trial_range*50,:,[8,12:17]) = NaN; 
    procResult.dc (:,:,12) = NaN; 
    save 116-sequence.nirs
    
    clear all
    trial_range = 5.1;
    load ('117-sequence.nirs','-mat')
    procResult.dc (:,:,[8,12,16,17]) = NaN; 
    save 117-sequence.nirs   
    
    clear all
    trial_range = 5.1;
    load ('119-sequence.nirs','-mat')
    [trials,~] = find (s);
    procResult.dc (trials(32)-50:trials (32)+trial_range*50,:,12) = NaN;
    procResult.dc (trials(66)-50:trials (66)+trial_range*50,:,14) = NaN;
    procResult.dc (trials(86)-50:trials (86)+trial_range*50,:,[13 15]) = NaN;
    procResult.dc (trials(93)-50:trials (93)+trial_range*50,:,[12 15 16]) = NaN;
    procResult.dc (:,:,8) = NaN;
    save 119-sequence.nirs
    
    clear all
    trial_range = 5.1;
    load ('120-sequence.nirs','-mat')
    procResult.dc (:,:,17) = NaN; 
    save 120-sequence.nirs  
    
    clear all
    trial_range = 5.1;
    load ('121-sequence.nirs','-mat')
    procResult.dc (:,:,[16,17]) = NaN; 
    save 121-sequence.nirs
%% Extract data (overlapped 15 trials)

    datatrials = 15;
    data_length = datatrials * trial_range; %length of the data with seconds
    Neutral_data (:,:,:,:,:) = zeros (17,floor((153-data_length)/5),50*data_length+1,3,24); %sub*trials*pnts*hbo*channels
    Happy_data (:,:,:,:,:) = zeros(17,floor((153-data_length)/5),50*data_length+1,3,24); % better to let the pnts ahead of trials!!
    N2H_data (:,:,:,:,:) = zeros (17,floor((153-data_length)/5),50*data_length+1,3,24);
    H2N_data (:,:,:,:,:) = zeros (17,floor((153-data_length)/5),50*data_length+1,3,24);
  
  for subi = [105:112,114:121] %113 excluded,but still subnumber is 17. the 9th row will be empty
      load ([num2str(subi),'-sequence.nirs'],'-mat')
      [trials,conds] = find (s);
      switch conds(2);
          case 1
              for condi = 2:(150-data_length)/5+1;
                  Neutral_data (subi-104,condi-1,:,:,:) = procResult.dc (trials(condi)-50:trials(condi)+round((data_length-1)*50),:,:).*1e6;
              end
          case 2
              for condi = 2:(150-data_length)/5+1;
                  Happy_data (subi-104,condi-1,:,:,:) = procResult.dc (trials(condi)-50:trials(condi)+round((data_length-1)*50),:,:).*1e6;
              end
          case 3
              for condi = 2:(150-data_length)/5+1;
                  N2H_data (subi-104,condi-1,:,:,:) = procResult.dc (trials(condi)-50:trials(condi)+round((data_length-1)*50),:,:).*1e6;
              end
          case 4
              for condi = 2:(150-data_length)/5+1;
                  H2N_data (subi-104,condi-1,:,:,:) = procResult.dc (trials(condi)-50:trials(condi)+round((data_length-1)*50),:,:).*1e6;
              end                 
      end
      
      switch conds(32);
          case 1
              for cond2i = 32:(150-data_length)/5+31;
                  Neutral_data (subi-104,cond2i-31,:,:,:) = procResult.dc (trials(cond2i)-50:trials(cond2i)+round((data_length-1)*50),:,:).*1e6;
              end
          case 2
              for cond2i = 32:(150-data_length)/5+31;
                  Happy_data (subi-104,cond2i-31,:,:,:) = procResult.dc (trials(cond2i)-50:trials(cond2i)+round((data_length-1)*50),:,:).*1e6;
              end
          case 3
              for cond2i = 32:(150-data_length)/5+31;
                  N2H_data (subi-104,cond2i-31,:,:,:) = procResult.dc (trials(cond2i)-50:trials(cond2i)+round((data_length-1)*50),:,:).*1e6;
              end
          case 4
              for cond2i = 32:(150-data_length)/5+31;
                  H2N_data (subi-104,cond2i-31,:,:,:) = procResult.dc (trials(cond2i)-50:trials(cond2i)+round((data_length-1)*50),:,:).*1e6;
              end                 
      end      
      
      switch conds(62);
          case 1
              for cond3i = 62:(150-data_length)/5+61;
                  Neutral_data (subi-104,cond3i-61,:,:,:) = procResult.dc (trials(cond3i)-50:trials(cond3i)+round((data_length-1)*50),:,:).*1e6;
              end
          case 2
              for cond3i = 62:(150-data_length)/5+61;
                  Happy_data (subi-104,cond3i-61,:,:,:) = procResult.dc (trials(cond3i)-50:trials(cond3i)+round((data_length-1)*50),:,:).*1e6;
              end
          case 3
              for cond3i = 62:(150-data_length)/5+61;
                  N2H_data (subi-104,cond3i-61,:,:,:) = procResult.dc (trials(cond3i)-50:trials(cond3i)+round((data_length-1)*50),:,:).*1e6;
              end
          case 4
              for cond3i = 62:(150-data_length)/5+61;
                  H2N_data (subi-104,cond3i-61,:,:,:) = procResult.dc (trials(cond3i)-50:trials(cond3i)+round((data_length-1)*50),:,:).*1e6;
              end                 
      end
      
      switch conds(92);
          case 1
              for cond4i = 92:(150-data_length)/5+91;
                  Neutral_data (subi-104,cond4i-91,:,:,:) = procResult.dc (trials(cond4i)-50:trials(cond4i)+round((data_length-1)*50),:,:).*1e6;
              end
          case 2
              for cond4i = 92:(150-data_length)/5+91;
                  Happy_data (subi-104,cond4i-91,:,:,:) = procResult.dc (trials(cond4i)-50:trials(cond4i)+round((data_length-1)*50),:,:).*1e6;
              end
          case 3
              for cond4i = 92:(150-data_length)/5+91;
                  N2H_data (subi-104,cond4i-91,:,:,:) = procResult.dc (trials(cond4i)-50:trials(cond4i)+round((data_length-1)*50),:,:).*1e6;
              end
          case 4
              for cond4i = 92:(150-data_length)/5+91;
                  H2N_data (subi-104,cond4i-91,:,:,:) = procResult.dc (trials(cond4i)-50:trials(cond4i)+round((data_length-1)*50),:,:).*1e6;
              end                 
      end
 end
timerange = -1:1/50:75.5;

    Neutral_data (9,:,:,:,:) = NaN; 
	Happy_data (9,:,:,:,:) = NaN;
	N2H_data (9,:,:,:,:) = NaN;
    H2N_data (9,:,:,:,:) = NaN;

save  Raw_data timerange Neutral_data Happy_data N2H_data H2N_data 
%%  Plot the time course of the fNIRS data
%     clear all
%     clc
%     close all
%     
%     load Raw_data
    channels = [8,13,12,16,14,15,17];
    
 for Ch = 1:7;
     
     figure(Ch)
     set (gcf,'color','w');
     
     SigHbO1 = squeeze(nanmean (nanmean(Neutral_data (:,:,:,1,channels(Ch)),2)));
     SigHbO2 = squeeze(nanmean (nanmean(Happy_data (:,:,:,1,channels(Ch)),2)));
     SigHbO3 = squeeze(nanmean (nanmean(N2H_data (:,:,:,1,channels(Ch)),2)));
     SigHbO4 = squeeze(nanmean (nanmean(H2N_data (:,:,:,1,channels(Ch)),2)));
     
     plot(timerange,(SigHbO1 + SigHbO2 + SigHbO3 + SigHbO4)./4,'-k','linewidth',3),
     
     %      plot(timerange,SigHbO1,'-k','linewidth',3),
     %      hold on
     %      plot(timerange,SigHbO2,'-','color',[.5, .5,.5],'linewidth',3),
     %      hold on
     %      plot(timerange,SigHbO3,'-b','linewidth',3),
     %      hold on
     %      plot(timerange,SigHbO4,'-r','linewidth',3),
     
     set(gca,'xlim',[-1 30],'xtick',0:5:80,'ylim',[-2e-3 2e-3]);
     set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
     title (['Ch 0', num2str(Ch)],'FontSize',28,'fontweight','bold','fontname','arial black')
%      xlabel ('Time (s)','FontSize',28,'fontweight','bold','fontname','arial black')
%      ylabel (['\Delta ','HbO(micro-mole)'],'FontSize',28,'fontweight','bold','fontname','arial black')
     grid on
 end
%  legend('Neutral','Happy','N2H','H2N')



% %%  Extract data from a single trial-------not good
%     clear all
%     clc
%     cd ('F:\Fusion\Face-sequence\fNIRS\Rawdata');
%     datatrials = 1;
%     data_length = datatrials * trial_range; %length of the data with seconds
%     Neutral_data (:,:,:,:,:) = zeros (17,round(50*data_length+1),30,3,24); %sub*pnts*trials*hbo*channels
%     Happy_data (:,:,:,:,:) = zeros(17,round(50*data_length+1),30,3,24);
%     N2H_data (:,:,:,:,:) = zeros (17,round(50*data_length+1),30,3,24);
%     H2N_data (:,:,:,:,:) = zeros (17,round(50*data_length+1),30,3,24);
% 
%   
%   for subi = [105:112,114:121] %113 excluded,but still subnumber is 17. the 9th row will be empty
%       load ([num2str(subi),'-sequence.nirs'],'-mat')
%       [trials,conds] = find (s);
%       switch conds(2);
%           case 1
%               for condi = 1:30;
%                   Neutral_data (subi-104,:,condi,:,:) = procResult.dc (trials(condi)-50:trials(condi)+round(4.1*50),:,:).*1e6;
%               end
%           case 2
%               for condi = 1:30;
%                   Happy_data (subi-104,:,condi,:,:) = procResult.dc (trials(condi)-50:trials(condi)+round(4.1*50),:,:).*1e6;
%               end
%           case 3
%               for condi = 1:30;
%                   N2H_data (subi-104,:,condi,:,:) = procResult.dc (trials(condi)-50:trials(condi)+round(4.1*50),:,:).*1e6;
%               end
%           case 4
%               for condi = 1:30;
%                   H2N_data (subi-104,:,condi,:,:) = procResult.dc (trials(condi)-50:trials(condi)+round(4.1*50),:,:).*1e6;
%               end                 
%       end
%       
%       switch conds(32);
%           case 1
%               for cond2i = 31:60;
%                   Neutral_data (subi-104,:,cond2i-30,:,:) = procResult.dc (trials(cond2i)-50:trials(cond2i)+round(4.1*50),:,:).*1e6;
%               end
%           case 2
%               for cond2i = 31:60;
%                   Happy_data (subi-104,:,cond2i-30,:,:) = procResult.dc (trials(cond2i)-50:trials(cond2i)+round(4.1*50),:,:).*1e6;
%               end
%           case 3
%               for cond2i = 31:60;
%                   N2H_data (subi-104,:,cond2i-30,:,:) = procResult.dc (trials(cond2i)-50:trials(cond2i)+round(4.1*50),:,:).*1e6;
%               end
%           case 4
%               for cond2i = 31:60;
%                   H2N_data (subi-104,:,cond2i-30,:,:) = procResult.dc (trials(cond2i)-50:trials(cond2i)+round(4.1*50),:,:).*1e6;
%               end                 
%       end      
%       
%       switch conds(62);
%           case 1
%               for cond3i = 61:90;
%                   Neutral_data (subi-104,:,cond3i-60,:,:) = procResult.dc (trials(cond3i)-50:trials(cond3i)+round(4.1*50),:,:).*1e6;
%               end
%           case 2
%               for cond3i =61:90;
%                   Happy_data (subi-104,:,cond3i-60,:,:) = procResult.dc (trials(cond3i)-50:trials(cond3i)+round(4.1*50),:,:).*1e6;
%               end
%           case 3
%               for cond3i = 61:90;
%                   N2H_data (subi-104,:,cond3i-60,:,:) = procResult.dc (trials(cond3i)-50:trials(cond3i)+round(4.1*50),:,:).*1e6;
%               end
%           case 4
%               for cond3i = 61:90;
%                   H2N_data (subi-104,:,cond3i-60,:,:) = procResult.dc (trials(cond3i)-50:trials(cond3i)+round(4.1*50),:,:).*1e6;
%               end                 
%       end
%       
%       switch conds(92);
%           case 1
%               for cond4i = 91:120;
%                   Neutral_data (subi-104,:,cond4i-90,:,:) = procResult.dc (trials(cond4i)-50:trials(cond4i)+ round(4.1*50),:,:).*1e6;
%               end
%           case 2
%               for cond4i = 91:120;
%                   Happy_data (subi-104,:,cond4i-90,:,:) = procResult.dc (trials(cond4i)-50:trials(cond4i)+round(4.1*50),:,:).*1e6;
%               end
%           case 3
%               for cond4i = 91:120;
%                   N2H_data (subi-104,:,cond4i-90,:,:) = procResult.dc (trials(cond4i)-50:trials(cond4i)+round(4.1*50),:,:).*1e6;
%               end
%           case 4
%               for cond4i = 91:120;
%                   H2N_data (subi-104,:,cond4i-90,:,:) = procResult.dc (trials(cond4i)-50:trials(cond4i)+round(4.1*50),:,:).*1e6;
%               end                 
%       end
%  end
% timerange = procResult.tHRF;
% 
%     Neutral_data (9,:,:,:,:) = NaN; 
% 	Happy_data (9,:,:,:,:) = NaN;
% 	N2H_data (9,:,:,:,:) = NaN;
%     H2N_data (9,:,:,:,:) = NaN;
% 
% save  Raw_data_single_trial timerange Neutral_data Happy_data N2H_data H2N_data 
% %%  Plot the time course of the fNIRS data
% 
% 
%  channels = [8,13,12,16,14,15,17];
%  for Ch = 1:7;
%      
%      figure(Ch)
%      set (gcf,'color','w');
%      
% %      SigHbO1 = squeeze(nanmean (nanmean(Neutral_data (:,:,:,1,channels(Ch)),2)));
% %      SigHbO2 = squeeze(nanmean (nanmean(Happy_data (:,:,:,1,channels(Ch)),2)));
% %      SigHbO3 = squeeze(nanmean (nanmean(N2H_data (:,:,:,1,channels(Ch)),2)));
% %      SigHbO4 = squeeze(nanmean (nanmean(H2N_data (:,:,:,1,channels(Ch)),2)));
%      
%      SigHbO1 = squeeze(Neutral_data (1,1,:,1,channels(Ch)));
%      SigHbO2 = squeeze(Happy_data (1,1,:,1,channels(Ch)));
%      SigHbO3 = squeeze(N2H_data (1,1,:,1,channels(Ch)));
%      SigHbO4 = squeeze(H2N_data (1,1,:,1,channels(Ch)));
%      
%      plot(timerange,SigHbO1,'-k','linewidth',3),
%      hold on
%      plot(timerange,SigHbO2,'-','color',[.5, .5,.5],'linewidth',3),
%      hold on
%      plot(timerange,SigHbO3,'-b','linewidth',3),
%      hold on
%      plot(timerange,SigHbO4,'-r','linewidth',3),
%      
% %      set(gca,'xlim',[-1 4],'xtick',-1:1:4,'ylim',[-2e-3,2e-3]);
%           set(gca,'xlim',[-1 4],'xtick',-1:1:4);
% 
%      set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%      title (['Ch', num2str(Ch)],'FontSize',28,'fontweight','bold','fontname','arial black')
% %      xlabel ('Time (s)','FontSize',28,'fontweight','bold','fontname','arial black')
% %      ylabel (['\Delta ','HbO(micro-mole)'],'FontSize',28,'fontweight','bold','fontname','arial black')
%      grid on
%  end
% %  legend('Neutral','Happy','N2H','H2N')