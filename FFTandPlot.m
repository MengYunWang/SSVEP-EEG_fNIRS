% Fourier transform for fnirs data

% created by M-Y Wang

%% fft for overlapped trials

    clear all
    clc
    cd  F:\Fusion\Face-sequence\fNIRS\

    load Raw_data
    channels = [8,13,12,16,14,15,17]; 

    HbO1_raw = squeeze(nanmean(Neutral_data (:,:,:,1,channels),2)); 
    HbO2_raw = squeeze(nanmean(Happy_data (:,:,:,1,channels),2));   
    HbO3_raw = squeeze(nanmean(N2H_data (:,:,:,1,channels),2));     
    HbO4_raw = squeeze(nanmean(H2N_data (:,:,:,1,channels),2));     
 
    HbR1_raw = squeeze(nanmean(Neutral_data (:,:,:,2,channels),2)); 
    HbR2_raw = squeeze(nanmean(Happy_data (:,:,:,2,channels),2));   
    HbR3_raw = squeeze(nanmean(N2H_data (:,:,:,2,channels),2));     
    HbR4_raw = squeeze(nanmean(H2N_data (:,:,:,2,channels),2));     
    
    HbT1_raw = squeeze(nanmean(Neutral_data (:,:,:,3,channels),2)); 
    HbT2_raw = squeeze(nanmean(Happy_data (:,:,:,3,channels),2));   
    HbT3_raw = squeeze(nanmean(N2H_data (:,:,:,3,channels),2));     
    HbT4_raw = squeeze(nanmean(H2N_data (:,:,:,3,channels),2));     
    
    datasave(1).name = {'Neutral_pow_HbO','Happy_pow_HbO','N2H_pow_HbO','H2N_pow_HbO'};    
    datasave(2).name = {'Neutral_pow_HbR','Happy_pow_HbR','N2H_pow_HbR','H2N_pow_HbR'};    
    datasave(3).name = {'Neutral_pow_HbT','Happy_pow_HbT','N2H_pow_HbT','H2N_pow_HbT'}; 
    
    datasave(4).name = {'Neutral_amp_HbO','Happy_amp_HbO','N2H_amp_HbO','H2N_amp_HbO'};    
    datasave(5).name = {'Neutral_amp_HbR','Happy_amp_HbR','N2H_amp_HbR','H2N_amp_HbR'};    
    datasave(6).name = {'Neutral_amp_HbT','Happy_amp_HbT','N2H_amp_HbT','H2N_amp_HbT'};
%----------------------------------------
    samplerate = 50;                    
    data_n = size(HbO1_raw,2);
    data_N = 2^nextpow2(data_n); 
    data_hz = linspace(0,samplerate/2,floor(data_N/2+1));


for condi = 1:4;
    dataX = eval(['fft([','HbO', num2str(condi),'_raw],data_N,2)']);
    data_amp = abs(dataX)./data_N;
    data_amp = data_amp(:,1:data_N/2+1,:);
    data_amp (:,2:end-1,:) = 2.*data_amp(:,2:end-1,:);
    eval ([datasave(4).name{condi}, '= data_amp;']);    
    data_pow = data_amp.^2;
    eval ([datasave(1).name{condi}, '= data_pow;']);
end

for condi = 1:4;
    dataX = eval(['fft([','HbR', num2str(condi),'_raw],data_N,2)']);
    data_amp = abs(dataX)./data_N;
    data_amp = data_amp(:,1:data_N/2+1,:);
    data_amp (:,2:end-1,:) = 2.*data_amp(:,2:end-1,:);
    eval ([datasave(5).name{condi}, '= data_amp;']);      
    data_pow = data_amp.^2;
    eval ([datasave(2).name{condi}, '= data_pow;']);
end

for condi = 1:4;
    dataX = eval(['fft([','HbT', num2str(condi),'_raw],data_N,2)']);
    data_amp = abs(dataX)./data_N;
    data_amp = data_amp(:,1:data_N/2+1,:);
    data_amp (:,2:end-1,:) = 2.*data_amp(:,2:end-1,:);
    eval ([datasave(6).name{condi}, '= data_amp;']);  
    data_pow = data_amp.^2;
    eval ([datasave(3).name{condi}, '= data_pow;']);
end

save fft_several_trials.mat data_hz...
     Neutral_pow_HbO Happy_pow_HbO N2H_pow_HbO H2N_pow_HbO Neutral_amp_HbO Happy_amp_HbO N2H_amp_HbO H2N_amp_HbO...
     Neutral_pow_HbR Happy_pow_HbR N2H_pow_HbR H2N_pow_HbR Neutral_amp_HbR Happy_amp_HbR N2H_amp_HbR H2N_amp_HbR...
     Neutral_pow_HbT Happy_pow_HbT N2H_pow_HbT H2N_pow_HbT Neutral_amp_HbT Happy_amp_HbT N2H_amp_HbT H2N_amp_HbT

 %% Compute the SNR and Z-score
     clear all
     clc
     load fft_several_trials
%---------------------------------------------------HbO 
    dimen = size  (H2N_amp_HbO);
    Neutral_SNR_HbO = zeros (dimen(1),dimen(2),dimen(3));  Neutral_Z_HbO = zeros (dimen(1),dimen(2),dimen(3));
    Happy_SNR_HbO = zeros (dimen(1),dimen(2),dimen(3));    Happy_Z_HbO = zeros (dimen(1),dimen(2),dimen(3));
    N2H_SNR_HbO = zeros (dimen(1),dimen(2),dimen(3));      N2H_Z_HbO = zeros (dimen(1),dimen(2),dimen(3));
    H2N_SNR_HbO = zeros (dimen(1),dimen(2),dimen(3));      H2N_Z_HbO = zeros (dimen(1),dimen(2),dimen(3));
    % 
    for mm = 1:size (Neutral_amp_HbO,1); % subjects
        for ii = 11:length (data_hz)-11; % Hz
            neutral_temp = zeros (18,7); % channel * aside frex
            neutral_temp(1:9,:) = squeeze (Neutral_amp_HbO(mm,ii-10:ii-2,:)); neutral_temp (10:18,:) = squeeze (Neutral_amp_HbO(mm,ii+2:ii+10,:));
            [~,indx1] = max (neutral_temp,[],1); for jj = 1:7; neutral_temp (indx1(jj),jj) = NaN; end % exclude two max values
            [~,indx2] = max (neutral_temp,[],1); for jj = 1:7; neutral_temp (indx2(jj),jj) = NaN; end
            Neutral_SNR_HbO (mm,ii,:) = squeeze(Neutral_amp_HbO(mm,ii,:)) ./ nanmean(neutral_temp)';
            Neutral_Z_HbO (mm,ii,:)= (squeeze(Neutral_amp_HbO(mm,ii,:)) - nanmean(neutral_temp)') ./ nanstd(neutral_temp,[],1)';
            
            happy_temp = zeros (18,7); % channel * aside frex
            happy_temp(1:9,:) = squeeze (Happy_amp_HbO(mm,ii-10:ii-2,:)); happy_temp (10:18,:) = squeeze (Happy_amp_HbO(mm,ii+2:ii+10,:));
            [~,indx1] = max (happy_temp,[],1); for jj = 1:7; happy_temp (indx1(jj),jj) = NaN; end
            [~,indx2] = max (happy_temp,[],1); for jj = 1:7; happy_temp (indx2(jj),jj) = NaN; end
            Happy_SNR_HbO (mm,ii,:) = squeeze(Happy_amp_HbO(mm,ii,:)) ./ nanmean(happy_temp)';
            Happy_Z_HbO (mm,ii,:)= (squeeze(Happy_amp_HbO(mm,ii,:)) - nanmean(happy_temp)') ./ nanstd(happy_temp,[],1)';
            
            N2H_temp = zeros (18,7); % channel * aside frex
            N2H_temp(1:9,:) = squeeze (N2H_amp_HbO(mm,ii-10:ii-2,:)); N2H_temp (10:18,:) = squeeze (N2H_amp_HbO(mm,ii+2:ii+10,:));
            [~,indx1] = max (N2H_temp,[],1); for jj = 1:7; N2H_temp (indx1(jj),jj) = NaN; end
            [~,indx2] = max (N2H_temp,[],1); for jj = 1:7; N2H_temp (indx2(jj),jj) = NaN; end
            N2H_SNR_HbO (mm,ii,:) = squeeze(N2H_amp_HbO(mm,ii,:)) ./ nanmean(N2H_temp)';
            N2H_Z_HbO (mm,ii,:)= (squeeze(N2H_amp_HbO(mm,ii,:)) - nanmean(N2H_temp)') ./ nanstd(N2H_temp,[],1)';
            
            H2N_temp = zeros (18,7); % channel * aside frex
            H2N_temp(1:9,:) = squeeze (H2N_amp_HbO(mm,ii-10:ii-2,:)); H2N_temp (10:18,:) = squeeze (H2N_amp_HbO(mm,ii+2:ii+10,:));
            [~,indx1] = max (H2N_temp,[],1); for jj = 1:7; H2N_temp (indx1(jj),jj) = NaN; end
            [~,indx2] = max (H2N_temp,[],1); for jj = 1:7; H2N_temp (indx2(jj),jj) = NaN; end
            H2N_SNR_HbO (mm,ii,:) = squeeze(H2N_amp_HbO(mm,ii,:)) ./ nanmean(H2N_temp)';
            H2N_Z_HbO (mm,ii,:)= (squeeze(H2N_amp_HbO(mm,ii,:)) - nanmean(H2N_temp)') ./ nanstd(H2N_temp,[],1)';
        end
    end
    
    %---------------------------------------------------HbR 
    dimen = size  (H2N_amp_HbR);
    Neutral_SNR_HbR = zeros (dimen(1),dimen(2),dimen(3));  Neutral_Z_HbR = zeros (dimen(1),dimen(2),dimen(3));
    Happy_SNR_HbR = zeros (dimen(1),dimen(2),dimen(3));    Happy_Z_HbR = zeros (dimen(1),dimen(2),dimen(3));
    N2H_SNR_HbR = zeros (dimen(1),dimen(2),dimen(3));      N2H_Z_HbR = zeros (dimen(1),dimen(2),dimen(3));
    H2N_SNR_HbR = zeros (dimen(1),dimen(2),dimen(3));      H2N_Z_HbR = zeros (dimen(1),dimen(2),dimen(3));
    % 
    for mm = 1:size (Neutral_amp_HbR,1); % subjects
        for ii = 11:length (data_hz)-11; % Hz
            neutral_temp = zeros (18,7); % channel * aside frex
            neutral_temp(1:9,:) = squeeze (Neutral_amp_HbR(mm,ii-10:ii-2,:)); neutral_temp (10:18,:) = squeeze (Neutral_amp_HbR(mm,ii+2:ii+10,:));
            [~,indx1] = max (neutral_temp,[],1); for jj = 1:7; neutral_temp (indx1(jj),jj) = NaN; end % exclude two max values
            [~,indx2] = max (neutral_temp,[],1); for jj = 1:7; neutral_temp (indx2(jj),jj) = NaN; end
            Neutral_SNR_HbR (mm,ii,:) = squeeze(Neutral_amp_HbR(mm,ii,:)) ./ nanmean(neutral_temp)';
            Neutral_Z_HbR (mm,ii,:)= (squeeze(Neutral_amp_HbR(mm,ii,:)) - nanmean(neutral_temp)') ./ nanstd(neutral_temp,[],1)';
            
            happy_temp = zeros (18,7); % channel * aside frex
            happy_temp(1:9,:) = squeeze (Happy_amp_HbR(mm,ii-10:ii-2,:)); happy_temp (10:18,:) = squeeze (Happy_amp_HbR(mm,ii+2:ii+10,:));
            [~,indx1] = max (happy_temp,[],1); for jj = 1:7; happy_temp (indx1(jj),jj) = NaN; end
            [~,indx2] = max (happy_temp,[],1); for jj = 1:7; happy_temp (indx2(jj),jj) = NaN; end
            Happy_SNR_HbR (mm,ii,:) = squeeze(Happy_amp_HbR(mm,ii,:)) ./ nanmean(happy_temp)';
            Happy_Z_HbR (mm,ii,:)= (squeeze(Happy_amp_HbR(mm,ii,:)) - nanmean(happy_temp)') ./ nanstd(happy_temp,[],1)';
            
            N2H_temp = zeros (18,7); % channel * aside frex
            N2H_temp(1:9,:) = squeeze (N2H_amp_HbR(mm,ii-10:ii-2,:)); N2H_temp (10:18,:) = squeeze (N2H_amp_HbR(mm,ii+2:ii+10,:));
            [~,indx1] = max (N2H_temp,[],1); for jj = 1:7; N2H_temp (indx1(jj),jj) = NaN; end
            [~,indx2] = max (N2H_temp,[],1); for jj = 1:7; N2H_temp (indx2(jj),jj) = NaN; end
            N2H_SNR_HbR (mm,ii,:) = squeeze(N2H_amp_HbR(mm,ii,:)) ./ nanmean(N2H_temp)';
            N2H_Z_HbR (mm,ii,:)= (squeeze(N2H_amp_HbR(mm,ii,:)) - nanmean(N2H_temp)') ./ nanstd(N2H_temp,[],1)';
            
            H2N_temp = zeros (18,7); % channel * aside frex
            H2N_temp(1:9,:) = squeeze (H2N_amp_HbR(mm,ii-10:ii-2,:)); H2N_temp (10:18,:) = squeeze (H2N_amp_HbR(mm,ii+2:ii+10,:));
            [~,indx1] = max (H2N_temp,[],1); for jj = 1:7; H2N_temp (indx1(jj),jj) = NaN; end
            [~,indx2] = max (H2N_temp,[],1); for jj = 1:7; H2N_temp (indx2(jj),jj) = NaN; end
            H2N_SNR_HbR (mm,ii,:) = squeeze(H2N_amp_HbR(mm,ii,:)) ./ nanmean(H2N_temp)';
            H2N_Z_HbR (mm,ii,:)= (squeeze(H2N_amp_HbR(mm,ii,:)) - nanmean(H2N_temp)') ./ nanstd(H2N_temp,[],1)';
        end
    end
    
    
    %---------------------------------------------------HbT 
    dimen = size  (H2N_amp_HbT);
    Neutral_SNR_HbT = zeros (dimen(1),dimen(2),dimen(3));  Neutral_Z_HbT = zeros (dimen(1),dimen(2),dimen(3));
    Happy_SNR_HbT = zeros (dimen(1),dimen(2),dimen(3));    Happy_Z_HbT = zeros (dimen(1),dimen(2),dimen(3));
    N2H_SNR_HbT = zeros (dimen(1),dimen(2),dimen(3));      N2H_Z_HbT = zeros (dimen(1),dimen(2),dimen(3));
    H2N_SNR_HbT = zeros (dimen(1),dimen(2),dimen(3));      H2N_Z_HbT = zeros (dimen(1),dimen(2),dimen(3));
    % 
    for mm = 1:size (Neutral_amp_HbT,1); % subjects
        for ii = 11:length (data_hz)-11; % Hz
            neutral_temp = zeros (18,7); % channel * aside frex
            neutral_temp(1:9,:) = squeeze (Neutral_amp_HbT(mm,ii-10:ii-2,:)); neutral_temp (10:18,:) = squeeze (Neutral_amp_HbT(mm,ii+2:ii+10,:));
            [~,indx1] = max (neutral_temp,[],1); for jj = 1:7; neutral_temp (indx1(jj),jj) = NaN; end % exclude two max values
            [~,indx2] = max (neutral_temp,[],1); for jj = 1:7; neutral_temp (indx2(jj),jj) = NaN; end
            Neutral_SNR_HbT (mm,ii,:) = squeeze(Neutral_amp_HbT(mm,ii,:)) ./ nanmean(neutral_temp)';
            Neutral_Z_HbT (mm,ii,:)= (squeeze(Neutral_amp_HbT(mm,ii,:)) - nanmean(neutral_temp)') ./ nanstd(neutral_temp,[],1)';
            
            happy_temp = zeros (18,7); % channel * aside frex
            happy_temp(1:9,:) = squeeze (Happy_amp_HbT(mm,ii-10:ii-2,:)); happy_temp (10:18,:) = squeeze (Happy_amp_HbT(mm,ii+2:ii+10,:));
            [~,indx1] = max (happy_temp,[],1); for jj = 1:7; happy_temp (indx1(jj),jj) = NaN; end
            [~,indx2] = max (happy_temp,[],1); for jj = 1:7; happy_temp (indx2(jj),jj) = NaN; end
            Happy_SNR_HbT (mm,ii,:) = squeeze(Happy_amp_HbT(mm,ii,:)) ./ nanmean(happy_temp)';
            Happy_Z_HbT (mm,ii,:)= (squeeze(Happy_amp_HbT(mm,ii,:)) - nanmean(happy_temp)') ./ nanstd(happy_temp,[],1)';
            
            N2H_temp = zeros (18,7); % channel * aside frex
            N2H_temp(1:9,:) = squeeze (N2H_amp_HbT(mm,ii-10:ii-2,:)); N2H_temp (10:18,:) = squeeze (N2H_amp_HbT(mm,ii+2:ii+10,:));
            [~,indx1] = max (N2H_temp,[],1); for jj = 1:7; N2H_temp (indx1(jj),jj) = NaN; end
            [~,indx2] = max (N2H_temp,[],1); for jj = 1:7; N2H_temp (indx2(jj),jj) = NaN; end
            N2H_SNR_HbT (mm,ii,:) = squeeze(N2H_amp_HbT(mm,ii,:)) ./ nanmean(N2H_temp)';
            N2H_Z_HbT (mm,ii,:)= (squeeze(N2H_amp_HbT(mm,ii,:)) - nanmean(N2H_temp)') ./ nanstd(N2H_temp,[],1)';
            
            H2N_temp = zeros (18,7); % channel * aside frex
            H2N_temp(1:9,:) = squeeze (H2N_amp_HbT(mm,ii-10:ii-2,:)); H2N_temp (10:18,:) = squeeze (H2N_amp_HbT(mm,ii+2:ii+10,:));
            [~,indx1] = max (H2N_temp,[],1); for jj = 1:7; H2N_temp (indx1(jj),jj) = NaN; end
            [~,indx2] = max (H2N_temp,[],1); for jj = 1:7; H2N_temp (indx2(jj),jj) = NaN; end
            H2N_SNR_HbT (mm,ii,:) = squeeze(H2N_amp_HbT(mm,ii,:)) ./ nanmean(H2N_temp)';
            H2N_Z_HbT (mm,ii,:)= (squeeze(H2N_amp_HbT(mm,ii,:)) - nanmean(H2N_temp)') ./ nanstd(H2N_temp,[],1)';
        end
    end
save fft_SNRZ  Neutral_SNR_HbT  Neutral_Z_HbT  Happy_SNR_HbT Happy_Z_HbT  N2H_SNR_HbT N2H_Z_HbT H2N_SNR_HbT H2N_Z_HbT data_hz...
               Neutral_SNR_HbR  Neutral_Z_HbR  Happy_SNR_HbR Happy_Z_HbR  N2H_SNR_HbR N2H_Z_HbR H2N_SNR_HbR H2N_Z_HbR...
               Neutral_SNR_HbO  Neutral_Z_HbO  Happy_SNR_HbO Happy_Z_HbO  N2H_SNR_HbO N2H_Z_HbO H2N_SNR_HbO H2N_Z_HbO

%% Plot SNR for each fNRIS channels
    close all
    clear all
    clc

    load fft_SNRZ.mat
%----------------------------------HbO
for chani = 1:7;
    figure, clf, 
    set (gcf,'color','w')
         
%     plot (data_hz, squeeze (nanmean(Neutral_Z_HbO(:,:,chani))),'--k','linewidth',2.5);
%     hold on
%     plot (data_hz, squeeze (nanmean(Happy_Z_HbO(:,:,chani))),'-.','color',[.5,.5,.5],'linewidth',2.5);
%     hold on
%     plot (data_hz, squeeze (nanmean(N2H_Z_HbO(:,:,chani))),'-b','linewidth',2.5');
%     hold on
%     plot (data_hz, squeeze (nanmean(H2N_Z_HbO(:,:,chani))),'-r','linewidth',2.5);
%     hold on
    plot (data_hz, (squeeze (nanmean(Neutral_SNR_HbO(:,:,chani))) + squeeze (nanmean(Happy_SNR_HbO(:,:,chani)))...
        + squeeze (nanmean(N2H_SNR_HbO(:,:,chani))) + squeeze (nanmean(H2N_SNR_HbO(:,:,chani))))./4,'-k','linewidth',2.5);
%        set (gca,'xlim',[0.15 0.25],'xtick',0:0.2:0.6,'linewidth',2.5) 
    set (gca,'xlim',[0.15 0.45],'xtick',0:0.2:0.6,'ylim',[0 15],'ytick',0:5:15,'linewidth',2.5)
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel ('SNR','FontSize',28,'fontweight','bold','fontname','arial black')
    title (['Ch 0',num2str(chani)],'Fontsize',28,'fontweight','bold','fontname','arial black')
end
% legend Neutral Happy N2H H2N
%----------------------------------HbR
for chani = 1:7;
    figure, clf, 
    set (gcf,'color','w')
         
    plot (data_hz, squeeze (nanmean(Neutral_Z_HbR(:,:,chani))),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(Happy_Z_HbR(:,:,chani))),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(N2H_Z_HbR(:,:,chani))),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (nanmean(H2N_Z_HbR(:,:,chani))),'-r','linewidth',2.5);
       set (gca,'xlim',[0.15 0.25],'xtick',0:0.2:0.6,'linewidth',2.5) 
%     set (gca,'xlim',[0.15 0.45],'xtick',0:0.2:0.6,'ylim',[0 15],'ytick',0:5:15,'linewidth',2.5)
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('SNR','FontSize',28,'fontweight','bold','fontname','arial black')
    title (['Ch',num2str(chani)],'Fontsize',28,'fontweight','bold','fontname','arial black')
end

%----------------------------------HbT
for chani = 1:7;
    figure, clf, 
    set (gcf,'color','w')
         
    plot (data_hz, squeeze (nanmean(Neutral_Z_HbT(:,:,chani))),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(Happy_Z_HbT(:,:,chani))),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(N2H_Z_HbT(:,:,chani))),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (nanmean(H2N_Z_HbT(:,:,chani))),'-r','linewidth',2.5);
       set (gca,'xlim',[0.15 0.25],'xtick',0:0.2:0.6,'linewidth',2.5) 
%     set (gca,'xlim',[0.15 0.45],'xtick',0:0.2:0.6,'ylim',[0 15],'ytick',0:5:15,'linewidth',2.5)
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('SNR','FontSize',28,'fontweight','bold','fontname','arial black')
    title (['Ch',num2str(chani)],'Fontsize',28,'fontweight','bold','fontname','arial black')
end

 %% Plot SNR across all fNIRS channels
 % across all conditions and subjects
     figure, clf
     SNR_all = zeros (4,2049,7);
     SNR_all(1,:,:) = squeeze (nanmean(Neutral_SNR_HbO(:,:,:)));
     SNR_all(2,:,:) = squeeze (nanmean(Happy_SNR_HbO(:,:,:)));
     SNR_all(3,:,:)= squeeze (nanmean(N2H_SNR_HbO(:,:,:)));
     SNR_all(4,:,:) = squeeze (nanmean(H2N_SNR_HbO(:,:,:)));

%      set (gcf,'color','w')
%      plot (data_hz,nanmean(SNR_all),'-k','linewidth',2.5);
%     set (gca,'xlim',[0.15 0.45],'xtick',0:0.2:0.6,'ylim',[0 15],'ytick',0:5:15,'linewidth',2.5)
%     set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel ('SNR','FontSize',28,'fontweight','bold','fontname','arial black')
%     title (['Ch',num2str(chani)],'Fontsize',28,'fontweight','bold','fontname','arial black')

% for left brain
    figure, clf
    
    set (gcf,'color','w')
    plot (data_hz, squeeze (nanmean(SNR_all(1,:,[3 4 7]),3)),'-k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(SNR_all(2,:,[3 4 7]),3)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(SNR_all(3,:,[3 4 7]),3)),'-b','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(SNR_all(4,:,[3 4 7]),3)),'-r','linewidth',2.5);
    set (gca,'xlim',[0.14 0.26],'xtick',0:0.05:0.6,'ylim',[-0.75 12],'ytick',0:5:15,'linewidth',2.5)
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel ('SNR','FontSize',28,'fontweight','bold','fontname','arial black')
    title ('Left','Fontsize',28,'fontweight','bold','fontname','arial black')
%     legend Neutral Happy N2H H2N
 
% for right brain
    figure, clf
    
    set (gcf,'color','w')
    plot (data_hz, squeeze (nanmean(SNR_all(1,:,[1 2 5]),3)),'-k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(SNR_all(2,:,[1 2 5]),3)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(SNR_all(3,:,[1 2 5]),3)),'-b','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(SNR_all(4,:,[1 2 5]),3)),'-r','linewidth',2.5);
    set (gca,'xlim',[0.12 0.28],'xtick',0:0.05:0.6,'ylim',[0.75 15],'ytick',0:5:15,'linewidth',2.5)
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    title ('Right','Fontsize',28,'fontweight','bold','fontname','arial black')

%% fft for single trial
%     clear all
%     close all
%     cd  F:\Fusion\Face-sequence\fNIRS\
% 
%     load Raw_data_single_trial
%     channels = [8,13,12,16,14,15,17]; 
% 
%     HbO1_raw = squeeze(Neutral_data (:,:,:,1,channels));   HbO1_raw = reshape (HbO1_raw,17,30*256,7);  % ×ª»»ÓÐ´íÎó!!!!
%     HbO2_raw = squeeze(Happy_data (:,:,:,1,channels));     HbO2_raw = reshape (HbO2_raw,17,30*256,7);
%     HbO3_raw = squeeze(N2H_data (:,:,:,1,channels));       HbO3_raw = reshape (HbO3_raw,17,30*256,7);
%     HbO4_raw = squeeze(H2N_data (:,:,:,1,channels));       HbO4_raw = reshape (HbO4_raw,17,30*256,7);
%  
%     HbR1_raw = squeeze(Neutral_data (:,:,:,2,channels)); HbR1_raw = reshape (HbR1_raw,17,30*256,7);
%     HbR2_raw = squeeze(Happy_data (:,:,:,2,channels));   HbR2_raw = reshape (HbR2_raw,17,30*256,7);
%     HbR3_raw = squeeze(N2H_data (:,:,:,2,channels));     HbR3_raw = reshape (HbR3_raw,17,30*256,7);
%     HbR4_raw = squeeze(H2N_data (:,:,:,2,channels));     HbR4_raw = reshape (HbR4_raw,17,30*256,7);
%     
%     HbT1_raw = squeeze(Neutral_data (:,:,:,3,channels)); HbT1_raw = reshape (HbT1_raw,17,30*256,7);
%     HbT2_raw = squeeze(Happy_data (:,:,:,3,channels));   HbT2_raw = reshape (HbT2_raw,17,30*256,7);
%     HbT3_raw = squeeze(N2H_data (:,:,:,3,channels));     HbT3_raw = reshape (HbT3_raw,17,30*256,7);
%     HbT4_raw = squeeze(H2N_data (:,:,:,3,channels));     HbT4_raw = reshape (HbT4_raw,17,30*256,7);
%     
%     datasave(1).name = {'Neutral_pow_HbO','Happy_pow_HbO','N2H_pow_HbO','H2N_pow_HbO'};    
%     datasave(2).name = {'Neutral_pow_HbR','Happy_pow_HbR','N2H_pow_HbR','H2N_pow_HbR'};    
%     datasave(3).name = {'Neutral_pow_HbT','Happy_pow_HbT','N2H_pow_HbT','H2N_pow_HbT'}; 
%     
%     datasave(4).name = {'Neutral_amp_HbO','Happy_amp_HbO','N2H_amp_HbO','H2N_amp_HbO'};    
%     datasave(5).name = {'Neutral_amp_HbR','Happy_amp_HbR','N2H_amp_HbR','H2N_amp_HbR'};    
%     datasave(6).name = {'Neutral_amp_HbT','Happy_amp_HbT','N2H_amp_HbT','H2N_amp_HbT'};
% %----------------------------------------
%     samplerate = 50;                    
%     data_n = size(HbO1_raw,2);
%     data_N = 2^nextpow2(data_n); 
%     data_hz = linspace(0,samplerate/2,floor(data_N/2+1));
% 
% 
% for condi = 1:4;
%     dataX = eval(['fft([','HbO', num2str(condi),'_raw],data_N,2)']);
%     data_amp = abs(dataX)./data_N;
%     data_amp = data_amp(:,1:data_N/2+1,:);
%     data_amp (:,2:end-1,:) = 2.*data_amp(:,2:end-1,:);
%     eval ([datasave(4).name{condi}, '= data_amp;']);    
%     data_pow = data_amp.^2;
%     eval ([datasave(1).name{condi}, '= data_pow;']);
% end
% 
% for condi = 1:4;
%     dataX = eval(['fft([','HbR', num2str(condi),'_raw],data_N,2)']);
%     data_amp = abs(dataX)./data_N;
%     data_amp = data_amp(:,1:data_N/2+1,:);
%     data_amp (:,2:end-1,:) = 2.*data_amp(:,2:end-1,:);
%     eval ([datasave(5).name{condi}, '= data_amp;']);      
%     data_pow = data_amp.^2;
%     eval ([datasave(2).name{condi}, '= data_pow;']);
% end
% 
% for condi = 1:4;
%     dataX = eval(['fft([','HbT', num2str(condi),'_raw],data_N,2)']);
%     data_amp = abs(dataX)./data_N;
%     data_amp = data_amp(:,1:data_N/2+1,:);
%     data_amp (:,2:end-1,:) = 2.*data_amp(:,2:end-1,:);
%     eval ([datasave(6).name{condi}, '= data_amp;']);  
%     data_pow = data_amp.^2;
%     eval ([datasave(3).name{condi}, '= data_pow;']);
% end
% 
% save fft_sigle_trial.mat data_hz...
%      Neutral_pow_HbO Happy_pow_HbO N2H_pow_HbO H2N_pow_HbO Neutral_amp_HbO Happy_amp_HbO N2H_amp_HbO H2N_amp_HbO...
%      Neutral_pow_HbR Happy_pow_HbR N2H_pow_HbR H2N_pow_HbR Neutral_amp_HbR Happy_amp_HbR N2H_amp_HbR H2N_amp_HbR...
%      Neutral_pow_HbT Happy_pow_HbT N2H_pow_HbT H2N_pow_HbT Neutral_amp_HbT Happy_amp_HbT N2H_amp_HbT H2N_amp_HbT
% %% Plot power of each channels
% 
% close all
% for chi = 1:7;
% figure (chi), clf
% set (gcf,'color','w')
% plot (data_hz, squeeze (nanmean(Neutral_amp_HbO(:,:,chi))),'-k','linewidth',3);
% hold on
% plot (data_hz, squeeze (nanmean(Happy_amp_HbO(:,:,chi))),'-','color',[.6,.6,.6],'linewidth',3);
% hold on
% plot (data_hz, squeeze (nanmean(N2H_amp_HbO(:,:,chi))),'-b','linewidth',3);
% hold on
% plot (data_hz, squeeze (nanmean(H2N_amp_HbO(:,:,chi))),'-r','linewidth',3);
% set (gca,'xlim',[0 0.3],'xtick',0:0.1:0.5,'linewidth',3)
% set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
% title (['Ch' num2str(chi)],'FontSize',28,'fontweight','bold','fontname','arial black')
% % xlabel ('frequency (Hz)','FontSize',16,'fontweight','bold','fontname','arial black')
% % ylabel (['Amplitude(','\mu','V)'],'FontSize',16,'fontweight','bold','fontname','arial black')
% % legend Neutral Happy N2H H2N
% end
% 
% % close all
% % 
% % data2plot (:,1) = squeeze (nanmean(Neutral_amp_HbO(:,:,chi)))';
% % data2plot (:,2) = squeeze (nanmean(Happy_amp_HbO(:,:,chi)))';
% % data2plot (:,3) = squeeze (nanmean(N2H_amp_HbO(:,:,chi)))';
% % data2plot (:,4) = squeeze (nanmean(H2N_amp_HbO(:,:,chi)))';
% % 
% % for chi = 1:7;
% %     set (gcf,'color','w')
% %     plot3 (repmat(data_hz',1,4), repmat([1,2,3,4],2049,1),data2plot,'linewidth',2.5);
% %     set (gca,'xlim',[.1 1],'xtick',.1:.1:1,'ytick',1:3,'zlim',[0 2*1e-03])
% %     grid on
% %         
% % end
 
 