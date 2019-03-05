clear all
clc
%% Extract the whole frontal brain data
load fft_SNRZ

amp_all_2 = zeros (17,4);
amp_all_4 = zeros (17,4);
%--------------------------------0.2 Hz--------------
amp_all_2 (:,1) = (squeeze(nanmean (Neutral_Z_HbO(:,17,:),3))); 
amp_all_2 (:,2) = (squeeze(nanmean (Happy_Z_HbO(:,17,:),3)));
amp_all_2 (:,3) = (squeeze(nanmean (N2H_Z_HbO(:,17,:),3)));
amp_all_2 (:,4) = (squeeze(nanmean (H2N_Z_HbO(:,17,:),3)));

% xlswrite ('amp_2',amp_all_2); 

%--------------------------------0.4 Hz--------------
amp_all_4 (:,1) = (squeeze(nanmean (Neutral_Z_HbO(:,33,:),3)));
amp_all_4 (:,2) = (squeeze(nanmean (Happy_Z_HbO(:,33,:),3)));
amp_all_4 (:,3) = (squeeze(nanmean (N2H_Z_HbO(:,33,:),3)));
amp_all_4 (:,4) = (squeeze(nanmean (H2N_Z_HbO(:,33,:),3)));

% xlswrite ('amp_4',amp_all_4); 
%% Extract the left and right hemisphere data
load fft_SNRZ

amp_L_2 = zeros (17,4);
amp_R_2 = zeros (17,4);
%--------------------------------0.2 Hz--------------
amp_L_2 (:,1) = (squeeze(nanmean (Neutral_Z_HbO(:,17,[3 4 7]),3))); 
amp_L_2 (:,2) = (squeeze(nanmean (Happy_Z_HbO(:,17,[3 4 7]),3)));
amp_L_2 (:,3) = (squeeze(nanmean (N2H_Z_HbO(:,17,[3 4 7]),3)));
amp_L_2 (:,4) = (squeeze(nanmean (H2N_Z_HbO(:,17,[3 4 7]),3)));

amp_R_2 (:,1) = (squeeze(nanmean (Neutral_Z_HbO(:,17,[1 2 5]),3))); 
amp_R_2 (:,2) = (squeeze(nanmean (Happy_Z_HbO(:,17,[1 2 5]),3)));
amp_R_2 (:,3) = (squeeze(nanmean (N2H_Z_HbO(:,17,[1 2 5]),3)));
amp_R_2 (:,4) = (squeeze(nanmean (H2N_Z_HbO(:,17,[1 2 5]),3)));

%% extract info from each channel
clear all
clc

load fft_SNRZ

amp_each_2 = zeros (17,4,7);
% amp_each_4 = zeros (17,4,7);
%--------------------------------0.2 Hz--------------
amp_each_2 (:,1,:) = (squeeze(Neutral_Z_HbO(:,17,:))); 
amp_each_2 (:,2,:) = (squeeze(Happy_Z_HbO(:,17,:)));
amp_each_2 (:,3,:) = (squeeze(N2H_Z_HbO(:,17,:)));
amp_each_2 (:,4,:) = (squeeze(H2N_Z_HbO(:,17,:)));

% %--------------------------------0.4 Hz--------------
% amp_each_4 (:,1,:) = (squeeze(Neutral_Z_HbO(:,33,:)));
% amp_each_4 (:,2,:) = (squeeze(Happy_Z_HbO(:,33,:)));
% amp_each_4 (:,3,:) = (squeeze(N2H_Z_HbO(:,33,:)));
% amp_each_4 (:,4,:) = (squeeze(H2N_Z_HbO(:,33,:)));


%% Extract amplitude of raw data of 0.2 and 0.4 Hz
load data(.08-.5)_fft.mat
%-----------------------------------------------------------------HbO
amp_HbO_2 = zeros (17,7,4);
amp_HbO_4 = zeros (17,7,4);
%--------------------------------0.2 Hz--------------
amp_HbO_2 (:,:,1) = squeeze(Neutral_amp_HbO(:,18,:));
amp_HbO_2 (:,:,2) = squeeze(Happy_amp_HbO(:,18,:));
amp_HbO_2 (:,:,3) = squeeze(N2H_amp_HbO(:,18,:));
amp_HbO_2 (:,:,4) = squeeze(H2N_amp_HbO(:,18,:));
%--------------------------------0.4 Hz--------------
amp_HbO_4 (:,:,1) = squeeze(Neutral_amp_HbO(:,34,:));
amp_HbO_4 (:,:,2) = squeeze(Happy_amp_HbO(:,34,:));
amp_HbO_4 (:,:,3) = squeeze(N2H_amp_HbO(:,34,:));
amp_HbO_4 (:,:,4) = squeeze(H2N_amp_HbO(:,34,:));

%-----------------------------------------------------------------HbR
amp_HbR_2 = zeros (17,7,4);
amp_HbR_4 = zeros (17,7,4);
%--------------------------------0.2 Hz--------------
amp_HbR_2 (:,:,1) = squeeze(Neutral_amp_HbR(:,18,:));
amp_HbR_2 (:,:,2) = squeeze(Happy_amp_HbR(:,18,:));
amp_HbR_2 (:,:,3) = squeeze(N2H_amp_HbR(:,18,:));
amp_HbR_2 (:,:,4) = squeeze(H2N_amp_HbR(:,18,:));
%--------------------------------0.4 Hz--------------
amp_HbR_4 (:,:,1) = squeeze(Neutral_amp_HbR(:,34,:));
amp_HbR_4 (:,:,2) = squeeze(Happy_amp_HbR(:,34,:));
amp_HbR_4 (:,:,3) = squeeze(N2H_amp_HbR(:,34,:));
amp_HbR_4 (:,:,4) = squeeze(H2N_amp_HbR(:,34,:));

%-----------------------------------------------------------------HbT
amp_HbT_2 = zeros (17,7,4);
amp_HbT_4 = zeros (17,7,4);
%--------------------------------0.2 Hz--------------
amp_HbT_2 (:,:,1) = squeeze(Neutral_amp_HbT(:,18,:));
amp_HbT_2 (:,:,2) = squeeze(Happy_amp_HbT(:,18,:));
amp_HbT_2 (:,:,3) = squeeze(N2H_amp_HbT(:,18,:));
amp_HbT_2 (:,:,4) = squeeze(H2N_amp_HbT(:,18,:));
%--------------------------------0.4 Hz--------------
amp_HbT_4 (:,:,1) = squeeze(Neutral_amp_HbT(:,34,:));
amp_HbT_4 (:,:,2) = squeeze(Happy_amp_HbT(:,34,:));
amp_HbT_4 (:,:,3) = squeeze(N2H_amp_HbT(:,34,:));
amp_HbT_4 (:,:,4) = squeeze(H2N_amp_HbT(:,34,:));

save amp_raw amp_HbO_2 amp_HbO_4  amp_HbR_2 amp_HbR_4 amp_HbT_2 amp_HbT_4

%% Extract amplitude and power of deconvolution data 0.2 and 0.4 Hz
load data_deconv_fft
%-----------------------------------------------------------------HbO
amp_HbO_deconv_2 = zeros (17,7,4);
amp_HbO_deconv_4 = zeros (17,7,4);
%--------------------------------0.2 Hz--------------
amp_HbO_deconv_2 (:,:,1) = squeeze(Neutral_amp_HbO_deconv(:,18,:));
amp_HbO_deconv_2 (:,:,2) = squeeze(Happy_amp_HbO_deconv(:,18,:));
amp_HbO_deconv_2 (:,:,3) = squeeze(N2H_amp_HbO_deconv(:,18,:));
amp_HbO_deconv_2 (:,:,4) = squeeze(H2N_amp_HbO_deconv(:,18,:));
%--------------------------------0.4 Hz--------------
amp_HbO_deconv_4 (:,:,1) = squeeze(Neutral_amp_HbO_deconv(:,34,:));
amp_HbO_deconv_4 (:,:,2) = squeeze(Happy_amp_HbO_deconv(:,34,:));
amp_HbO_deconv_4 (:,:,3) = squeeze(N2H_amp_HbO_deconv(:,34,:));
amp_HbO_deconv_4 (:,:,4) = squeeze(H2N_amp_HbO_deconv(:,34,:));

%-----------------------------------------------------------------HbR
amp_HbR_deconv_2 = zeros (17,7,4);
amp_HbR_deconv_4 = zeros (17,7,4);
%--------------------------------0.2 Hz--------------
amp_HbR_deconv_2 (:,:,1) = squeeze(Neutral_amp_HbR_deconv(:,18,:));
amp_HbR_deconv_2 (:,:,2) = squeeze(Happy_amp_HbR_deconv(:,18,:));
amp_HbR_deconv_2 (:,:,3) = squeeze(N2H_amp_HbR_deconv(:,18,:));
amp_HbR_deconv_2 (:,:,4) = squeeze(H2N_amp_HbR_deconv(:,18,:));
%--------------------------------0.4 Hz--------------
amp_HbR_deconv_4 (:,:,1) = squeeze(Neutral_amp_HbR_deconv(:,34,:));
amp_HbR_deconv_4 (:,:,2) = squeeze(Happy_amp_HbR_deconv(:,34,:));
amp_HbR_deconv_4 (:,:,3) = squeeze(N2H_amp_HbR_deconv(:,34,:));
amp_HbR_deconv_4 (:,:,4) = squeeze(H2N_amp_HbR_deconv(:,34,:));

%-----------------------------------------------------------------HbT
amp_HbT_deconv_2 = zeros (17,7,4);
amp_HbT_deconv_4 = zeros (17,7,4);
%--------------------------------0.2 Hz--------------
amp_HbT_deconv_2 (:,:,1) = squeeze(Neutral_amp_HbT_deconv(:,18,:));
amp_HbT_deconv_2 (:,:,2) = squeeze(Happy_amp_HbT_deconv(:,18,:));
amp_HbT_deconv_2 (:,:,3) = squeeze(N2H_amp_HbT_deconv(:,18,:));
amp_HbT_deconv_2 (:,:,4) = squeeze(H2N_amp_HbT_deconv(:,18,:));
%--------------------------------0.4 Hz--------------
amp_HbT_deconv_4 (:,:,1) = squeeze(Neutral_amp_HbT_deconv(:,34,:));
amp_HbT_deconv_4 (:,:,2) = squeeze(Happy_amp_HbT_deconv(:,34,:));
amp_HbT_deconv_4 (:,:,3) = squeeze(N2H_amp_HbT_deconv(:,34,:));
amp_HbT_deconv_4 (:,:,4) = squeeze(H2N_amp_HbT_deconv(:,34,:));

save amp_deconv amp_HbO_deconv_2 amp_HbO_deconv_4  amp_HbR_deconv_2 amp_HbR_deconv_4 amp_HbT_deconv_2 amp_HbT_deconv_4


