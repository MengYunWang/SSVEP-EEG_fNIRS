clear all
clc
%% Extract amplitude and power 

amp_all_10 = zeros (17,8);
amp_all_20 = zeros (17,8);
%--------------------------------10 Hz--------------
amp_all_10 (:,1) = (squeeze(mean (Neutral_amp([26:30,63,64],20,:))));
amp_all_10 (:,2) = (squeeze(mean (Happy_amp([26:30,63,64],20,:))));
amp_all_10 (:,3) = (squeeze(mean (N2H_amp([26:30,63,64],20,:))));
amp_all_10 (:,4) = (squeeze(mean (H2N_amp([26:30,63,64],20,:))));
amp_all_10 (:,5) = (squeeze(mean (Neutral_amp([1 3 33 34 36 37],20,:))));
amp_all_10 (:,6) = (squeeze(mean (Happy_amp([1 3 33 34 36 37],20,:))));
amp_all_10 (:,7) = (squeeze(mean (N2H_amp([1 3 33 34 36 37],20,:))));
amp_all_10 (:,8) = (squeeze(mean (H2N_amp([1 3 33 34 36 37],20,:))));
xlswrite ('amp_10',amp_all_10); 

%--------------------------------20 Hz--------------
amp_all_20 (:,1) = (squeeze(mean (Neutral_amp([26:30,63,64],39,:))));
amp_all_20 (:,2) = (squeeze(mean (Happy_amp([26:30,63,64],39,:))));
amp_all_20 (:,3) = (squeeze(mean (N2H_amp([26:30,63,64],39,:))));
amp_all_20 (:,4) = (squeeze(mean (H2N_amp([26:30,63,64],39,:))));
amp_all_20 (:,5) = (squeeze(mean (Neutral_amp([1 3 33 34 36 37],39,:))));
amp_all_20 (:,6) = (squeeze(mean (Happy_amp([1 3 33 34 36 37],39,:))));
amp_all_20 (:,7) = (squeeze(mean (N2H_amp([1 3 33 34 36 37],39,:))));
amp_all_20 (:,8) = (squeeze(mean (H2N_amp([1 3 33 34 36 37],39,:))));
xlswrite ('amp_20',amp_all_20); 


%%
load fft_SNRZ

amp_all_10 = zeros (17,8);
amp_all_20 = zeros (17,8);
%--------------------------------10 Hz--------------
amp_all_10 (:,1) = (squeeze(mean (Neutral_Z([26:30,63,64],20,:)))); % occipital
amp_all_10 (:,2) = (squeeze(mean (Happy_Z([26:30,63,64],20,:))));
amp_all_10 (:,3) = (squeeze(mean (N2H_Z([26:30,63,64],20,:))));
amp_all_10 (:,4) = (squeeze(mean (H2N_Z([26:30,63,64],20,:))));
amp_all_10 (:,5) = (squeeze(mean (Neutral_Z([1 3 33 34 36 37],20,:))));
amp_all_10 (:,6) = (squeeze(mean (Happy_Z([1 3 33 34 36 37],20,:))));
amp_all_10 (:,7) = (squeeze(mean (N2H_Z([1 3 33 34 36 37],20,:))));
amp_all_10 (:,8) = (squeeze(mean (H2N_Z([1 3 33 34 36 37],20,:))));
xlswrite ('amp_10',amp_all_10); 

%--------------------------------20 Hz--------------
amp_all_20 (:,1) = (squeeze(mean (Neutral_Z([26:30,63,64],39,:))));
amp_all_20 (:,2) = (squeeze(mean (Happy_Z([26:30,63,64],39,:))));
amp_all_20 (:,3) = (squeeze(mean (N2H_Z([26:30,63,64],39,:))));
amp_all_20 (:,4) = (squeeze(mean (H2N_Z([26:30,63,64],39,:))));
amp_all_20 (:,5) = (squeeze(mean (Neutral_Z([1 3 33 34 36 37],39,:))));
amp_all_20 (:,6) = (squeeze(mean (Happy_Z([1 3 33 34 36 37],39,:))));
amp_all_20 (:,7) = (squeeze(mean (N2H_Z([1 3 33 34 36 37],39,:))));
amp_all_20 (:,8) = (squeeze(mean (H2N_Z([1 3 33 34 36 37],39,:))));
xlswrite ('amp_20',amp_all_20); 

save ampZ_10+20 amp_all_10 amp_all_20


%%
load fft_SNRZ

amp_L_10 = zeros (17,8); amp_R_10 = zeros (17,8);
amp_L_20 = zeros (17,8); amp_R_20 = zeros (17,8);
%--------------------------------10 Hz--------------
amp_L_10 (:,1) = (squeeze(mean (Neutral_Z(25:27,20,:)))); % occipital
amp_L_10 (:,2) = (squeeze(mean (Happy_Z(25:27,20,:))));
amp_L_10 (:,3) = (squeeze(mean (N2H_Z(25:27,20,:))));
amp_L_10 (:,4) = (squeeze(mean (H2N_Z(25:27,20,:))));
amp_L_10 (:,5) = (squeeze(mean (Neutral_Z(1:3,20,:))));
amp_L_10 (:,6) = (squeeze(mean (Happy_Z(1:3,20,:))));
amp_L_10 (:,7) = (squeeze(mean (N2H_Z(1:3,20,:))));
amp_L_10 (:,8) = (squeeze(mean (H2N_Z(1:3,20,:))));


amp_R_10 (:,1) = (squeeze(mean (Neutral_Z(62:64,20,:)))); % occipital
amp_R_10 (:,2) = (squeeze(mean (Happy_Z(62:64,20,:))));
amp_R_10 (:,3) = (squeeze(mean (N2H_Z(62:64,20,:))));
amp_R_10 (:,4) = (squeeze(mean (H2N_Z(62:64,20,:))));
amp_R_10 (:,5) = (squeeze(mean (Neutral_Z(34:36,20,:))));
amp_R_10 (:,6) = (squeeze(mean (Happy_Z(34:36,20,:))));
amp_R_10 (:,7) = (squeeze(mean (N2H_Z(34:36,20,:))));
amp_R_10 (:,8) = (squeeze(mean (H2N_Z(34:36,20,:))));

%--------------------------------20 Hz--------------
amp_L_20 (:,1) = (squeeze(mean (Neutral_Z(25:27,39,:)))); % occipital
amp_L_20 (:,2) = (squeeze(mean (Happy_Z(25:27,39,:))));
amp_L_20 (:,3) = (squeeze(mean (N2H_Z(25:27,39,:))));
amp_L_20 (:,4) = (squeeze(mean (H2N_Z(25:27,39,:))));
amp_L_20 (:,5) = (squeeze(mean (Neutral_Z(1:3,39,:))));
amp_L_20 (:,6) = (squeeze(mean (Happy_Z(1:3,39,:))));
amp_L_20 (:,7) = (squeeze(mean (N2H_Z(1:3,39,:))));
amp_L_20 (:,8) = (squeeze(mean (H2N_Z(1:3,39,:))));


amp_R_20 (:,1) = (squeeze(mean (Neutral_Z(62:64,39,:)))); % occipital
amp_R_20 (:,2) = (squeeze(mean (Happy_Z(62:64,39,:))));
amp_R_20 (:,3) = (squeeze(mean (N2H_Z(62:64,39,:))));
amp_R_20 (:,4) = (squeeze(mean (H2N_Z(62:64,39,:))));
amp_R_20 (:,5) = (squeeze(mean (Neutral_Z(34:36,39,:))));
amp_R_20 (:,6) = (squeeze(mean (Happy_Z(34:36,39,:))));
amp_R_20 (:,7) = (squeeze(mean (N2H_Z(34:36,39,:))));
amp_R_20 (:,8) = (squeeze(mean (H2N_Z(34:36,39,:))));
