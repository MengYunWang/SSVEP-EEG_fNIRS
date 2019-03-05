%% Pre-processing
% The preprocessing consisted of two stages:
% Stage I: preICA and run ICA
%   preICA: Import data> Channel Location > Filter > resample the data>  
%           rereference the data > extract all epoches > reject bad trials
%   run ICA: runica > exclude artifacts with ADjust
%   Baseline correction

% Stage II: extract conditions

% If you do not konw which function you should run one of your data with EEGlab
% and then type EEG.history.
% created by M.-Y. Wang
% 02-06-2019

%% Stage1-preICA
clear all
clc
cd ('F:\Fusion\Face-sequence\EEG\Raw')

for subi = 105:121;
    EEG = pop_biosig([num2str(subi),'-sequence.bdf'], 'channels',1:64,'ref',47);
    EEG.setname = [num2str(subi),'_preICA'];
    EEG = pop_chanedit(EEG, 'lookup','D:\\Program Files\\MATLAB\\R2014a\\matlabtoolbox\\eeglab14_1_2b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');
    EEG =pop_eegfiltnew(EEG, 0.5,45,13518,0,[],0);
%     EEG = pop_eegfiltnew(EEG,1,30,6760,0,[],0);
%     EEG = pop_eegfiltnew(EEG, [],0.5,13518,1,[],0);
%     EEG = pop_eegfiltnew(EEG, [],1,6760,1,[],0);
    EEG = pop_resample( EEG, 250);
    EEG = pop_reref( EEG, []);
    EEG = pop_epoch( EEG, {  '1'  '2' '3' '4'}, [-0.7 2.3], 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-200 0]);
    EEG.setname = [num2str(subi),'_preICA'];
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), '_preICA.set'],'filepath','F:\Fusion\Face-sequence\EEG\Pre_processing\PreICA');
end
%% disgard bad trials



%% Stage1-RunICA
clear all
clc
cd ('F:\Fusion\face-sequence\EEG\Pre_processing\PreICA')

for subi=105:121;
    EEG = pop_loadset('filename',[num2str(subi),'_PreICA.set'],'filepath','F:\Fusion\face-sequence\EEG\Pre_processing\PreICA\');
    EEG = pop_reref( EEG, []);
    EEG = pop_runica(EEG,'icatype','fastica','approach','symm','numofIC',35);
     EEG.setname = [num2str(subi),'_runICA'];
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), '_runICA.set'],'filepath','F:\Fusion\face-sequence\EEG\Pre_processing\ICA\runica');
end

%% Stage1-ICA exclude
% clear all
% clc
% cd ('F:\Fusion\Face-sequence\EEG\\Pre_processing\ICA')
% 
% for subi=105:121;
%     EEG = pop_loadset('filename',[num2str(subi),'_runICA.set'],'filepath','F:\Fusion\Face-sequence\EEG\\Pre_processing\ICA');
%     EEG = interface_ADJ (EEG, [num2str(subi),'.txt']);
% end

%% run it again
clear all
clc

for subi= 105:121;
    EEG = pop_loadset('filename',[num2str(subi),'_pruned with ICA.set'],'filepath','F:\Fusion\face-sequence\EEG\Pre_processing\ICA\pruned');
    EEG = pop_reref( EEG, []);
    EEG = pop_rmbase( EEG, [-200 0]);
    EEG = pop_runica(EEG,'icatype','fastica','approach','symm','numofIC',10);
    EEG.setname = [num2str(subi),'_pruned with ICA'];
    EEG = pop_saveset( EEG, 'filename',[num2str(subi),'_pruned with ICA.set'],'filepath','F:\Fusion\face-sequence\EEG\Pre_processing\ICA\pruned');
end
%% Stage2-select the conditions
clear all
clc
cd ('F:\Fusion\Face-sequence\EEG\\Pre_processing\ICA\pruned')

data(1).name = '_Neutral.set';data(2).name = '_Happy.set';data(3).name = '_N2H.set';data(4).name = '_H2N.set';
data(1).file = 'Condition1_Neutral'; data(2).file = 'Condition2_Happy';data(3).file = 'Condition3_N2H';data(4).file = 'Condition4_H2N';
for subi = 105:121;
    EEG = pop_loadset('filename',[num2str(subi),'_pruned with ICA.set'],'filepath','F:\Fusion\face-sequence\EEG\Pre_processing\ICA\pruned');
    EEG = pop_selectevent( EEG, 'type',1,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(1).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(1).name],'filepath',['F:\Fusion\face-sequence\EEG\Pre_processing\Conditions\',data(1).file]);
end

for subi = 105:121;
    EEG = pop_loadset('filename',[num2str(subi),'_pruned with ICA.set'],'filepath','F:\Fusion\Face-sequence\EEG\Pre_processing\ICA\pruned');
    EEG = pop_selectevent( EEG, 'type',2,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(2).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(2).name],'filepath',['F:\Fusion\Face-sequence\EEG\\Pre_processing\Conditions\',data(2).file]);
end

for subi = 105:121;
    EEG = pop_loadset('filename',[num2str(subi),'_pruned with ICA.set'],'filepath','F:\Fusion\Face-sequence\EEG\\Pre_processing\ICA\pruned');
    EEG = pop_selectevent( EEG, 'type',3,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(3).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(3).name],'filepath',['F:\Fusion\Face-sequence\EEG\\Pre_processing\Conditions\',data(3).file]);
end

for subi = 105:121;
    EEG = pop_loadset('filename',[num2str(subi),'_pruned with ICA.set'],'filepath','F:\Fusion\Face-sequence\EEG\\Pre_processing\ICA\pruned');
    EEG = pop_selectevent( EEG, 'type',4,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(4).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(4).name],'filepath',['F:\Fusion\Face-sequence\EEG\\Pre_processing\Conditions\',data(4).file]);
end
