%% Preparation
clear
close all

%%
addpath(genpath('../../dependencies/ECGdeli'));
addpath(genpath('../../algorithms/'));
input=load('sim_signals.mat');

% Define some important parameters
samplerate=1000;
nfeatures=18;


% Define all sweeping parameters
noiseLvl=[10,20,30]; %dB
lowPassFg=[20,40,50,60,70,80,100,150,200];%Hz
highPassFg=[0.05,0.1,0.2,0.3,0.4,0.5];%Hz
combinationsFg=combvec(lowPassFg,highPassFg)';


for dBlvl=1:1:size(noiseLvl)

% Study variables
noisedB=noiseLvl(dBlvl);
nRep=10;
numSetups=length(input.ECGs);
ecgNoise=cell(numSetups,1);
ecgMaxQRS=cell(numSetups,1);
ecg12Ch=cell(numSetups,1);
featureMatrix_clean=cell(numSetups,1);
featureMatrix_cleanFilt=cell(numSetups,size(combinationsFg,1));
featureMatrix_noise=cell(numSetups,size(combinationsFg,1),size(input.ECGs(1).signal,2));
FPT=cell(numSetups,1);


% Evaluate all ECGs
parfor setupsNo=1:length(input.ECGs)
    %% Clean experiments
    disp(['Starting with Signal ' num2str(setupsNo)])
    ecgMatrixUnique=input.ECGs(setupsNo).signal;
    samplerate=input.ECGs(setupsNo).fs;
    
    %Calculate features
    [~,FPT{setupsNo,1}]=getFPTFromSimulations(ecgMatrixUnique,samplerate);
    
    [featureMatrix_clean{setupsNo,1},featureNames]=calculateFeaturesOneBeat(ecgMatrixUnique,samplerate,1,FPT{setupsNo,1});
    
    %% Noisy
    featureMatrix_noise_tmp=cell(1,size(combinationsFg,1),size(input.ECGs(setupsNo).signal,2));
    featureMatrix_cleanFilt_tmp=cell(1,size(combinationsFg,1),size(input.ECGs(setupsNo).signal,2));
    for fg=1:size(combinationsFg,1)
        fgHigh=combinationsFg(fg,2);
        fgLow=combinationsFg(fg,1);
        
        % Evaluate only filtering without noise
        ecg_cleanFilt=ECG_High_Low_Filter(ecgMatrixUnique,samplerate,fgHigh,fgLow);
        [~,FPTcellfilt]=getFPTFromSimulations(ecgMatrixUnique,samplerate);
        featureMatrix_cleanFilt_tmp{1,fg}=calculateFeaturesOneBeat(ecg_cleanFilt,samplerate,1,FPTcellfilt);
        
        % Evaluate only filtering with noise
        featureMatrix_noise_tmptmp=cell(1,1,size(input.ECGs(setupsNo).signal,2));
        for ld=1:size(input.ECGs(setupsNo).signal,2)
            featRep=nan(nfeatures,nRep);
            for rep=1:nRep
                tmp=repmat(ecgMatrixUnique(:,ld),3,1)+wgn(3*size(ecgMatrixUnique(:,ld),1),1,var(repmat(ecgMatrixUnique(:,ld),3,1))/10^(noisedB/10),'linear');
                tmp=ECG_High_Low_Filter(tmp,samplerate,fgHigh,fgLow);
                ecgNoise=tmp(size(ecgMatrixUnique,1)+1:size(ecgMatrixUnique,1)*2);
                %Calculate features
                featRep(:,rep)=calculateFeaturesOneBeat(ecgNoise,samplerate,1,FPTcellfilt{ld,1});
            end
            featureMatrix_noise_tmptmp{1,1,ld}=featRep;
        end
        featureMatrix_noise_tmp(1,fg,:)=featureMatrix_noise_tmptmp;
    end
    featureMatrix_cleanFilt(setupsNo,:)=featureMatrix_cleanFilt_tmp;
    featureMatrix_noise(setupsNo,:,:)=featureMatrix_noise_tmp;
end

save(strcat('Results_',num2str(noisedB),'dB'),'featureMatrix_noise','featureMatrix_clean','featureMatrix_cleanFilt','combinationsFg','nfeatures');
end


% %% Visualize the results
% 
% %% Calculate errors
% featAll_clean=cell(size(featureMatrix_noise,3),1);
% featAll_cleanFilt=cell(size(featureMatrix_noise,3),1);
% featAll_noise=cell(size(featureMatrix_noise,3),1);
% err_noise=cell(size(featureMatrix_noise,3),1);
% err_cleanFilt=cell(size(featureMatrix_noise,3),1);
% relErr_noise=cell(size(featureMatrix_noise,3),1);
% relErr_cleanFilt=cell(size(featureMatrix_noise,3),1);
% 
% numComb=size(combinationsFg,1);
% for ld=1:1:size(featureMatrix_noise,3)
%     featAll_clean{ld,1}=cell2mat(cellfun(@(x) x(:,ld),featureMatrix_clean(:,:),'UniformOutput',0));
%     featAll_cleanFilt{ld,1}=cell2mat(cellfun(@(x) x(:,ld),featureMatrix_cleanFilt(:,:),'UniformOutput',0));
%     featAll_noise{ld,1}=cell2mat(cellfun(@(x)reshape(x,numel(x),1),featureMatrix_noise(:,:,ld),'UniFormOutput',0));
%     
%     err_noise{ld,1}=featAll_noise{ld,1}-repmat(featAll_clean{ld,1},nRep,numComb);
%     relErr_noise{ld,1}=err_noise{ld,1}./repmat(featAll_clean{ld,1},nRep,numComb);
%     
%     err_cleanFilt{ld,1}=featAll_cleanFilt{ld,1}-repmat(featAll_clean{ld,1},1,numComb);
%     relErr_cleanFilt{ld,1}=err_cleanFilt{ld,1}./repmat(featAll_clean{ld,1},1,numComb);
% end
% %save(['./FeatureStudy_Effects_' num2str(noiseLvl(dBlvl)) '.mat'])
% 
% 
% % Get the error for different filtering
% meanErrNoise=nan(size(featureMatrix_noise,3),nfeatures);
% stdErrNois=nan(size(featureMatrix_noise,3),nfeatures);
% for ld=1:1:8 % leads
%     for feat=1:1:18 % features
%         meanErrNoise(ld,feat)=mean(relErr_noise{ld,1}(feat:nfeatures:end))*100; %relative error in %
%         stdErrNois(ld,feat)=std(relErr_noise{ld,1}(feat:nfeatures:end))*100; %relative error in %
%     end
% end
% 
% 
% meanErrFilt=nan(size(featureMatrix_noise,3),nfeatures);
% stdErrFilt=nan(size(featureMatrix_noise,3),nfeatures);
% for ld=1:1:8 % leads
%     for feat=1:1:18 % features
%         meanErrFilt(ld,feat)=mean(relErr_noise{ld,1}(feat:nfeatures:end))*100; %relative error in %
%         stdErrFilt(ld,feat)=std(relErr_noise{ld,1}(feat:nfeatures:end))*100; %relative error in %
%     end
% end



