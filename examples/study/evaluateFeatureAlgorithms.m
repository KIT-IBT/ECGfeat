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


for dBlvl=1:1:size(noiseLvl,2)
    
    % Study variables
    noisedB=noiseLvl(dBlvl);
    nRep=50;
    numSetups=length(input.ECGs);
    ecgNoise=cell(numSetups,1);
    ecgMaxQRS=cell(numSetups,1);
    ecg12Ch=cell(numSetups,1);
    featureMatrix_clean=cell(numSetups,1);
    featureMatrix_cleanFilt=cell(numSetups,size(combinationsFg,1));
    featureMatrix_noise=cell(numSetups,size(combinationsFg,1),size(input.ECGs(1).signal,2));
    FPT=cell(numSetups,1);
    
    if ~exist(strcat('Results_',num2str(noisedB),'dB.mat'),'file')
        
        % Evaluate all ECGs
        for setupsNo=1:length(input.ECGs)
            %% Clean experiments
            disp(['Starting with Signal ' num2str(setupsNo)])
            ecgMatrixUnique=input.ECGs(setupsNo).signal;
            samplerate=input.ECGs(setupsNo).fs;
            
            %Calculate features
            [~,FPT{setupsNo,1}]=getFPTFromSimulations(ecgMatrixUnique,samplerate);
            
            [featureMatrix_clean{setupsNo,1},featureNames]=calculateFeaturesOneBeat(ecgMatrixUnique,samplerate,0,FPT{setupsNo,1});
            
            %% Noisy
            featureMatrix_noise_tmp=cell(1,size(combinationsFg,1),size(input.ECGs(setupsNo).signal,2));
            featureMatrix_cleanFilt_tmp=cell(1,size(combinationsFg,1));
            for fg=1:size(combinationsFg,1)
                fgHigh=combinationsFg(fg,2);
                fgLow=combinationsFg(fg,1);
                
                % Evaluate only filtering without noise
                ecg_cleanFilt=ECG_High_Low_Filter(ecgMatrixUnique,samplerate,fgHigh,fgLow);
                [~,FPTcellfilt]=getFPTFromSimulations(ecgMatrixUnique,samplerate);
                featureMatrix_cleanFilt_tmp{1,fg}=calculateFeaturesOneBeat(ecg_cleanFilt,samplerate,0,FPTcellfilt);
                
                % Evaluate only filtering with noise
                featureMatrix_noise_tmptmp=cell(1,1,size(input.ECGs(setupsNo).signal,2));
                for ld=1:size(input.ECGs(setupsNo).signal,2)
                    featRep=nan(nfeatures,nRep);
                    for rep=1:nRep
                        tmp=repmat(ecgMatrixUnique(:,ld),3,1)+wgn(3*size(ecgMatrixUnique(:,ld),1),1,var(repmat(ecgMatrixUnique(:,ld),3,1))/10^(noisedB/10),'linear');
                        tmp=ECG_High_Low_Filter(tmp,samplerate,fgHigh,fgLow);
                        ecgNoise=tmp(size(ecgMatrixUnique,1)+1:size(ecgMatrixUnique,1)*2);
                        %Calculate features
                        featRep(:,rep)=calculateFeaturesOneBeat(ecgNoise,samplerate,0,FPTcellfilt{ld,1});
                    end
                    featureMatrix_noise_tmptmp{1,1,ld}=featRep;
                end
                featureMatrix_noise_tmp(1,fg,:)=featureMatrix_noise_tmptmp;
            end
            featureMatrix_cleanFilt(setupsNo,:)=featureMatrix_cleanFilt_tmp;
            featureMatrix_noise(setupsNo,:,:)=featureMatrix_noise_tmp;
        end
        
        save(strcat('Results_',num2str(noisedB),'dB.mat'),'featureMatrix_noise','featureMatrix_clean','featureMatrix_cleanFilt','combinationsFg','nfeatures','nRep');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize the results
%%%%%%%%%%%%%%%%%%%%%%%%%

for dBlvl=1:1:size(noiseLvl,2)
    noisedB=noiseLvl(dBlvl);
    load(strcat('Results_',num2str(noisedB),'dB.mat'))
    disp(['Evaluating results with SNR of ',num2str(noisedB),'dB'])
    
    % Calculate errors and relative errors
    featAll_clean=cell(size(featureMatrix_noise,3),1);
    featAll_cleanFilt=cell(size(featureMatrix_noise,3),1);
    featAll_noise=cell(size(featureMatrix_noise,3),1);
    featAll_noise_mean=cell(size(featureMatrix_noise,3),1);
    err_noise_mean=cell(size(featureMatrix_noise,3),1);
    err_noise=cell(size(featureMatrix_noise,3),1);
    err_cleanFilt=cell(size(featureMatrix_noise,3),1);
    relErr_noise=cell(size(featureMatrix_noise,3),1);
    relErr_cleanFilt=cell(size(featureMatrix_noise,3),1);
    relErr_noise_mean=cell(size(featureMatrix_noise,3),1);
    
    numComb=size(combinationsFg,1);
    for ld=1:1:size(featureMatrix_noise,3)
        featAll_clean{ld,1}=cell2mat(cellfun(@(x) x(:,ld),featureMatrix_clean(:,:),'UniformOutput',0));
        featAll_cleanFilt{ld,1}=cell2mat(cellfun(@(x) x(:,ld),featureMatrix_cleanFilt(:,:),'UniformOutput',0));
        featAll_noise{ld,1}=cell2mat(cellfun(@(x)reshape(x,numel(x),1),featureMatrix_noise(:,:,ld),'UniFormOutput',0));
        featAll_noise_mean{ld,1}=cell2mat(cellfun(@(x)mean(x,2),featureMatrix_noise(:,:,ld),'UniFormOutput',0));
        
        err_noise{ld,1}=featAll_noise{ld,1}-repmat(featAll_clean{ld,1},nRep,numComb);
        relErr_noise{ld,1}=err_noise{ld,1}./repmat(featAll_clean{ld,1},nRep,numComb);
        
        err_noise_mean{ld,1}=featAll_noise_mean{ld,1}-repmat(featAll_clean{ld,1},1,numComb);
        relErr_noise_mean{ld,1}=err_noise_mean{ld,1}./repmat(featAll_clean{ld,1},1,numComb);
        
        err_cleanFilt{ld,1}=featAll_cleanFilt{ld,1}-repmat(featAll_clean{ld,1},1,numComb);
        relErr_cleanFilt{ld,1}=err_cleanFilt{ld,1}./repmat(featAll_clean{ld,1},1,numComb);
    end
    
    
    %% Build the table for filtering recommendations
    % get the lowest lowpass frequency without disturbing the features
    idxLP=find(combinationsFg(:,2)==min(combinationsFg(:,2)));
    meanErrFilt=nan(size(featureMatrix_noise,3),nfeatures,length(idxLP));
    for ld=1:1:8 % leads
        for lp=1:1:length(idxLP)
            for feat=1:1:18 % features
                if feat==17
                    meanErrFilt(ld,feat,lp)=sum(abs(sign(err_cleanFilt{ld,1}(feat:nfeatures:end,lp)-1)))/sum(~isnan(err_cleanFilt{ld,1}(feat:nfeatures:end,lp)));
                else
                    meanErrFilt(ld,feat,lp)=nanmean(relErr_cleanFilt{ld,1}(feat:nfeatures:end,lp))*100; %relative error in %
                end
            end
        end
    end
    % get the frequency showing maximum 5% error
    tmp=squeeze(max(abs(meanErrFilt),[],1));
    tmp(17,:)=squeeze(min(abs(meanErrFilt(:,17,:)),[],1));
    bsm1prc=nan(nfeatures,1);
    bsm1prcValLP=nan(nfeatures,1);
    for ft=1:1:nfeatures
        if ft==17
            bsm1prc(ft,1)=find(tmp(ft,:)==max(tmp(ft,:)),1,'first');
            bsm1prcValLP(ft,1)=tmp(ft,bsm1prc(ft,1));
        else
            if ~isempty(find(tmp(ft,:)<=5,1,'first'))
                bsm1prc(ft,1)=find(tmp(ft,:)<=5,1,'first');
                bsm1prcValLP(ft,1)=tmp(ft,bsm1prc(ft,1));
            end
        end
    end
    coffFreqLP=combinationsFg(bsm1prc,1);
    
    
    % get the highest highpass frequency without disturbing the features
    idxLP=find(combinationsFg(:,1)==max(combinationsFg(:,1)));
    meanErrFilt=nan(size(featureMatrix_noise,3),nfeatures,length(idxLP));
    for ld=1:1:8 % leads
        for lp=1:1:length(idxLP)
            for feat=1:1:18 % features
                if feat==17
                    meanErrFilt(ld,feat,lp)=sum(abs(sign(err_cleanFilt{ld,1}(feat:nfeatures:end,lp)-1)))/sum(~isnan(err_cleanFilt{ld,1}(feat:nfeatures:end,lp)));
                else
                    meanErrFilt(ld,feat,lp)=nanmean(relErr_cleanFilt{ld,1}(feat:nfeatures:end,lp))*100; %relative error in %
                end
            end
        end
    end
    
    % get the frequency showing maximum 5% error
    tmp=squeeze(max(abs(meanErrFilt),[],1));
    bsm1prc=nan(nfeatures,1);
    bsm1prcValHP=nan(nfeatures,1);
    for ft=1:1:nfeatures
        if ft==17
            bsm1prc(ft,1)=find(tmp(ft,:)==max(tmp(ft,:)),1,'first');
            bsm1prcValHP(ft,1)=tmp(ft,bsm1prc(ft,1));
        else
            if ~isempty(find(tmp(ft,:)<=5,1,'first'))
                bsm1prc(ft,1)=find(tmp(ft,:)<=5,1,'first');
                bsm1prcValHP(ft,1)=tmp(ft,bsm1prc(ft,1));
            end
        end
    end
    coffFreqHP=combinationsFg(idxLP(bsm1prc),2);
    
    
    table1=cell(5,1);
    for i=1:1:length(bsm1prc)
        table1{1,1}=[table1{1,1},'Feat',num2str(i,'%02i'),char(9)];
        table1{2,1}=[table1{2,1},num2str(coffFreqLP(i,1),'%02i'),char(9)];
        table1{3,1}=[table1{3,1},num2str(bsm1prcValLP(i,1),'%02.2f'),char(9)];
        table1{4,1}=[table1{4,1},num2str(coffFreqHP(i,1),'%02.2f'),char(9)];
        table1{5,1}=[table1{5,1},num2str(bsm1prcValHP(i,1),'%02.2f'),char(9)];
    end
    disp(table1)
    
    
    
    
    %% Extract the errors with the noisy signals and the filtering recommendations

    medErrNoise=nan(18,1);
    iqrErrNoise=nan(18,1);
    for feat=1:1:nfeatures
        hplptmp=find(combinationsFg(:,1)==coffFreqLP(feat,1)&combinationsFg(:,2)==coffFreqHP(feat,1));
        if feat==17
            relErr_noise2=cell2mat(cellfun(@(x) x(feat:nfeatures:end,hplptmp),err_noise,'UniformOutput',0));
            acc=sum(abs(sign(relErr_noise2-1)))/sum(~isnan(relErr_noise2));
            medErrNoise(feat,1)=acc;
            iqrErrNoise(feat,1)=0;
        else
            relErr_noise2=cell2mat(cellfun(@(x) x(feat:nfeatures:end,hplptmp),relErr_noise,'UniformOutput',0));
            medErrNoise(feat,1)=nanmedian(abs(relErr_noise2))*100; %relative error in %
            iqrErrNoise(feat,1)=iqr(abs(relErr_noise2(~isnan(relErr_noise2))))*100; %relative error in %
        end
    end
    
    table3=cell(2,1);
    for feat=1:1:nfeatures
        table3{1,1}=[table3{1,1},'Feat',num2str(feat,'%02i'),char(9)];
        table3{2,1}=[table3{2,1},num2str(medErrNoise(feat,1),'%3.1f'),char(177),num2str(iqrErrNoise(feat,1),'%02.1f') ,char(9)];
    end
    disp(table3)
    
    
    %% Evaluate the noisy signals again but average beforehand the feature values on the noisy data
    
    medErrNoise=nan(18,1);
    iqrErrNoise=nan(18,1);
    for feat=1:1:nfeatures
        lphp=find(combinationsFg(:,1)==coffFreqLP(feat,1)&combinationsFg(:,2)==coffFreqHP(feat,1));
        if feat==17
            relErr_noise2=cell2mat(cellfun(@(x) x(feat:nfeatures:end,hplptmp),err_noise_mean,'UniformOutput',0));
            acc=sum(abs(sign(relErr_noise2-1)))/sum(~isnan(relErr_noise2));
            medErrNoise(feat,1)=acc;
            iqrErrNoise(feat,1)=0;
        else
            relErr_noise2=cell2mat(cellfun(@(x) x(feat:nfeatures:end,lphp),relErr_noise_mean,'UniformOutput',0));
            medErrNoise(feat,1)=nanmean(abs(relErr_noise2))*100; %relative error in %
            iqrErrNoise(feat,1)=iqr(abs(relErr_noise2(~isnan(relErr_noise2))))*100; %relative error in %
        end
    end
    
    table4=cell(2,1);
    for feat=1:1:nfeatures
        table4{1,1}=[table4{1,1},'Feat',num2str(feat,'%02i'),char(9)];
        table4{2,1}=[table4{2,1},num2str(medErrNoise(feat,1),'%3.1f'), char(177),num2str(iqrErrNoise(feat,1),'%02.1f') ,char(9)];
    end
    disp(table4)
end