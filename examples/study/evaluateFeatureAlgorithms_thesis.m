%% Preparation
clear
close all

%%
addpath(genpath('../../dependencies/ECGdeli'));
addpath(genpath('../../algorithms/'));
input=load('sim_signals.mat');

% create a local cluster object
pc = parcluster('local');

% get the number of dedicated cores from environment
num_workers = str2num(getenv('SLURM_NPROCS'));

% explicitly set the JobStorageLocation to the tmp directory that is unique to each cluster job (and is on local, fast scratch)
parpool_tmpdir = [getenv('TMP'),'/.matlab/local_cluster_jobs/slurm_jobID_',getenv('SLURM_JOB_ID')];
mkdir(parpool_tmpdir);
pc.JobStorageLocation = parpool_tmpdir;

% start the parallel pool
parpool(pc, num_workers, 'IdleTimeout',Inf);


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

for dBlvl=1:1:size(noiseLvl,2)
    load(strcat('Results_',num2str(noisedB),'dB.mat'))
    
    %% Visualize the results
    nRep=50;
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
    % Get the median errors for all leads and all LP-filterings (1-9)
    idxLP=find(combinationsFg(:,2)==min(combinationsFg(:,2)));
    meanErrNoise=nan(nfeatures,length(idxLP));
    stdErrNoise=nan(nfeatures,length(idxLP));
    
    for lp=1:1:length(idxLP)
        for feat=1:1:18 % features
            relErr_noise2=cell2mat(cellfun(@(x) x(feat:nfeatures:end,lp),relErr_noise,'UniformOutput',0));
            meanErrNoise(feat,lp)=nanmedian(relErr_noise2)*100; %relative error in %
            stdErrNoise(feat,lp)=iqr(relErr_noise2)*100; %relative error in %
        end
    end
    
    latextbl2=cell(3,1);
    for feat=1:1:nfeatures/3
        latextbl2{1,1}=[latextbl2{1,1},num2str(meanErrNoise(feat,9),'%3.1f'), '$\pm$',num2str(stdErrNoise(feat,9),'%02.1f') ,'&'];
    end
    for feat=nfeatures/3+1:1:2*nfeatures/3
        latextbl2{2,1}=[latextbl2{2,1},num2str(meanErrNoise(feat,9),'%3.1f'), '$\pm$',num2str(stdErrNoise(feat,9),'%02.1f') ,'&'];
    end
    for feat=2*nfeatures/3+1:1:nfeatures
        latextbl2{3,1}=[latextbl2{3,1},num2str(meanErrNoise(feat,9),'%3.1f'), '$\pm$',num2str(stdErrNoise(feat,9),'%02.1f') ,'&'];
    end
    disp(latextbl2)
    
    
    
    
    %
    latextbl=cell(8,1);
    % get the lowest LP frequency without disturbing the features
    idxLP=find(combinationsFg(:,2)==min(combinationsFg(:,2)));
    meanErrFilt=nan(size(featureMatrix_noise,3),nfeatures,length(idxLP));
    for ld=1:1:8 % leads
        for lp=1:1:length(idxLP)
            for feat=1:1:18 % features
                meanErrFilt(ld,feat,lp)=nanmean(relErr_cleanFilt{ld,1}(feat:nfeatures:end,lp))*100; %relative error in %
            end
        end
    end
    
    tmp=squeeze(max(abs(meanErrFilt),[],1));
    bsm1prc=nan(nfeatures,1);
    bsm1prcValLP=nan(nfeatures,1);
    for ft=1:1:nfeatures
        if ~isempty(find(tmp(ft,:)<=5,1,'first'))
            bsm1prc(ft,1)=find(tmp(ft,:)<=5,1,'first');
            bsm1prcValLP(ft,1)=tmp(ft,bsm1prc(ft,1));
        end
    end
    %bsm1prc(12,1)=1;
    bsm1prc(17,1)=1;
    coffFreqLP=combinationsFg(bsm1prc,1);
    
    for i=1:1:length(bsm1prc)/2
        latextbl{1,1}=[latextbl{1,1},num2str(coffFreqLP(i,1),'%02i'),'&'];
        latextbl{2,1}=[latextbl{2,1},num2str(bsm1prcValLP(i,1),'%02.2f'),'&'];
    end
    
    for i=length(bsm1prc)/2:1:length(bsm1prc)
        latextbl{5,1}=[latextbl{5,1},num2str(coffFreqLP(i,1),'%02i'),'&'];
        latextbl{6,1}=[latextbl{6,1},num2str(bsm1prcValLP(i,1),'%02.2f'),'&'];
    end
    
    % get the highest HP frequency without disturbing the features
    idxLP=find(combinationsFg(:,1)==max(combinationsFg(:,1)));
    meanErrFilt=nan(size(featureMatrix_noise,3),nfeatures,length(idxLP));
    for ld=1:1:8 % leads
        for lp=1:1:length(idxLP)
            for feat=1:1:18 % features
                meanErrFilt(ld,feat,lp)=nanmean(relErr_cleanFilt{ld,1}(feat:nfeatures:end,lp))*100; %relative error in %
            end
        end
    end
    tmp=squeeze(max(abs(meanErrFilt),[],1));
    bsm1prc=nan(nfeatures,1);
    bsm1prcValHP=nan(nfeatures,1);
    for ft=1:1:nfeatures
        if ~isempty(find(tmp(ft,:)<=5,1,'first'))
            bsm1prc(ft,1)=find(tmp(ft,:)<=5,1,'first');
            bsm1prcValHP(ft,1)=tmp(ft,bsm1prc(ft,1));
        end
    end
    bsm1prc(17,1)=1;
    coffFreqHP=combinationsFg(idxLP(bsm1prc),2);
    for i=1:1:length(bsm1prc)/2
        latextbl{3,1}=[latextbl{3,1},num2str(coffFreqHP(i,1),'%02.2f'),'&'];
        latextbl{4,1}=[latextbl{4,1},num2str(bsm1prcValHP(i,1),'%02.2f'),'&'];
    end
    
    for i=length(bsm1prc)/2:1:length(bsm1prc)
        latextbl{7,1}=[latextbl{7,1},num2str(coffFreqHP(i,1),'%02.2f'),'&'];
        latextbl{8,1}=[latextbl{8,1},num2str(bsm1prcValHP(i,1),'%02.2f'),'&'];
    end
    disp(latextbl)
    
    
    
    
    %
    
    nleads=8;
    % get the lowest LP frequency without disturbing the features
    idxLP=find(combinationsFg(:,2)==min(combinationsFg(:,2)));
    meanErrFilt=nan(size(featureMatrix_noise,3),nfeatures,length(idxLP));
    for ld=1:1:nleads % leads
        for lp=1:1:length(idxLP)
            for feat=1:1:18 % features
                meanErrFilt(ld,feat,lp)=nanmean(relErr_cleanFilt{ld,1}(feat:nfeatures:end,lp))*100; %relative error in %
            end
        end
    end
    
    tmp=squeeze(max(abs(meanErrFilt),[],1));
    bsm1prcLP=nan(nfeatures,1);
    for ft=1:1:nfeatures
        if ~isempty(find(tmp(ft,:)<=5,1,'first'))
            bsm1prcLP(ft,1)=find(tmp(ft,:)<=5,1,'first');
        end
    end
    %bsm1prc(12,1)=1;
    bsm1prcLP(17,1)=1;
    coffFreqLP=combinationsFg(idxLP(bsm1prcLP),1);
    
    
    % get the highest HP frequency without disturbing the features
    idxHP=find(combinationsFg(:,1)==max(combinationsFg(:,1)));
    meanErrFilt=nan(size(featureMatrix_noise,3),nfeatures,length(idxHP));
    for ld=1:1:nleads % leads
        for lp=1:1:length(idxHP)
            for feat=1:1:18 % features
                meanErrFilt(ld,feat,lp)=nanmean(relErr_cleanFilt{ld,1}(feat:nfeatures:end,lp))*100; %relative error in %
            end
        end
    end
    tmp=squeeze(max(abs(meanErrFilt),[],1));
    bsm1prc=nan(nfeatures,1);
    bsm1prcValHP=nan(nfeatures,1);
    for ft=1:1:nfeatures
        if ~isempty(find(tmp(ft,:)<=5,1,'first'))
            bsm1prc(ft,1)=find(tmp(ft,:)<=5,1,'first');
            bsm1prcValHP(ft,1)=tmp(ft,bsm1prc(ft,1));
        end
    end
    bsm1prc(17,1)=1;
    coffFreqHP=combinationsFg(idxHP(bsm1prc),2);
    
    meanErrNoise=nan(18,1);
    stdErrNoise=nan(18,1);
    for feat=1:1:nfeatures
        lp(feat)=find(combinationsFg(:,1)==coffFreqLP(feat,1)&combinationsFg(:,2)==coffFreqHP(feat,1));
        relErr_noise2=cell2mat(cellfun(@(x) x(feat:nfeatures:end,lp(feat)),relErr_noise,'UniformOutput',0));
        meanErrNoise(feat,1)=nanmean(abs(relErr_noise2))*100; %relative error in %
        stdErrNoise(feat,1)=nanstd(abs(relErr_noise2))*100; %relative error in %
    end
    
    latextbl3=cell(3,1);
    for feat=1:1:nfeatures/3
        latextbl3{1,1}=[latextbl3{1,1},num2str(meanErrNoise(feat,1),'%3.1f'), '$\pm$',num2str(stdErrNoise(feat,1),'%02.1f') ,'&'];
    end
    for feat=nfeatures/3+1:1:2*nfeatures/3
        latextbl3{2,1}=[latextbl3{2,1},num2str(meanErrNoise(feat,1),'%3.1f'), '$\pm$',num2str(stdErrNoise(feat,1),'%02.1f') ,'&'];
    end
    for feat=2*nfeatures/3+1:1:nfeatures
        latextbl3{3,1}=[latextbl3{3,1},num2str(meanErrNoise(feat,1),'%3.1f'), '$\pm$',num2str(stdErrNoise(feat,1),'%02.1f') ,'&'];
    end
    disp(latextbl3)
    
    
    
    
    
    
    
    
    % again with mean value of every iteration
    
    nleads=8;
    % get the lowest LP frequency without disturbing the features
    idxLP=find(combinationsFg(:,2)==min(combinationsFg(:,2)));
    meanErrFilt=nan(size(featureMatrix_noise,3),nfeatures,length(idxLP));
    for ld=1:1:nleads % leads
        for lp=1:1:length(idxLP)
            for feat=1:1:18 % features
                meanErrFilt(ld,feat,lp)=nanmean(relErr_cleanFilt{ld,1}(feat:nfeatures:end,lp))*100; %relative error in %
            end
        end
    end
    
    tmp=squeeze(max(abs(meanErrFilt),[],1));
    bsm1prcLP=nan(nfeatures,1);
    for ft=1:1:nfeatures
        if ~isempty(find(tmp(ft,:)<=5,1,'first'))
            bsm1prcLP(ft,1)=find(tmp(ft,:)<=5,1,'first');
        end
    end
    %bsm1prc(12,1)=1;
    bsm1prcLP(17,1)=1;
    coffFreqLP=combinationsFg(idxLP(bsm1prcLP),1);
    
    
    % get the highest HP frequency without disturbing the features
    idxHP=find(combinationsFg(:,1)==max(combinationsFg(:,1)));
    meanErrFilt=nan(size(featureMatrix_noise,3),nfeatures,length(idxHP));
    for ld=1:1:nleads % leads
        for lp=1:1:length(idxHP)
            for feat=1:1:18 % features
                meanErrFilt(ld,feat,lp)=nanmean(relErr_cleanFilt{ld,1}(feat:nfeatures:end,lp))*100; %relative error in %
            end
        end
    end
    tmp=squeeze(max(abs(meanErrFilt),[],1));
    bsm1prc=nan(nfeatures,1);
    bsm1prcValHP=nan(nfeatures,1);
    for ft=1:1:nfeatures
        if ~isempty(find(tmp(ft,:)<=5,1,'first'))
            bsm1prc(ft,1)=find(tmp(ft,:)<=5,1,'first');
            bsm1prcValHP(ft,1)=tmp(ft,bsm1prc(ft,1));
        end
    end
    bsm1prc(17,1)=1;
    coffFreqHP=combinationsFg(idxHP(bsm1prc),2);
    
    meanErrNoise=nan(18,1);
    stdErrNoise=nan(18,1);
    for feat=1:1:nfeatures
        lp(feat)=find(combinationsFg(:,1)==coffFreqLP(feat,1)&combinationsFg(:,2)==coffFreqHP(feat,1));
        relErr_noise2=cell2mat(cellfun(@(x) x(feat:nfeatures:end,lp(feat)),relErr_noise_mean,'UniformOutput',0));
        meanErrNoise(feat,1)=nanmean(abs(relErr_noise2))*100; %relative error in %
        stdErrNoise(feat,1)=nanstd(abs(relErr_noise2))*100; %relative error in %
    end
    
    latextbl4=cell(3,1);
    for feat=1:1:nfeatures/3
        latextbl4{1,1}=[latextbl4{1,1},num2str(meanErrNoise(feat,1),'%3.1f'), '$\pm$',num2str(stdErrNoise(feat,1),'%02.1f') ,'&'];
    end
    for feat=nfeatures/3+1:1:2*nfeatures/3
        latextbl4{2,1}=[latextbl4{2,1},num2str(meanErrNoise(feat,1),'%3.1f'), '$\pm$',num2str(stdErrNoise(feat,1),'%02.1f') ,'&'];
    end
    for feat=2*nfeatures/3+1:1:nfeatures
        latextbl4{3,1}=[latextbl4{3,1},num2str(meanErrNoise(feat,1),'%3.1f'), '$\pm$',num2str(stdErrNoise(feat,1),'%02.1f') ,'&'];
    end
    disp(latextbl4)
end