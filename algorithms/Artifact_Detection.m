
function [FPT]=Artifact_Detection(Filtered_Values,samplerate,FPT)
%This is a function used to detect if beats classified as ectopic are
%actually artifacts.
%Inputs:
%Filtered_Values: The ECG signal that was used to classify the beats. It
%should be a good looking ecg signal.
%samplerate: Sample Frequency of the signal
%FPT: Fidutial Point Table corresponding the the ECG signal. The beats must
%be classified.
%Outputs:
%FPT: Fidutial Point Table containing the reclassified beats as artifacts. 

disp('Detecting Artifacts...');

%% Downsampling
%Downsampling to accelerate artifact detection. The sample frequency used
%is 400Hz
fdownsample=400;
flagdownsample=false;
if samplerate>fdownsample;
    length_signal=round(length(Filtered_Values)*fdownsample/samplerate);
    k=round(linspace(1,length(Filtered_Values),length_signal)); %sample values for downsampling
    Filtered_Values=Filtered_Values(k); %Downsampling
    RPOS=FPT(:,6); %Save position of R peaks before downsampling for restitution in the end
    FPT(:,6)=round(FPT(:,6)*fdownsample/samplerate);
    samplerate=fdownsample;
    [Template_ecg,posRpeak,ampRpeak]=Create_Template(Filtered_Values,samplerate,FPT,'ECG'); %Create Template for comparison later on
    L1=posRpeak-1;%Length of Template before Rpeak
    L2=length(Template_ecg)-posRpeak; %Length of Template after R peak
    MTA=ampRpeak; 
    flagdownsample=true;
else
    [Template_ecg,posRpeak,ampRpeak]=Create_Template(Filtered_Values,samplerate,FPT,'ECG'); %Create Template for comparison later on
    L1=posRpeak-1;
    L2=length(Template_ecg)-posRpeak;
    MTA=ampRpeak;
end

%% Variables
%Create variables for later use.
ES=find(FPT(:,13)==2 | FPT(:,13)==3);
NE_vector=zeros(length(ES),1);
SE_vector=zeros(length(ES),1);

ML_matrix=zeros(length(ES),2);
ML_vector=zeros(length(ES),1);
MC_matrix=zeros(size(ML_matrix));
MC_vector=zeros(length(ES),1);

%Create thresholds for the artifact detection
SE_threshold=80;
NE_threshold=2.75;
ML_threshold=0.8;
MC_threshold=0.8;

%% Artifact Detection
%Check if ectopic beats are actually artifacts
if ~isempty(Template_ecg)
for i=1:length(ES)
       
    %Every beat classified as ectopic is segmented. T contains the RR
    %intervals before and after every beat classified as ectopic.
    T=[FPT(ES(i)-1,6);FPT(ES(i),6);FPT(ES(i)+1,6)];
    dT=diff(T);
    L3=round(dT(1)*1/3);%Length before the QRS complex
    L4=round(dT(2)*2/3);%Length after the QRS complex
    signal_ecg=Filtered_Values(FPT(ES(i),6)-L3:FPT(ES(i),6)+L4);%segmented beat marked as ectopic
    
    [signal_ecg]=Isoline_Correction(signal_ecg);%remove offset from beat marked as ectopic
    signal_ecg=signal_ecg/MTA; %normalize amplitude
    
    %Equalize the lengths of the Template and the segmented beat
    if L3<L1
        Fake_ecg1=Template_ecg(L1-L3+1:1:L1+1);
        signal_ecg1=signal_ecg(1:L3+1);
    elseif L3>L1
        Fake_ecg1=Template_ecg(1:L1+1);
        signal_ecg1=signal_ecg(L3-L1+1:L3+1);
    elseif L3==L1
        Fake_ecg1=Template_ecg(1:L1+1);
        signal_ecg1=signal_ecg(1:L3+1);
    end
    if L4<L2
        Fake_ecg2=Template_ecg(L1+2:L1+1+L4);
        signal_ecg2=signal_ecg(L3+2:L3+1+L4);
    elseif L4>L2
        Fake_ecg2=Template_ecg(L1+2:L1+1+L2);
        signal_ecg2=signal_ecg(L3+2:L3+1+L2);
    elseif L4==L2
        signal_ecg2=signal_ecg(L3+2:L3+1+L4);
        Fake_ecg2=Template_ecg(L1+2:L1+1+L2);
    end
    %Create Template and beat marked as ectopic with same length
    Fake_ecg=[Fake_ecg1;Fake_ecg2];
    signal_ecg=[signal_ecg1;signal_ecg2];
    
    %power of beat compared to power of template
    SE=sum(signal_ecg.^2)/sum(Fake_ecg.^2);
    SE_vector(i)=SE;%Power ratio
    
    %Number of extrema in the template compared to the number of extrema in
    %the beat marked as ectopic
    DFake_ecg=diff(Fake_ecg);
    k=1:length(DFake_ecg)-1;
    NE=(DFake_ecg(k)>=0 & DFake_ecg(k+1)<0) | (DFake_ecg(k)<0 & DFake_ecg(k+1)>=0);
    NE=sum(NE); %Number of extrema in the template        
    
    Dsignal_ecg=diff(signal_ecg);
    NEs=(Dsignal_ecg(k)>=0 & Dsignal_ecg(k+1)<0) | (Dsignal_ecg(k)<0 & Dsignal_ecg(k+1)>=0);
    NEs=sum(NEs); %Number of extrema in the beat marked as ectopic
        
    NE_vector(i)=NEs/NE;%Ratio of number of extrema
    
    %If the beat marked as ectopic is classified as artifact at this point,
    %no more investigations are needed
    if SE_vector(i)>=SE_threshold || NE_vector(i)>=NE_threshold
        continue
    end
    
    %Investigation of neighboring normal beats
    for j=1:2
        %Find first and last normal beat before and after the beat marked
        %as ectopic
        K=[-1,1];
        K1flag=false;
        K2flag=false;
        while FPT(ES(i)+K(1),13)~=1 
            if ES(i)+K(1)==1 %special consideration for the beginning of the signal
                K1flag=true;
                break
            else
                K(1)=K(1)-1;
            end
        end
        while FPT(ES(i)+K(2),13)~=1 
            if ES(i)+K(2)==size(FPT,1) %special consideration for the end of the signal
                K2flag=true;
                break
            else
                K(2)=K(2)+1;
            end
        end
        %For beats marked as ectopic at the beginning and end of signal
        %only one neighboring beat is used
        if K1flag
            K(1)=K(2);
        end
        if K2flag
            K(2)=K(1);
        end
        if K1flag && K2flag
            continue
        end
        
        %Lengths for segmentation of the neighboring normal beat
        T=[FPT(ES(i)+K(j)-1,6);FPT(ES(i)+K(j),6);FPT(ES(i)+K(j)+1,6)];
        dT=diff(T);
        L3=round(dT(1)*0.35);%Length before the QRS complex
        L4=round(dT(2)*0.55);%Length after the QRS complex
        signal_ecg=Filtered_Values(FPT(ES(i)+K(j),6)-L3:FPT(ES(i)+K(j),6)+L4);%neighboring normal beat
        
        [signal_ecg]=Isoline_Correction(signal_ecg);%remove offset from normal beat neighboring the ectopic
        signal_ecg=signal_ecg/MTA;%normalize amplitude
        
        %Equalize the lengths of the Template and the segmented neighboring normal beat
        if L3<L1
            Fake_ecg1=Template_ecg(L1-L3+1:1:L1+1);
            signal_ecg1=signal_ecg(1:L3+1);
        elseif L3>L1
            Fake_ecg1=Template_ecg(1:L1+1);
            signal_ecg1=signal_ecg(L3-L1+1:L3+1);
        elseif L3==L1
            Fake_ecg1=Template_ecg(1:L1+1);
            signal_ecg1=signal_ecg(1:L3+1);
        end
        if L4<L2
            Fake_ecg2=Template_ecg(L1+2:L1+1+L4);
            signal_ecg2=signal_ecg(L3+2:L3+1+L4);
        elseif L4>L2
            Fake_ecg2=Template_ecg(L1+2:L1+1+L2);
            signal_ecg2=signal_ecg(L3+2:L3+1+L2);
        elseif L4==L2
            signal_ecg2=signal_ecg(L3+2:L3+1+L4);
            Fake_ecg2=Template_ecg(L1+2:L1+1+L2);
        end
        
        %Create Template and beat marked as ectopic with same lengths.
        %Every signal is created by joining left and right sides
        Fake_ecg=[Fake_ecg1;Fake_ecg2];
        signal_ecg=[signal_ecg1;signal_ecg2];        
        
        %power of diference signal between normal beat and template 
        %compared to power of template   
        ML=l_operator(signal_ecg,Fake_ecg);  
        
        %correlation between normal beat and template 
        MC=corr(signal_ecg,Fake_ecg);
        
        %Matrices containing DSE and MC for the normal beat before and
        %after the beat marked as ectopic
        ML_matrix(i,j)=ML;
        MC_matrix(i,j)=MC;
        
    end
    
    %The worst values are chose for the classification
    ML_vector(i)=min(ML_matrix(i,:));
    MC_vector(i)=min(MC_matrix(i,:));
    
end

%% Classification 
%The beats are labeled as artifacts if they do not fulfill the following
%threschholds 
index1=find(SE_vector>=SE_threshold | NE_vector>=NE_threshold);
index2=find(ML_vector<=ML_threshold |  MC_vector<=MC_threshold);
index=[index1;index2];
FPT(ES(index),13)=20;
disp(['Number of Artifacts Detected: ',num2str(length(unique(index))),'/',num2str(length(ES))]);

%Restitute positions of Rpeaks with original sample rate
if flagdownsample
    FPT(:,6)=RPOS;
end
else
    disp('Please check the signal...')
    FPT(:,13)=20;
end
disp('Done');
