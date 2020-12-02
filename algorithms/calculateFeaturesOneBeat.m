% -------------------------------------------------------
%
%    calculateFeaturesOneBeat.m  - Calculates several features from single
%    beat ECGs
%
%    Ver. 1.0.0
%
%    Created:           Nicolas Pilia (29.11.2018)
%    Last modified:     Nicolas Pilia (27.11.2020)
%
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2020 - All rights reserved.
%
% ------------------------------------------------------
%   INPUTS:
%       ecg_matrix: matrix (NxL), where N is
%       the temporal dimension and L
%           is the lead dimension
%       samplerate: samplerate of the signal
%       checkForBiphasic: flag  1: the template is scanned for biphasic T waves
%                               0: T waves are treated as they were
%                               monophasic
%       FPT: Fiducial points of the template; Here, the position of the R
%       peak and of the T peak are needed in the form of the FPT table like
%       in ECGdeli. If the algorithm should check for biphasic waves, the
%       position of P peak is also needed.
%
%	OUTPUTS:
%       featureMatrix: Values of the features
%       featureMatrix: Names of the features
%

function [featureMatrix,featureNames,RPos]=calculateFeaturesOneBeat(ecg_matrix,samplerate,checkForBiphasic,FPT)
featureMatrix=zeros(18,length(ecg_matrix(1,:)));
flagBi=nan(size(ecg_matrix,2),1);
if size(ecg_matrix,2)==1 && ~iscell(FPT)
    FPTcell=cell(1,1);
    FPTcell{1,1}=FPT;
elseif size(ecg_matrix,2) == size(FPT,1)
    FPTcell=FPT;
else
    error('FPT dimension is not in accordance with lead number!')
end


%% If there are more leads, go through them
for setupLead=1:1:size(ecg_matrix,2)
    % Preprocessing: Get T peak and R peak
    
    TPos=FPTcell{setupLead}(1,11);
    RPos=FPTcell{setupLead}(1,6);
    TAmp=ecg_matrix(TPos,setupLead);
    if checkForBiphasic==1
        [flagBi(setupLead,1),Tpos2]=checkForBiphasicT(ecg_matrix(:,setupLead),samplerate,FPTcell{setupLead,1});
    else
        flagBi(setupLead,1)=0;
    end
    
    
    if flagBi(setupLead)==1
        TPosBi=[TPos;Tpos2];
        [TPosBi,idx]=sort(TPosBi);
        TAmpBi=[TAmp;ecg_matrix(Tpos2,setupLead)];
        TAmpBi=TAmpBi(idx);
        [featureMatrix(:,setupLead),featureNames]=calculateFeaturesOneBeat_biphasic(ecg_matrix(:,setupLead),samplerate,TPosBi,RPos,TAmpBi);
    else
        [featureMatrix(:,setupLead),featureNames]=calculateFeaturesOneBeat_monophasic(ecg_matrix(:,setupLead),samplerate,TPos,RPos,TAmp);
    end
    
end
end

%% Function for the detection of biphasic T waves
function [flagBi,Tpos2]=checkForBiphasicT(signal,samplerate,FPT)
flagBi=0;
Tpos2=[];
extendedSignal=repmat(signal,3,1);
filteredExtendedSignal=ECG_High_Low_Filter(extendedSignal,samplerate,0.6,30,'B');
filteredExtendedSignal=filteredExtendedSignal(size(signal,1)+1:size(signal,1)*2);

% Detect PR
[~,baseline]=ECG_Baseline_Removal(filteredExtendedSignal,samplerate,0.10,0.5);
constSec=(abs(diff(baseline)));
[~,PR]=min(constSec(FPT(1,2)+round(0.03*samplerate):FPT(1,6)-0.05*samplerate));
PR=PR+FPT(1,2)+round(0.03*samplerate)-1;
if isempty(PR)
    PR=11;
end
filteredExtendedSignal=filteredExtendedSignal-mean(filteredExtendedSignal(PR-10:PR+10));

% Check for bisphasic waves
[~,maxPos,~,promMax]=findpeaks(filteredExtendedSignal);
[~,minPos,~,promMin]=findpeaks(-filteredExtendedSignal);

if min(abs(maxPos-FPT(:,11))) > min(abs(minPos-FPT(:,11))) %i.e. T peak is a minimum -> we have to look for a maximum
    [~,filteredPeakPos]=min(abs(minPos-FPT(:,11)));
    filteredPeakProm=promMin(filteredPeakPos);
    filteredPeakPos=minPos(filteredPeakPos);
    
    
    Tcand=maxPos(maxPos<=filteredPeakPos+0.2*samplerate & maxPos>=max(filteredPeakPos-0.2*samplerate,round(FPT(:,6)+0.05*samplerate)));
    TcandProm=promMax(maxPos<=filteredPeakPos+0.2*samplerate & maxPos>=max(filteredPeakPos-0.2*samplerate,round(FPT(:,6)+0.05*samplerate)));
    [~,T2peakPos]=max(filteredExtendedSignal(Tcand));
    T2peakProm=TcandProm(T2peakPos);
    T2peakPos=Tcand(T2peakPos);
   
    
    if ~isempty(T2peakPos)
        [pAmp,pPos,~,pProm]=findpeaks(signal);
        pPos(pPos==filteredPeakPos)=[];
        tmp=find(abs(pPos-T2peakPos)<0.01*samplerate);
        [~,realT2peakPos]=max(pAmp(tmp));
        %realT2peakProm=pProm(tmp(realT2peakPos));
        realT2peakPos=pPos(tmp(realT2peakPos));
    else
        %realT2peakProm=[];
        realT2peakPos=[];
    end
    
else %i.e. T peak is a maximum -> we have to look for a minimum
    [~,filteredPeakPos]=min(abs(maxPos-FPT(:,11)));
    filteredPeakProm=promMax(filteredPeakPos);
    filteredPeakPos=maxPos(filteredPeakPos);
    
    Tcand=minPos(minPos<=filteredPeakPos+0.2*samplerate & minPos>=max(filteredPeakPos-0.2*samplerate,round(FPT(:,6)+0.05*samplerate)));
    TcandProm=promMin(minPos<=filteredPeakPos+0.2*samplerate & minPos>=max(filteredPeakPos-0.2*samplerate,round(FPT(:,6)+0.05*samplerate)));
    [~,T2peakPos]=min(filteredExtendedSignal(Tcand));
    T2peakProm=TcandProm(T2peakPos);
    T2peakPos=Tcand(T2peakPos);
    
    if ~isempty(T2peakPos)
        [pAmp,pPos,~,pProm]=findpeaks(-signal);
        pPos(pPos==filteredPeakPos)=[];
        tmp=find(abs(pPos-T2peakPos)<0.05*samplerate);
        [~,realT2peakPos]=max(pAmp(tmp));
        %realT2peakProm=pProm(tmp(realT2peakPos));
        realT2peakPos=pPos(tmp(realT2peakPos));
    else
        realT2peakPos=[];
    end
    
end
if ~isempty(realT2peakPos) && (sign(filteredExtendedSignal(filteredPeakPos,1)) ~= sign(filteredExtendedSignal(T2peakPos,1))) %% vielleicht korrigieren?!?
    if abs(filteredExtendedSignal(filteredPeakPos,1)/filteredExtendedSignal(T2peakPos,1))<2.5 && abs(filteredExtendedSignal(filteredPeakPos,1)/filteredExtendedSignal(T2peakPos,1))>0.4
        if T2peakProm/filteredPeakProm>0.5 && T2peakProm/filteredPeakProm<2 
        flagBi=1;
        Tpos2=realT2peakPos;
        end
    end
end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monophasic implememtation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [featureMatrix,featureNames]=calculateFeaturesOneBeat_monophasic(ecg_matrix,samplerate,TPos,RPos,TAmp)
featureNames={'center'; 'variance'; 'skewness'; 'curtosis'; 'RT distance'; 'RTmid distance';...
    'Peakness T wave'; 'T amplitude'; 'T upslope'; 'T downslope'; 'left half T wave area ratio';...
    'right half T wave area ratio'; 'R Amplitude';'R power';'ratio R power/R amplitude';'ST elevation';'Biphasic';'AUC T Wave/AUC ECG Beat'};
%%
signal=ecg_matrix;
featureMatrix=nan(size(featureNames,1),1);

%% Calculate statistical features: Center, Variance, Skewness and Curtosis
GIdxLeft=max(1,TPos-round(150*10^(-3)*samplerate));
GIdxRight=min(TPos+round(150*10^(-3)*samplerate),size(signal,1));
t=(0:1/samplerate:(length(GIdxLeft:GIdxRight)-1)/samplerate)';
fittedT=signal(GIdxLeft:GIdxRight).^2/trapz(signal(GIdxLeft:GIdxRight).^2);
featureMatrix(1,1)=trapz(t,t.*fittedT); %center
featureMatrix(2,1)=trapz(t,(t-featureMatrix(1,1)).^2.*fittedT); %variance
featureMatrix(3,1)=trapz(t,(t-featureMatrix(1,1)).^3.*fittedT); %Schiefe/skewness
featureMatrix(4,1)=trapz(t,(t-featureMatrix(1,1)).^4.*fittedT); %Wölbung/curtosis
featureMatrix(18,1)=trapz(signal(GIdxLeft:GIdxRight).^2)/trapz(signal.^2); %AUC T wave

%% Calculate temporal features: R-T-distance
% For monophasic waves, RT and RTmid are equal
featureMatrix(5,1)=TPos-RPos;
featureMatrix(6,1)=TPos-RPos;

%% Fitting two Gaussians as a surrogate for the peakedness
boundLeft=max(1,TPos(1)-round(150*10^(-3)*samplerate));
boundRight=min(TPos(end)+round(150*10^(-3)*samplerate),size(signal,1));

fitsig=signal(boundLeft:boundRight,1); % normalisation
fitsig=fitsig-min(fitsig); % normalisation
fitsig=fitsig/max(fitsig); % normalisation
xfit=linspace(0,1,length(fitsig)); % normalised time vector

pkPos=xfit(round(150*10^(-3)*samplerate)+1);
if TPos(1)-round(150*10^(-3)*samplerate)>0
    npkPos=round(round(150*10^(-3)*samplerate))+1;
else
    npkPos=TPos(1);
end
    
% First Gaussian on the descending part (right part)
xexcl=ones(length(xfit),1);
xexcl((round(npkPos-20*10^(-3)*samplerate):round(npkPos+round(120*10^(-3)*samplerate)))',1)=0;
xexcl=find(xexcl);
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[-1,-2.001,0,-2],...
    'Upper',[1,2,1,2],...
    'StartPoint',[0.001,sign(signal(TPos)),pkPos,0],...
    'Algorithm','Trust-Region',...
    'Exclude',xexcl);
ft = fittype('(k/sigma^2*exp(-(x-mu)^2/(2*sigma^2)))+h','independent','x','coefficients',{'sigma','k','mu','h'});
curve1 = fit(xfit',fitsig,ft,fo);

% Second Gaussian on the ascending part (left part)
xexcl=ones(length(xfit),1);
xexcl((round(npkPos-round(120*10^(-3)*samplerate)):round(npkPos+20*10^(-3)*samplerate))',1)=0;
xexcl=find(xexcl);
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[-1,-2.001,0,-2],...
    'Upper',[1,2,1,2],...
    'StartPoint',[0.001,sign(signal(TPos)),pkPos,0],...
    'Algorithm','Trust-Region',...
    'Exclude',xexcl);
ft = fittype('(k/sigma^2*exp(-(x-mu)^2/(2*sigma^2)))+h','independent','x','coefficients',{'sigma','k','mu','h'});
curve2 = fit(xfit',fitsig,ft,fo);

% feature is the sum of both
featureMatrix(7,1)=curve1.sigma*length(fitsig)+curve2.sigma*length(fitsig);
featureMatrix(16,1)=curve2.h-curve1.h;


% T Amplitude
featureMatrix(8,1)=TAmp;

% Calculate upslope
% cut a piece of the wave
boundRight=min(length(signal),TPos(1)+round(15*10^(-3)*samplerate));
boundLeft=max(1,TPos(1)-round(80*10^(-3)*samplerate));
xsl=(boundLeft:boundRight)';
% fit a polynomial of 4th order
ysl=signal(xsl);
xsl=xsl-xsl(round(length(xsl)/2));
Xsl=[xsl.^4,xsl.^3,xsl.^2,xsl,ones(length(xsl),1)];
up4thGrade = lsqminnorm(Xsl,ysl);
yslE=Xsl*up4thGrade;
minNorm=min(yslE(1:length(yslE)-round(15*10^(-3)*samplerate)));
maxNorm=max(yslE(1:length(yslE)-round(15*10^(-3)*samplerate)));
yslEnorm=(yslE-minNorm)./(maxNorm-minNorm);
yidx=find(yslEnorm<0.9 & yslEnorm>=0.3);
yidx=yidx(yidx<length(yslE)-round(15*10^(-3)*samplerate));
featureMatrix(9,1)=mean(diff(yslE(yidx,1)));

% Calculate downslope
% cut a piece of the wave
boundRight=min(TPos(end)+round(80*10^(-3)*samplerate),size(signal,1));
boundLeft=max(1,TPos(end)-15);
% fit a polynomial of 4th order
xsr=(boundLeft:boundRight)';
ysr=signal(xsr);
xsr=xsr-xsr(round(length(xsr)/2));
Xsr=[xsr.^4,xsr.^3,xsr.^2,xsr,ones(length(xsr),1)];
down4thGrade = lsqminnorm(Xsr,ysr);
ysrE=Xsr*down4thGrade;
minNorm=min(ysrE(round(15*10^(-3)*samplerate):end));
maxNorm=max(ysrE(round(15*10^(-3)*samplerate):end));
ysrEnorm=(ysrE-minNorm)./(maxNorm-minNorm);
yidx=find(ysrEnorm<0.9 & ysrEnorm>=0.3);
yidx=yidx(yidx>round(15*10^(-3)*samplerate));
featureMatrix(10,1)=mean(diff(ysrE(yidx,1)));

%% Fraction of left/right half at the area of the T wave
leftIdx=max(1,TPos-round(150*10^(-3)*samplerate));
rightIdx=min(TPos+round(150*10^(-3)*samplerate),size(signal,1));
featureMatrix(11,1)=sum(signal(leftIdx:TPos-1))/sum(signal(leftIdx:rightIdx));
featureMatrix(12,1)=sum(signal(TPos+1:rightIdx))/sum(signal(leftIdx:rightIdx));

%% QRS features: Amplitude, Energy, Ratio
featureMatrix(13,1)=signal(RPos);
RIdxLeft=max(1,RPos-0.035*samplerate);
RIdxRight=min(RPos+0.035*samplerate,size(signal,1));
featureMatrix(14,1)=sum(signal(RIdxLeft:RIdxRight).^2);
featureMatrix(15,1)=featureMatrix(14,1)/featureMatrix(13,1);

%% Biphasic parameter is 0
featureMatrix(17,1)=0;
end


%%%%%%%%%%%%%%%%%%%%%
% Biphasic features %
%%%%%%%%%%%%%%%%%%%%%
function [featureMatrix,featureNames]=calculateFeaturesOneBeat_biphasic(ecg_matrix,samplerate,TPos,RPos,TAmp)
featureNames={'center'; 'variance'; 'skewness'; 'curtosis'; 'RT distance'; 'RTmid distance';...
    'Peakness T wave'; 'T amplitude'; 'T upslope'; 'T downslope'; 'left half T wave area ratio';...
    'right half T wave area ratio'; 'R Amplitude';'R power';'ratio R power/R amplitude';'ST elevation';'Biphasic';'AUC T Wave/AUC ECG Beat'};

%%
featureMatrix=nan(18,1);
signal=ecg_matrix(:,1);

TmaxAmp=TAmp(1,1);
TminAmp=TAmp(2,1);
TmaxAmpPos=TPos(1,1);
TminAmpPos=TPos(2,1);

% find the inflection point between 2 peaks
boundRight=min(TPos(2)+round(15*10^(-3)*samplerate),size(signal,1));
boundLeft=max(RPos+round(0.05*samplerate),TPos(1)-round(15*10^(-3)*samplerate));
% fit a polynomial of 4th order
xsr=(boundLeft:boundRight)';
ysr=signal(xsr);
xsr=xsr-xsr(round(length(xsr)/2));
Xsr=[xsr.^4,xsr.^3,xsr.^2,xsr,ones(length(xsr),1)];
down4thGrade = lsqminnorm(Xsr,ysr);
ysrE=Xsr*down4thGrade;
d2T=diff(diff(ysrE));
TzeroPos=find(d2T(1:end-1,1)>0 & d2T(2:end,1)<=0)+TPos(1)-round(15*10^(-3)*samplerate)+1;
if isempty(TzeroPos)
    TzeroPos=find(d2T(1:end-1,1)<0 & d2T(2:end,1)>=0)+TPos(1)-round(15*10^(-3)*samplerate)+1;
end

if isempty(TzeroPos) %Obviously, this wave is not a biphasic wave...
    [TAmp,tidx]=max(abs(TAmp));
    TPos=TPos(tidx);
    [featureMatrix,featureNames]=calculateFeaturesOneBeat_monophasic(ecg_matrix,samplerate,QRSmaxTransform,TPos,RPos,TAmp);
else
    
    %% Fit Gaussian function for the calculation of features
    if length(TPos)>1
        TPosTmp=round((TPos(2)-TPos(1))/2+TPos(1));
    else
        GIdxLeft=max(1,TPos-round(150*10^(-3)*samplerate));
        GIdxRight=min(TPos+round(150*10^(-3)*samplerate),size(signal,1));
        leftM=find(abs(signal(GIdxLeft:TPos-1))<abs(0.5*TAmp(1)),1,'last')+TPos-round(150*10^(-3)*samplerate)-1;
        rightM=find(abs(signal(TPos+1:GIdxRight))<abs(0.5*TAmp(1)),1,'first')+TPos;
        
        if ~(isempty(leftM)||isempty(rightM))
            TPosTmp=round((rightM-leftM)/2+leftM);
        else
            TPosTmp=TPos(1);
        end
    end
    
    
    % Calculate statistical features: Center, Variance, Skewness and Curtosis
    GIdxLeft=max(1,TPosTmp-round(150*10^(-3)*samplerate));
    GIdxRight=min(TPosTmp+round(150*10^(-3)*samplerate),size(signal,1));
    t=(0:1/samplerate:(length(GIdxLeft:GIdxRight)-1)/samplerate)';
    fittedT=signal(GIdxLeft:GIdxRight).^2/trapz(signal(GIdxLeft:GIdxRight).^2);
    featureMatrix(1,1)=trapz(t,t.*fittedT); %center
    featureMatrix(2,1)=trapz(t,(t-featureMatrix(1,1)).^2.*fittedT); %variance
    featureMatrix(3,1)=trapz(t,(t-featureMatrix(1,1)).^3.*fittedT); %Schiefe/skewness
    featureMatrix(4,1)=trapz(t,(t-featureMatrix(1,1)).^4.*fittedT); %Wölbung/curtosis
    featureMatrix(18,1)=trapz(signal(GIdxLeft:GIdxRight).^2)/trapz(signal.^2); %AUC T wave
    
    
    % R-T-distance
    featureMatrix(5,1)=min(TzeroPos)-RPos;
    
    %R-T-distance with the middle of T wave
    if length(TPos)>1
        featureMatrix(6,1)=abs((TPos(2)-TPos(1))/2+TPos(1))-RPos;
    else
        leftM=find(abs(signal(TPos-round(150*10^(-3)*samplerate):TPos-1))<abs(0.5*TAmp(1)),1,'last')+TPos-round(150*10^(-3)*samplerate)-1;
        rightM=find(abs(signal(TPos+1:TPos+round(150*10^(-3)*samplerate)))<abs(0.5*TAmp(1)),1,'first')+TPos;
        featureMatrix(6,1)=(rightM-leftM)/2+leftM-RPos;
    end
    
    % Determine the two halfwaves
    if TminAmpPos<TmaxAmpPos
        TP=TminAmpPos;
        TP2=TmaxAmpPos;
        Vz=-1;
    else
        TP=TmaxAmpPos;
        TP2=TminAmpPos;
        Vz=1;
    end
    
    % Fitting two Gaussians as a surrogate for the peakedness of the first wave
    boundLeft=max(1,TP-round(150*10^(-3)*samplerate));
    boundRight=min(TP+round(150*10^(-3)*samplerate),size(signal,1));
    
    fitsig=signal(boundLeft:boundRight,1); % normalisation
    fitsig=fitsig-min(fitsig); % normalisation
    fitsig=fitsig/max(fitsig); % normalisation
    xfit=linspace(0,1,length(fitsig)); % normalised time vector
    
    pkPos=xfit(round(150*10^(-3)*samplerate)+1);
    if TP-round(150*10^(-3)*samplerate)>0
        npkPos=round(150*10^(-3)*samplerate)+1;
    else
        npkPos=TP;
    end
    lengthbw=length(signal(TP:TP2));
    
    % First Gaussian on the left part
    xexcl=ones(length(xfit),1);
    xexcl((round(npkPos-round(120*10^(-3)*samplerate)):round(npkPos+0.2*lengthbw))',1)=0;
    xexcl=find(xexcl);
    
    fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1,-2.001,0,-2],...
        'Upper',[1,2,1,2],...
        'StartPoint',[0.001,sign(Vz),pkPos,0],...
        'Algorithm','Trust-Region',...
        'Exclude',xexcl);
    ft = fittype('(k/sigma^2*exp(-(x-mu)^2/(2*sigma^2)))+h','independent','x','coefficients',{'sigma','k','mu','h'});
    curve1 = fit(xfit',fitsig,ft,fo);
    
    % Second Gaussian on the right part
    xexcl=ones(length(xfit),1);
    xexcl((round(npkPos-20*10^(-3)*samplerate):round(npkPos+0.6*lengthbw))',1)=0;
    xexcl=find(xexcl);
    fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1,-2.001,0,-2],...
        'Upper',[1,2,1,2],...
        'StartPoint',[0.001,sign(Vz),pkPos,0],...
        'Algorithm','Trust-Region',...
        'Exclude',xexcl);
    ft = fittype('(k/sigma^2*exp(-(x-mu)^2/(2*sigma^2)))+h','independent','x','coefficients',{'sigma','k','mu','h'});
    curve2 = fit(xfit',fitsig,ft,fo);
    
    % Fitting two Gaussians as a surrogate for the peakedness of the second wave
    boundLeft=max(1,TP2-round(150*10^(-3)*samplerate));
    boundRight=min(TP2+round(150*10^(-3)*samplerate),size(signal,1));
    
    fitsig2=signal(boundLeft:boundRight,1); % normalisation
    fitsig2=fitsig2-min(fitsig2); % normalisation
    fitsig2=fitsig2/max(fitsig2); % normalisation
    xfit=linspace(0,1,length(fitsig2)); % normalised time vector
    
    pkPos=xfit(round(150*10^(-3)*samplerate)+1);
    if TP2-round(150*10^(-3)*samplerate)>0
        npkPos=round(150*10^(-3)*samplerate)+1;
    else
        npkPos=TP2;
    end
    lengthbw=length(signal(TP:TP2));
    
    % First Gaussian on the right part
    xexcl=ones(length(xfit),1);
    xexcl((npkPos-round(0.2*lengthbw):round(npkPos+round(120*10^(-3)*samplerate)))',1)=0;
    xexcl=find(xexcl);
    fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1,-2.001,0,-2],...
        'Upper',[1,2,1,2],...
        'StartPoint',[0.001,sign(Vz),pkPos,0],...
        'Algorithm','Trust-Region',...
        'Exclude',xexcl);
    ft = fittype('(k/sigma^2*exp(-(x-mu)^2/(2*sigma^2)))+h','independent','x','coefficients',{'sigma','k','mu','h'});
    curve3 = fit(xfit',fitsig2,ft,fo);
    
    % Second Gaussian on left part
    xexcl=ones(length(xfit),1);
    xexcl((npkPos-round(0.6*lengthbw):round(npkPos+20*10^(-3)*samplerate))',1)=0;
    xexcl=find(xexcl);
    fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1,-2.001,0,-2],...
        'Upper',[1,2,1,2],...
        'StartPoint',[0.001,sign(Vz),pkPos,0],...
        'Algorithm','Trust-Region',...
        'Exclude',xexcl);
    ft = fittype('(k/sigma^2*exp(-(x-mu)^2/(2*sigma^2)))+h','independent','x','coefficients',{'sigma','k','mu','h'});
    curve4 = fit(xfit',fitsig2,ft,fo);
    
    
    % feature is the sum of both
    featureMatrix(7,1)=(curve1.sigma+curve2.sigma)*length(fitsig)/((curve3.sigma+curve4.sigma)*length(fitsig2));
    featureMatrix(16,1)=curve3.h-curve1.h;
    
    % T Amplitude
    featureMatrix(8,1)=abs(TmaxAmp)+abs(TminAmp);

    % Calculate upslope
    % cut a piece of the wave
    boundRight=TP+round(15*10^(-3)*samplerate);
    boundLeft=max(1,TP-round(80*10^(-3)*samplerate));
    xsl=(boundLeft:boundRight)';
    % fit a polynomial of 4th order
    ysl=signal(xsl);
    xsl=xsl-xsl(round(length(xsl)/2));
    Xsl=[xsl.^4,xsl.^3,xsl.^2,xsl,ones(length(xsl),1)];
    up4thGrade = lsqminnorm(Xsl,ysl);
    yslE=Xsl*up4thGrade;
    minNorm=min(yslE(1:length(yslE)-round(15*10^(-3)*samplerate)));
    maxNorm=max(yslE(1:length(yslE)-round(15*10^(-3)*samplerate)));
    yslEnorm=(yslE-minNorm)./(maxNorm-minNorm);
    yidx=find(yslEnorm<0.9 & yslEnorm>=0.3);
    yidx=yidx(yidx<length(yslE)-round(15*10^(-3)*samplerate));
    upslope=mean(diff(yslE(yidx,1)));
    
    
    % Calculate slope between Peaks
    boundRight=min(TP2+round(15*10^(-3)*samplerate),length(signal));
    boundLeft=max(1,TP-round(15*10^(-3)*samplerate));
    xsl=(boundLeft:boundRight)';
    % fit a polynomial of 4th order
    ysl=signal(xsl);
    xsl=xsl-xsl(round(length(xsl)/2));
    Xsl=[xsl.^4,xsl.^3,xsl.^2,xsl,ones(length(xsl),1)];
    up4thGrade = lsqminnorm(Xsl,ysl);
    ysmE=Xsl*up4thGrade;
    mEsample=(1:1:length(ysmE))';
    [~,posmx]=max(ysmE);
    [~,posmn]=min(ysmE);
    if posmn<posmx
        yidx=find(ysmE<0.7*max(ysmE) & ysmE>=(max(ysmE)-min(ysmE))*0.2+min(ysmE) & mEsample<posmx & mEsample>posmn);
    else
        yidx=find(ysmE<0.7*max(ysmE) & ysmE>=(max(ysmE)-min(ysmE))*0.2+min(ysmE) & mEsample>posmx & mEsample<posmn);
    end
    yidx=yidx(yidx<length(ysmE)-16);
    midslope=mean(diff(ysmE(yidx,1)));
    featureMatrix(9,1)=upslope/midslope;
    
    % Calculate downslope
    % cut a piece of the wave
    boundRight=min(TP2+round(80*10^(-3)*samplerate),size(signal,1));
    boundLeft=max(1,TP2-round(15*10^(-3)*samplerate));
    % fit a polynomial of 4th order
    xsr=(boundLeft:boundRight)';
    ysr=signal(xsr);
    xsr=xsr-xsr(round(length(xsr)/2));
    Xsr=[xsr.^4,xsr.^3,xsr.^2,xsr,ones(length(xsr),1)];
    down4thGrade = lsqminnorm(Xsr,ysr);
    ysrE=Xsr*down4thGrade;
    minNorm=min(ysrE(round(15*10^(-3)*samplerate):end));
    maxNorm=max(ysrE(round(15*10^(-3)*samplerate):end));
    ysrEnorm=(ysrE-minNorm)./(maxNorm-minNorm);
    yidx=find(ysrEnorm<0.9 & ysrEnorm>=0.3);
    yidx=yidx(yidx>round(15*10^(-3)*samplerate));
    downslope=mean(diff(ysrE(yidx,1)));
    featureMatrix(10,1)=downslope/midslope;
    
    % Fraction of left/right half at the area of the T wave
    leftIdx=max(1,min(TzeroPos)-round(150*10^(-3)*samplerate));
    rightIdx=min(min(TzeroPos)+round(150*10^(-3)*samplerate),size(signal,1));
    featureMatrix(11,1)=sum(signal(leftIdx:min(TzeroPos)-1))/sum(leftIdx:rightIdx);
    featureMatrix(12,1)=sum(signal(min(TzeroPos)+1:rightIdx))/sum(leftIdx:rightIdx);
    
    % QRS features: Amplitude, Energy, Ratio
    featureMatrix(13,1)=signal(RPos);
    RIdxLeft=max(1,RPos-0.035*samplerate);
    RIdxRight=min(RPos+0.035*samplerate,size(signal,1));
    featureMatrix(14,1)=sum(signal(RIdxLeft:RIdxRight).^2);
    featureMatrix(15,1)=featureMatrix(14,1)/featureMatrix(13,1);
    
    featureMatrix(17,1)=1;
end
end




