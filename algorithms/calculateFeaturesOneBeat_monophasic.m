function [featureMatrix,featureNames]=calculateFeaturesOneBeat_monophasic(ecg_matrix,samplerate,TPos,RPos,TAmp)
featureNames={'center'; 'variance'; 'skewness'; 'curtosis'; 'RT distance'; 'RTmid distance';...
    'Peakness T wave'; 'T amplitude'; 'T upslope'; 'T downslope'; 'left half T wave area ratio';...
    'right half T wave area ratio'; 'R Amplitude';'R power';'ratio R power/R amplitude';'ST elevation';'Biphasic';'AUC T Wave/AUC ECG Beat'};
%%
signal=ecg_matrix;
featureMatrix=nan(size(featureNames,1),1);

%% Calculate statistical features: Center, Variance, Skewness and Curtosis
GIdxLeft=max(1,TPos-150*10^(-3)*samplerate);
GIdxRight=min(TPos+150*10^(-3)*samplerate,size(signal,1));
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
boundLeft=max(1,TPos(1)-150*10^(-3)*samplerate);
boundRight=min(TPos(end)+150*10^(-3)*samplerate,size(signal,1));

fitsig=signal(boundLeft:boundRight,1); % normalisation
fitsig=fitsig-min(fitsig); % normalisation
fitsig=fitsig/max(fitsig); % normalisation
xfit=linspace(0,1,length(fitsig)); % normalised time vector

pkPos=xfit(150*10^(-3)*samplerate+1);
if TPos(1)-150*10^(-3)*samplerate>0
    npkPos=round(150*10^(-3)*samplerate)+1;
else
    npkPos=TPos(1);
end
    

% First Gaussian on the descending part (right part)
xexcl=ones(length(xfit),1);
xexcl((round(npkPos-20*10^(-3)*samplerate):round(npkPos+120*10^(-3)*samplerate))',1)=0;
xexcl=find(xexcl);
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[-1,-2.001,0,-2],...
    'Upper',[1,2,1,2],...
    'StartPoint',[0.001,sign(signal(TPos)),pkPos,0],...
    'Algorithm','Trust-Region',...
    'Exclude',xexcl);
ft = fittype('(k/sigma^2*exp(-(x-mu)^2/(2*sigma^2)))+h','independent','x','coefficients',{'sigma','k','mu','h'});
curve1 = fit(xfit',fitsig,ft,fo);
%figure; plot(curve1); hold on; plot(xfit',fitsig)

% Second Gaussian on the ascending part (left part)
xexcl=ones(length(xfit),1);
xexcl((round(npkPos-120*10^(-3)*samplerate):round(npkPos+20*10^(-3)*samplerate))',1)=0;
xexcl=find(xexcl);
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[-1,-2.001,0,-2],...
    'Upper',[1,2,1,2],...
    'StartPoint',[0.001,sign(signal(TPos)),pkPos,0],...
    'Algorithm','Trust-Region',...
    'Exclude',xexcl);
ft = fittype('(k/sigma^2*exp(-(x-mu)^2/(2*sigma^2)))+h','independent','x','coefficients',{'sigma','k','mu','h'});
curve2 = fit(xfit',fitsig,ft,fo);
%figure; plot(curve2); hold on; plot(xfit',fitsig)


% feature is the sum of both
featureMatrix(7,1)=curve1.sigma*length(fitsig)+curve2.sigma*length(fitsig);
%figure(sLead);
featureMatrix(16,1)=curve2.h-curve1.h;
%plot(xfit,fitsig, 'LineWidth', 2);hold on; plot(xfit,curve1(xfit)); %hold on; plot(xfit,curve2(xfit))


%% T Amplitude
featureMatrix(8,1)=TAmp;

%% Calculate upslope
% cut a piece of the wave
boundRight=min(length(signal),TPos(1)+round(15*10^(-3)*samplerate));
boundLeft=max(1,TPos(1)-round(80*10^(-3)*samplerate));
xsl=(boundLeft:boundRight)';
% fit a polynomial of 4th order
ysl=signal(xsl);
xsl=xsl-xsl(round(length(xsl)/2));
Xsl=[xsl.^4,xsl.^3,xsl.^2,xsl,ones(length(xsl),1)];
%down4thGrade=(transpose(X2)*X2)\transpose(X2)*y2;
up4thGrade = lsqminnorm(Xsl,ysl);
yslE=Xsl*up4thGrade;
minNorm=min(yslE(1:length(yslE)-round(15*10^(-3)*samplerate)));
maxNorm=max(yslE(1:length(yslE)-round(15*10^(-3)*samplerate)));
yslEnorm=(yslE-minNorm)./(maxNorm-minNorm);
yidx=find(yslEnorm<0.9 & yslEnorm>=0.3);
yidx=yidx(yidx<length(yslE)-round(15*10^(-3)*samplerate));
featureMatrix(9,1)=mean(diff(yslE(yidx,1)));

%% Calculate downslope
% cut a piece of the wave
boundRight=min(TPos(end)+80*10^(-3)*samplerate,size(signal,1));
boundLeft=max(1,TPos(end)-15);

xsr=(boundLeft:boundRight)';
ysr=signal(xsr);
xsr=xsr-xsr(round(length(xsr)/2));
Xsr=[xsr.^4,xsr.^3,xsr.^2,xsr,ones(length(xsr),1)];
%down4thGrade=(transpose(X2)*X2)\transpose(X2)*y2;
down4thGrade = lsqminnorm(Xsr,ysr);
ysrE=Xsr*down4thGrade;
minNorm=min(ysrE(round(15*10^(-3)*samplerate):end));
maxNorm=max(ysrE(round(15*10^(-3)*samplerate):end));
ysrEnorm=(ysrE-minNorm)./(maxNorm-minNorm);
yidx=find(ysrEnorm<0.9 & ysrEnorm>=0.3);
yidx=yidx(yidx>round(15*10^(-3)*samplerate));
featureMatrix(10,1)=mean(diff(ysrE(yidx,1)));


%% Fraction of left/right half at the area of the T wave
leftIdx=max(1,TPos-150*10^(-3)*samplerate);
rightIdx=min(TPos+150*10^(-3)*samplerate,size(signal,1));
featureMatrix(11,1)=sum(signal(leftIdx:TPos-1))/sum(signal(leftIdx:rightIdx));
featureMatrix(12,1)=sum(signal(TPos+1:rightIdx))/sum(signal(leftIdx:rightIdx));

%% QRS features: Amplitude, Energy, Ratio
featureMatrix(13,1)=signal(RPos);
RIdxLeft=max(1,RPos-0.035*samplerate);
RIdxRight=min(RPos+0.035*samplerate,size(signal,1));
%featureMatrix(14,1)=sum(signal(RIdxLeft:RIdxRight).^2);
featureMatrix(14,1)=sum(signal(RIdxLeft:RIdxRight).^2);
%     plot(1:QRSfinish,signal(1:QRSfinish),'o');
featureMatrix(15,1)=featureMatrix(14,1)/featureMatrix(13,1);

%% Biphasic parameter is 0
featureMatrix(17,1)=0;
end
