function [featureMatrix,featureNames]=calculateFeaturesOneBeat_biphasic(ecg_matrix,samplerate,TPos,RPos,TAmp)
featureNames={'center'; 'variance'; 'skewness'; 'curtosis'; 'RT distance'; 'RTmid distance';...
    'Peakness T wave'; 'T amplitude'; 'T upslope'; 'T downslope'; 'left half T wave area ratio';...
    'right half T wave area ratio'; 'R Amplitude';'R power';'ratio R power/R amplitude';'ST elevation';'Biphasic';'AUC T Wave/AUC ECG Beat'};

%%
featureMatrix=nan(18,1);
signal=ecg_matrix(:,1);

%estimate end of QRS complex
QRSend=0.2*samplerate;

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
%down4thGrade=(transpose(X2)*X2)\transpose(X2)*y2;
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
        GIdxLeft=max(1,TPos-150*10^(-3)*samplerate);
        GIdxRight=min(TPos+150*10^(-3)*samplerate,size(signal,1));
        leftM=find(abs(signal(GIdxLeft:TPos-1))<abs(0.5*TAmp(1)),1,'last')+TPos-150*10^(-3)*samplerate-1;
        rightM=find(abs(signal(TPos+1:GIdxRight))<abs(0.5*TAmp(1)),1,'first')+TPos;
        
        if ~(isempty(leftM)||isempty(rightM))
            TPosTmp=round((rightM-leftM)/2+leftM);
        else
            TPosTmp=TPos(1);
        end
    end
    %     GIdxLeft=max(1,TPosTmp-150*10^(-3)*samplerate);
    %     GIdxRight=min(TPosTmp+150*10^(-3)*samplerate,size(signal,1));
    %     x=(TPosTmp-GIdxLeft:TPosTmp+GIdxRight)';
    %     x=x(x<=length(ecg_matrix));
    %     x=x(x>0);
    %     fittedT=signal(x).^2/sum(signal(x).^2);
    %     featureMatrix(1,1)=sum((0:1:length(fittedT)-1)'.*fittedT); %center
    %     featureMatrix(2,1)=sum(((0:1:length(fittedT)-1)'-featureMatrix(1,1)).^2.*fittedT); %variance
    %     featureMatrix(3,1)=sum(((0:1:length(fittedT)-1)'-featureMatrix(1,1)).^3.*fittedT); %Schiefe/skewness
    %     featureMatrix(4,1)=sum(((0:1:length(fittedT)-1)'-featureMatrix(1,1)).^4.*fittedT); %Wlbung/curtosis
    
    
    %% Calculate statistical features: Center, Variance, Skewness and Curtosis
    GIdxLeft=max(1,TPosTmp-150*10^(-3)*samplerate);
    GIdxRight=min(TPosTmp+150*10^(-3)*samplerate,size(signal,1));
    t=(0:1/samplerate:(length(GIdxLeft:GIdxRight)-1)/samplerate)';
    fittedT=signal(GIdxLeft:GIdxRight).^2/trapz(signal(GIdxLeft:GIdxRight).^2);
    featureMatrix(1,1)=trapz(t,t.*fittedT); %center
    featureMatrix(2,1)=trapz(t,(t-featureMatrix(1,1)).^2.*fittedT); %variance
    featureMatrix(3,1)=trapz(t,(t-featureMatrix(1,1)).^3.*fittedT); %Schiefe/skewness
    featureMatrix(4,1)=trapz(t,(t-featureMatrix(1,1)).^4.*fittedT); %Wölbung/curtosis
    featureMatrix(18,1)=trapz(signal(GIdxLeft:GIdxRight).^2)/trapz(signal.^2); %AUC T wave
    
    
    %% R-T-distance
    %R-T-distance
    featureMatrix(5,1)=min(TzeroPos)-RPos;
    
    %R-T-distance with the middle of T wave
    if length(TPos)>1
        featureMatrix(6,1)=abs((TPos(2)-TPos(1))/2+TPos(1))-RPos;
    else
        leftM=find(abs(signal(TPos-150*10^(-3)*samplerate:TPos-1))<abs(0.5*TAmp(1)),1,'last')+TPos-150*10^(-3)*samplerate-1;
        rightM=find(abs(signal(TPos+1:TPos+150*10^(-3)*samplerate))<abs(0.5*TAmp(1)),1,'first')+TPos;
        featureMatrix(6,1)=(rightM-leftM)/2+leftM-RPos;
    end
    
    
    %% Fitting of a 2nd order polynomial for the peakness
    if TminAmpPos<TmaxAmpPos
        TP=TminAmpPos;
        TP2=TmaxAmpPos;
        Vz=-1;
    else
        TP=TmaxAmpPos;
        TP2=TminAmpPos;
        Vz=1;
    end
    
    %% Fitting two Gaussians as a surrogate for the peakedness of the first wave
    boundLeft=max(1,TP-150*10^(-3)*samplerate);
    boundRight=min(TP+150*10^(-3)*samplerate,size(signal,1));
    
    fitsig=signal(boundLeft:boundRight,1); % normalisation
    fitsig=fitsig-min(fitsig); % normalisation
    fitsig=fitsig/max(fitsig); % normalisation
    xfit=linspace(0,1,length(fitsig)); % normalised time vector
    
    pkPos=xfit(150*10^(-3)*samplerate+1);
    if TP-150*10^(-3)*samplerate>0
        npkPos=round(150*10^(-3)*samplerate)+1;
    else
        npkPos=TP;
    end
    lengthbw=length(signal(TP:TP2));
    
    
    % First Gaussian on the left part
    xexcl=ones(length(xfit),1);
    xexcl((round(npkPos-120*10^(-3)*samplerate):round(npkPos+0.2*lengthbw))',1)=0;
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
    
    %% Fitting two Gaussians as a surrogate for the peakedness of the second wave
    boundLeft=max(1,TP2-150*10^(-3)*samplerate);
    boundRight=min(TP2+150*10^(-3)*samplerate,size(signal,1));
    
    fitsig2=signal(boundLeft:boundRight,1); % normalisation
    fitsig2=fitsig2-min(fitsig2); % normalisation
    fitsig2=fitsig2/max(fitsig2); % normalisation
    xfit=linspace(0,1,length(fitsig2)); % normalised time vector
    
    pkPos=xfit(150*10^(-3)*samplerate+1);
    if TP2-150*10^(-3)*samplerate>0
        npkPos=round(150*10^(-3)*samplerate)+1;
    else
        npkPos=TP2;
    end
    lengthbw=length(signal(TP:TP2));
    
    % First Gaussian on the right part
    xexcl=ones(length(xfit),1);
    xexcl((npkPos-round(0.2*lengthbw):round(npkPos+120*10^(-3)*samplerate))',1)=0;
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
    
    %% T Amplitude
    featureMatrix(8,1)=abs(TmaxAmp)+abs(TminAmp);
    %     figure(4*sLead);
    %     plot(signal,'LineWidth',2);
    %     hold on;
    %     plot(TmaxAmpPos,TmaxAmp,'o');
    %     plot(TminAmpPos, TminAmp,'o');
    
    
    
    %% Calculate upslope
    % cut a piece of the wave
    boundRight=TP+round(15*10^(-3)*samplerate);
    boundLeft=max(1,TP-round(80*10^(-3)*samplerate));
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
    %featureMatrix(10,setupLead)=mean(diff(yslE(yidx,1)));
    upslope=mean(diff(yslE(yidx,1)));
    %   plot(yidx+boundLeft,yslE(yidx),'LineWidth',2)
    
    
    %% Calculate slope between Peaks
    boundRight=min(TP2+round(15*10^(-3)*samplerate),length(signal));
    boundLeft=max(1,TP-round(15*10^(-3)*samplerate));
    xsl=(boundLeft:boundRight)';
    % fit a polynomial of 4th order
    ysl=signal(xsl);
    xsl=xsl-xsl(round(length(xsl)/2));
    Xsl=[xsl.^4,xsl.^3,xsl.^2,xsl,ones(length(xsl),1)];
    %down4thGrade=(transpose(X2)*X2)\transpose(X2)*y2;
    up4thGrade = lsqminnorm(Xsl,ysl);
    ysmE=Xsl*up4thGrade;
    mEsample=(1:1:length(ysmE))';
    %yidx=find(ysmE<0.9*ysmE(end-16) & ysmE>=0.2*ysmE(end-16));
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
    %featureMatrix(11,setupLead)=mean(diff(ysmE(yidx,1)));
    %       plot(yidx+boundLeft,ysmE(yidx),'LineWidth',2)
    
    
    %% Calculate downslope
    % cut a piece of the wave
    boundRight=min(TP2+round(80*10^(-3)*samplerate),size(signal,1));
    boundLeft=max(1,TP2-round(15*10^(-3)*samplerate));
    % fit a polynomial of 4th order
    
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
    downslope=mean(diff(ysrE(yidx,1)));
    %featureMatrix(12,setupLead)=mean(diff(ysrE(yidx,1)));
    featureMatrix(10,1)=downslope/midslope;
    %       plot(yidx+boundLeft,ysrE(yidx),'LineWidth',2)
    
    
    %% Fraction of left/right half at the area of the T wave
    leftIdx=max(1,min(TzeroPos)-150*10^(-3)*samplerate);
    rightIdx=min(min(TzeroPos)+150*10^(-3)*samplerate,size(signal,1));
    featureMatrix(11,1)=sum(signal(leftIdx:min(TzeroPos)-1))/sum(leftIdx:rightIdx);
    featureMatrix(12,1)=sum(signal(min(TzeroPos)+1:rightIdx))/sum(leftIdx:rightIdx);
    %     figure(11);
    %     plot(signal);
    %     hold on;
    %     plot(leftIdx:min(TzeroPos)-1,signal(leftIdx:min(TzeroPos)-1),'o');
    
    %% QRS features: Amplitude, Energy, Ratio
    
    %% QRS features: Amplitude, Energy, Ratio
    featureMatrix(13,1)=signal(RPos);
    RIdxLeft=max(1,RPos-0.035*samplerate);
    RIdxRight=min(RPos+0.035*samplerate,size(signal,1));
    %featureMatrix(14,setupLead)=sum(signal(RIdxLeft:RIdxRight).^2);
    featureMatrix(14,1)=sum(signal(RIdxLeft:RIdxRight).^2);
    %     plot(1:QRSfinish,signal(1:QRSfinish),'o');
    featureMatrix(15,1)=featureMatrix(14,1)/featureMatrix(13,1);
    
    featureMatrix(17,1)=1;
end
end
