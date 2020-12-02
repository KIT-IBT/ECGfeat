function [FPT,FPTcell] = getFPTFromSimulations(signal,samplerate)

%Get the R peak

FPTcell=cell(size(signal,2),1);
for ld=1:1:size(signal,2)
    FPTsingle=nan(1,13);
    [pksRpos,locsRpos]=findpeaks(signal(:,ld));
    [pksRneg,locsRneg]=findpeaks(-signal(:,ld));
    [locsR,ordR]=sort([locsRpos;locsRneg]);
    pksR=[pksRpos;pksRneg];
    pksR=pksR(ordR);
    [~,tmp]=max(pksR);
    FPTsingle(1,6)=locsR(tmp,1);
    FPTcell{ld,1}=FPTsingle;
    
end
FPT=cell2mat(FPTcell);

% Get the end of QRS

if size(signal,2)>1
    [~,orderS]=sort(std(signal,1,2));
    pos=find(orderS>max(FPT(:,6)) & orderS<200,1,'first');
    QRSend=orderS(pos,1);
else
    QRSend=max(FPT(:,6))+0.1*samplerate;
end

% Get T peak

for ld=1:1:size(signal,2)
    [pksRpos,locsRpos]=findpeaks(signal(:,ld));
    [pksRneg,locsRneg]=findpeaks(-signal(:,ld));
    [locsR,ordR]=sort([locsRpos;locsRneg]);
    pksR=[pksRpos;pksRneg];
    pksR=pksR(ordR);
    locsT=locsR(locsR>QRSend);
    pksT=pksR(locsR>QRSend);
    [~,tmp]=max(pksT);
    FPTcell{ld,1}(1,11)=locsT(tmp,1); 
end

FPT=round(median(cell2mat(FPTcell),1));

% Get the end of T

if size(signal,2)>1
    [~,orderS]=sort(std(signal,1,2));
    pos=find(orderS>max(FPT(:,11)),1,'first');
    Tend=orderS(pos,1);
else
    Tend=max(FPT(:,11))+0.15*samplerate;
end

FPT(1,7)=QRSend;
FPT(1,12)=Tend;

% set non existent P wave to 1

for ld=1:1:size(signal,2)
    FPTcell{ld,1}(1,2)=1; 
end
FPT(1,2)=1;

end

