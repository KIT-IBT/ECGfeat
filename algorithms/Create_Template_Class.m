function [Template,posRpeak,ampRpeak,Template_matrix,FPTtemplate]=Create_Template_Class(signal,samplerate,FPT,Class,varargin)
%Function used to generate a Template for a specified type of beat. A
%previous classification of the beats is needed.
%The input variables are as follows:
%signal: ECG signal
%FPT: Fidutial Point Table containing position of the R-peaks an their
%corresponding classes.
%Class: The class for which the Template should be generated.
%varargin: The variable can be used to select a template of a whole ECG or
%just the QRS complex. Use 'ECG' for a whole ECG of 'QRS' for the
%QRS-complex. Additionally, you can give a desired length for the templates
%as vector [L1,L2], where L1 is the length before the R peak and L2 is the
%length after the R peak
%The output variables are as follows:
%Template: Signal contaning the calculated Template
%posRpeak: Position of the R-peak in the Template.
%ampRpeak: origianl amplitude of the R-peak.
%posPTpeaks: Position of P/T peak in the template

%Check if the given signal is a vector
if ~isvector(signal)
    error('The input ECG signal must be a vector.');
end

%Check if QRS complex or full ECG is required
if isempty(varargin)
    varargin={'ECG'};
    L1=round(min(1/3*quantile(diff(FPT(:,6)),0.8),500/1000*samplerate));
    L2=round(min(2/3*quantile(diff(FPT(:,6)),0.8),1000/1000*samplerate));
    flag_length=0;
else
    if length(varargin)>1
        L1=varargin{2}(1,1);
        L2=varargin{2}(1,2);
        flag_length=1;
    else
    L1=round(min(1/3*quantile(diff(FPT(:,6)),0.8),500/1000*samplerate));
    L2=round(min(2/3*quantile(diff(FPT(:,6)),0.8),1000/1000*samplerate));
        flag_length=0;
    end
end
L=L1+L2+1;

if strcmp(varargin{1},'ECG') || strcmp(varargin{1},'ecg') || strcmp(varargin{1},'QRS') || strcmp(varargin{1},'qrs')
    case_var=0;
    if strcmp(varargin{1},'QRS') || strcmp(varargin{1},'qrs')
        case_var=1;
    end
else
    error('Template type not recognized');
end

if size(signal,2)>1
    if size(signal,1)==1
        signal=signal';
    elseif size(signal,1)<size(signal,2)
        signal=signal(1,:)';
    end
elseif size(signal,2)>size(signal,1)
    signal=signal(:,1);
end

if ~(Class==1 || Class==2 || Class==3)
    error('Class type nor recognized. Only the values 1,2 and 3 are allowed.');
end

%index in the FPT of the beats classified as Class
class_index=find(FPT(:,13)==Class);
if isempty(class_index)
    display('Warning: No beats of the given class are present in the signal!');
    Template=0;
    posRpeak=0;
    ampRpeak=0;
    return
end

FPTtemplate=nan(1,13);

Preselected_Beats_Index=class_index;


%For the normal beats (class 1) a preselection is applied. Beats having
%strong deviations in their rhythmical properties are not considered. The
%deviations are measured using the Poincare plot. The procedure is repeated
%twice to ensure robustness.
if Class==1
    RR=diff(FPT(:,6));
    X=[RR(1:end-1),RR(2:end)];
    index=1:1:size(X,1);
    for i=1:2
        SCORE=(X-[mean(X(index,1))*ones(size(X,1),1),mean(X(index,2))*ones(size(X,1),1)])*1/sqrt(2)*[1,-1;1,1];
        D1=abs(SCORE(:,1));
        D2=abs(SCORE(:,2));
        Thl1=2.5*std(D1);
        Thl2=0.7*std(D2);
        index=D1<Thl1 & D2<Thl2;
        if all(~index)
            break
        end
    end
    Preselected_Beats_Index=find(index)+1;
    
    %In case too little beats are selected to build a template a new
    %approach is taken
    if length(Preselected_Beats_Index)<1/3*size(FPT,1)
        if all(~index)
            index=(SCORE(:,1)>=-Thl1 & SCORE(:,2)<=Thl2);
        else
            index=(SCORE(:,1)>=min(SCORE(index,1)) & SCORE(:,2)<=0);
        end
        Preselected_Beats_Index=find(index)+1;
        if ~flag_length
            L1=round(min(1/3*quantile(diff(FPT(:,6)),0.8),500/1000*samplerate));
            L2=round(min(2/3*quantile(X(index,2),0.8),1000/1000*samplerate));
            L=L1+L2+1;
        end
    end
end

%If the beats at the beginning and ending borders of ECG signal have
%smaller sizes than the template, they are removed
if FPT(Preselected_Beats_Index(1),6)-L1<1
    ind=find(FPT(Preselected_Beats_Index,6)-L1>=1,1,'first');
    Preselected_Beats_Index=Preselected_Beats_Index(ind:end);
end
if FPT(Preselected_Beats_Index(end),6)+L2>length(signal)
    ind=find(FPT(Preselected_Beats_Index,6)+L2<=length(signal),1,'last');
    Preselected_Beats_Index=Preselected_Beats_Index(1:ind);
end


%Beats having strong deviations in their morphology (maximal and minimal
%amplitudes) are not considered for now
Matrix_ecg=zeros(L,length(Preselected_Beats_Index));
MP=zeros(length(Preselected_Beats_Index),2);
for i=1:length(Preselected_Beats_Index)
    Matrix_ecg(:,i)=signal(FPT(Preselected_Beats_Index(i),6)-L1:FPT(Preselected_Beats_Index(i),6)+L2,1);
    MP(i,:)=[max(Matrix_ecg(:,i)),min(Matrix_ecg(:,i))];
end
Th11=quantile(MP(:,1),0.25);
Th12=quantile(MP(:,1),0.75);
Th21=quantile(MP(:,2),0.25);
Th22=quantile(MP(:,2),0.75);
Selected_Beats_Index=MP(:,1)>=Th11 & MP(:,1)<=Th12 & MP(:,2)>=Th21 & MP(:,2)<=Th22;
if ~sum(Selected_Beats_Index)
    Selected_Beats_Index=true(size(Preselected_Beats_Index));
end
Template=mean(Matrix_ecg(:,Selected_Beats_Index),2);

%Deviation in morphology is measured using the correlation coefficient.
%Beats having a maximum cross covariance lower than 0.75 are not considered.
%The procedure is repeated to ensure robustness
nit=3;
for j=1:nit
    T2=zeros(size(Matrix_ecg,2),1);
    for i=1:length(T2)
        T2(i)=l_operator(Matrix_ecg(:,i),Template);
    end
    Th=quantile(T2,0.15);
    
    if any(T2>=Th)
        Template=mean(Matrix_ecg(:,T2>=Th),2);
    else
        Template=mean(Matrix_ecg(:,Preselected_Beats_Index),2);
    end
    
    if j==nit
        Th=max(quantile(T2,0.2),0.90); %before gl121: 0.75
        if any(T2>=Th)
            Template=mean(Matrix_ecg(:,T2>=Th),2);
        end
    end
end

%Other important properties of the Template
Template_matrix=Matrix_ecg(:,T2>=Th);

if isempty(Template_matrix)
    disp('Threshold is too strict. Mean template is used...')
    Matrix_ecg=zeros(L,length(class_index));
    for i=1:length(class_index)
        Matrix_ecg(:,i)=signal(FPT(class_index(i),6)-L1:FPT(class_index(i),6)+L2,1);
    end
    Template=mean(Matrix_ecg(:,class_index),2);
    Template_matrix=Matrix_ecg(:,class_index);
end
%Offset removal
[~,offset]=Isoline_Correction(Template_matrix(:));
Template=Template-offset;
Template_matrix=Template_matrix-offset;

%Position of R peak
ampRpeak=Template(L1+1);
posRpeak=L1+1;


%%
%The QRS complex is segmented in case the user wants only the QRS complex.
if case_var
    n1=round(samplerate*0.08);
    n2=round(samplerate*0.11);
    if L1-n1<1
        n1=L1;
    end
    if n2>L2
        n2=L2-1;
    end
    Template=Template(L1+1-n1:L1+1+n2);
    Template_matrix=Template_matrix(L1+1-n1:L1+1+n2,:);
    ampRpeak=Template(n1+1);
    posRpeak=n1+1;
    % if the whole beat is taken, get P and T peak!
else
    FPTtemplate(1,6)=posRpeak;
    % Get T peak and P peak
    estRTdist=round(mean(FPT(Preselected_Beats_Index,11)-FPT(Preselected_Beats_Index,6)));
    estPRdist=round(mean(FPT(Preselected_Beats_Index,6)-FPT(Preselected_Beats_Index,2)));
    dTemplate=diff(abs(Template));
    maxPosAbsTemplate=find(dTemplate(1:end-1,1)>=0 & dTemplate(2:end,1)<0)+1;
    Pcand=maxPosAbsTemplate(maxPosAbsTemplate<(posRpeak-estPRdist*0.7) & maxPosAbsTemplate>(posRpeak-estPRdist*1.3));
    [~,PpeakPos]=max(abs(Template(Pcand)));
    FPTtemplate(1,2)=Pcand(PpeakPos);
    Tcand=maxPosAbsTemplate(maxPosAbsTemplate>(posRpeak+estRTdist*0.7) & maxPosAbsTemplate<(posRpeak+estRTdist*1.3));
    [~,TpeakPos]=max(abs(Template(Tcand)));
    FPTtemplate(1,11)=Tcand(TpeakPos);
end

if ampRpeak~=0
    Template=Template/ampRpeak;
    Template_matrix=Template_matrix/ampRpeak;
else
    [~,posRpeak]=max(abs(Template(:,1)));
    ampRpeak=Template(posRpeak);
    Template=Template/ampRpeak;
    Template_matrix=Template_matrix/ampRpeak;
end
