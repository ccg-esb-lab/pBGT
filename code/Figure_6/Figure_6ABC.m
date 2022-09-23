clc
clear all
close all

run('../matlab_code/lib/addpath_recurse');
addpath_recurse('../matlab_code/lib/');

%% LOAD DATA


dataPath='../../data/LB-AMP-LB/';


%***** Plasmid
pfileName1='pBGT_LB-AMP-LB_AMNIS_reread.csv';
pfileName='pBGT_LB-AMP-LB_AMNIS.csv';
cfileName='MGGT_LB-AMP-LB_AMNIS.csv';

numDays=3;
nreps=3;

colsP=[2,3,4,5,6,7,8,9,10,11,12];
AsP=[0,1/256, 1/128, 1/64, 1/32, 1/16,1/8,1/4,1/2,1];
strAsP={'AMP=0','1/256 MIC', '1/128 MIC','1/64 MIC', '1/32 MIC','1/16 MIC','1/8 MIC','1/4 MIC','1/2 MIC','MIC'};

colsC=[2,3,4,5,6,7,8,9,10,11];
AsC=[0,1/128, 1/64, 1/32,1/16,1/8,1/4,1/2,1];
strAsC={'AMP=0','1/128 MIC','1/64 MIC','1/32 MIC', '1/16 MIC','1/8 MIC','1/4 MIC','1/2 MIC','MIC'};



%% PLOT plasmids DISTRIBUTIONS

lim_aspect=[.2, 1];
lim_area=[25, 35];
xmax=.75e6;
xmin=.9e4;
xdelta=.5e4;
col=[0.3020    0.6863    0.2902];
meanGFPsP=zeros(numDays, length(AsP), nreps);

figure(); clf('reset');set(gcf,'DefaultLineLineWidth',1); set(gcf, 'color', 'white'); hold all

st={'-','-',':'};
stcol={'k',col,'k'};
e=xmin:xdelta:xmax;
iA=9;  %Plasmid

for day=1:numDays
    
    if(day==1)
        fileName=pfileName1;
    else
        fileName=pfileName;
    end
    
       
    filePath=[dataPath,fileName];
    data = readtable(filePath,'Delimiter', '\t');

    dayData=data(data.season==day & data.doseNumber==colsP(iA),:);
    %disp(head(dayData))
    repH=[];
    meanGFP=zeros(1,nreps);
    for irep=1:nreps
        
        repname=string(['rep-',sprintf('%d',irep)]);
        asrepData=dayData(dayData.replicate==repname,:);
        thisGFP=asrepData.Intensity_MC_Ch02(asrepData.Area_M02 > lim_area(1) & asrepData.Area_M02 < lim_area(2) & asrepData.AspectRatio_M02>lim_aspect(1) & asrepData.AspectRatio_M02<lim_aspect(2),:);
        
        h=hist(thisGFP(thisGFP>xmin & thisGFP<xmax), e);
        
        hnorm=h./max(h);
        %semilogx(e,hnorm,'LineWidth',1,'Color',[.8 .8 .8]); hold on;
        
        meanGFP(irep)=mean(thisGFP);
        repH=[repH; hnorm];
    
    end
    semilogx(e,mean(repH),st{day},'LineWidth',3,'Color',stcol{day}); hold on;
    axis([xmin xmax*.4 0 1]);
    
    ylabel('Frequency','FontSize',20);
    xlabel('GFP intensity (log scale)','FontSize',20);
end

set(gca,'fontsize',18);
legend({'Season 1: LB', 'Season 2: LB+AMP','Season 3: LB'},'FontSize',18);
legend boxoff

export_fig '../../figures/Fig-6A.pdf';

%% PLOT chromosome DISTRIBUTIONS

lim_aspect=[.2, 1];
lim_area=[1, 35];
xmax=7.5e3;
xmin=5e2;
xdelta=.5e2;
col=[0.2157    0.4941    0.7216];

figure(); clf('reset');set(gcf,'DefaultLineLineWidth',1); set(gcf, 'color', 'white'); hold all

st={'-','-',':'};
stcol={'k',col,'k'};

iA=9;  %Plasmid

for day=1:numDays
    
    fileName=cfileName;
          
    filePath=[dataPath,fileName];
    data = readtable(filePath,'Delimiter', '\t');

    dayData=data(data.season==day & data.doseNumber==colsC(iA),:);
    %disp(head(dayData));
    repH=[];
    meanGFP=zeros(1,nreps);
    for irep=1:nreps
        
        repname=string(['rep-',sprintf('%d',irep)]);
        asrepData=dayData(dayData.replicate==repname,:);
        %disp(head(asrepData));
        thisGFP=asrepData.Intensity_MC_Ch02(asrepData.Area_M02 > lim_area(1) & asrepData.Area_M02 < lim_area(2) & asrepData.AspectRatio_M02>lim_aspect(1) & asrepData.AspectRatio_M02<lim_aspect(2),:);
        L=length(thisGFP);
        if(L>1000)
            e=xmin:xdelta:xmax;
            h=hist(thisGFP(thisGFP>xmin & thisGFP<xmax), e);
            
            hnorm=h./max(h);
            %semilogx(e,hnorm,'LineWidth',1,'Color',[.8 .8 .8]); hold on;
            
            meanGFP(irep)=mean(thisGFP);
            repH=[repH; hnorm];
        end
    end
    
    semilogx(e,mean(repH),st{day},'LineWidth',3,'Color',stcol{day}); hold on;
    axis([xmin xmax*.4 0 1]);
    
    ylabel('Frequency','FontSize',20);
    xlabel('GFP intensity (log scale)','FontSize',20);
end

set(gca,'fontsize',18);
legend({'Season 1: LB', 'Season 2: LB+AMP','Season 3: LB'},'FontSize',18);
legend boxoff

export_fig '../../figures/Fig-6B.pdf';


%% 

%stdGFPs=zeros(numDays, length(As), length(reps));
meanGFPsP=zeros(numDays, length(AsP), nreps);
lim_area=[25, 35];
%cvGFPs=zeros(numDays, length(As), length(reps));
for day=1:numDays
    if(day==1)
        fileName=pfileName1;
    else
        fileName=pfileName;
    end
    
    filePath=[dataPath,fileName];
    data = readtable(filePath,'Delimiter', '\t');

    dayData=data(data.season==day ,:);
    for iA=1:length(AsP)
        dayAsData=dayData(dayData.doseNumber==colsP(iA),:);
         
        for irep=1:nreps
            repname=string(['rep-',sprintf('%d',irep)]);
            asrepData=dayAsData(dayAsData.replicate==repname,:);
            thisGFP=asrepData.Intensity_MC_Ch02(asrepData.Area_M02 > lim_area(1) & asrepData.Area_M02 < lim_area(2) & asrepData.AspectRatio_M02>lim_aspect(1) & asrepData.AspectRatio_M02<lim_aspect(2),:);           
            meanGFPsP(day, iA, irep)=mean(thisGFP);
            %cvGFPs(day, iA, irep)=getCV(thisGFP);
            %stdGFPs(day, iA, irep)=std(thisGFP);
        
        end
        
    end
end
meanGFPsC=zeros(numDays, length(AsC), nreps);
lim_area=[1, 35];
%cvGFPs=zeros(numDays, length(As), length(reps));
for day=1:numDays
     
    fileName=cfileName;
    
    filePath=[dataPath,fileName];
    data = readtable(filePath,'Delimiter', '\t');

    dayData=data(data.season==day ,:);
    for iA=1:length(AsC)
        dayAsData=dayData(dayData.doseNumber==colsC(iA),:);
         
        for irep=1:nreps
            repname=string(['rep-',sprintf('%d',irep)]);
            asrepData=dayAsData(dayAsData.replicate==repname,:);
            thisGFP=asrepData.Intensity_MC_Ch02(asrepData.Area_M02 > lim_area(1) & asrepData.Area_M02 < lim_area(2) & asrepData.AspectRatio_M02>lim_aspect(1) & asrepData.AspectRatio_M02<lim_aspect(2),:);           
            
            if(length(thisGFP)>1000)
                meanGFPsC(day, iA, irep)=mean(thisGFP);
            %cvGFPs(day, iA, irep)=getCV(thisGFP);
            %stdGFPs(day, iA, irep)=std(thisGFP);
            end
        end
        
    end
end




%% plot mean
figure();% clf('reset');
set(gcf,'DefaultLineLineWidth',1); set(gcf, 'color', 'white'); hold all

CTP=cbrewer('seq', 'YlOrRd', length(AsP));
for iA=1:length(AsP)
    time_meanGFP=zeros(1,numDays);
    for day=1:numDays
        meanRep=zeros(1,nreps);
        for irep=1:nreps
            meanRep(irep)=meanGFPsP(day, iA, irep);
        end
        time_meanGFP(day)=mean(meanRep);
    end
    plot(1:numDays, time_meanGFP,  'ko-','MarkerFaceColor', CTP(iA,:),'Color',CTP(iA,:),'MarkerSize',10);
    
end
CTC=cbrewer('seq', 'YlOrRd', length(AsC));
for iA=1:length(AsC)
    time_meanGFP=zeros(1,numDays);
    for day=1:numDays
        meanRep=zeros(1,nreps);
        for irep=1:nreps
            meanRep(irep)=meanGFPsC(day, iA, irep);
        end
        time_meanGFP(day)=mean(meanRep);
    end
    plot(1:numDays, time_meanGFP,  'ko-','MarkerFaceColor', CTC(iA,:),'Color',CTC(iA,:),'MarkerSize',10);
end

set(gca,'fontsize',18);
xlabel('Time (seasons)');
ylabel('Mean GFP intensity (per cell)');
colormap(CTP);

strAsP={'0','1/128','1/64','1/32', '1/16','1/8','1/4','1/2','1','2'};

cbh=colorbar();
cbtick=linspace(0,1,length(AsP)+1);
cbtick=cbtick(1:end-1)+.05;
set(cbh,'YTick',cbtick)
set(cbh,'YTickLabel',strAsP)
xticks([1,2,3]);

axis([.9 3.1 -1e4 1.5e5]);
%axis([.9 numDays+.1 0 1e4]);


ylabel(cbh,'Drug concentration (units of MIC)','FontSize',20);

export_fig '../../figures/Fig-6C.pdf';





