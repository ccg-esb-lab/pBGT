clc
clear all
close all

run('../matlab_code/lib/addpath_recurse');
addpath_recurse('../matlab_code/lib/');


%% LOAD DATA
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
dataPath0='../../data/LB-AMP-LB/';

%***** Plasmid


fileNameR='pBGT_LB-AMP-LB_AMNIS_reread.csv';
fileName='pBGT_LB-AMP-LB_AMNIS.csv';
strain='plasmid';
reps={'A','B','F'};
cols=[2,3,4,5,6,7,8,9,10,11,12];
As=[0,1/256, 1/128, 1/64, 1/32, 1/16,1/8,1/4,1/2,1];
strAs={'AMP=0','1/256 MIC', '1/128 MIC','1/64 MIC', '1/32 MIC','1/16 MIC','1/8 MIC','1/4 MIC','1/2 MIC','MIC'};
lim_aspect=[.2, 1];
lim_area=[25, 35];
xmax=.75e6;
xmin=.9e4;
xdelta=.5e4;
col=[0.3020    0.6863    0.2902];

numDays=3;

day=2;
filePath=[dataPath0,fileName];
data = readtable(filePath,'Delimiter', '\t');
disp(head(data))
seasondata=data(data.season==day ,:);
disp(head(seasondata));


%% dfdf
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
clc

figure(1); clf('reset');set(gcf,'DefaultLineLineWidth',1); set(gcf, 'color', 'white'); hold all
set(gcf,'Units','normalized','Position',[0.1 .1 .2 .8])


CT2=cbrewer('seq', 'Greens', 100);

xmax2=2.5e5;
ymax=1;


for iA=1:length(As)
    iAx=double(iA)+1;
       
    asData=seasondata(seasondata.doseNumber==iAx,:);
    
    repH=[];
    
    for irep=1:length(reps)
        repname=string(['rep-',sprintf('%d',irep)]);
        asrepData=asData(asData.replicate==repname,:);
        
        %subaxis(length(As),numDays, day, iA,  'SpacingVert',0.01, 'SpacingHoriz',0.02,'PaddingBottom',.0);
        thisGFP=asrepData.Intensity_MC_Ch02(asrepData.Area_M02 > lim_area(1) & asrepData.Area_M02 < lim_area(2) & asrepData.AspectRatio_M02>lim_aspect(1) & asrepData.AspectRatio_M02<lim_aspect(2),:);
        
        disp([num2str(As(iA)),'MIC   Day',num2str(day),'   ', fileName,': ',num2str(length(thisGFP)),' cells']);
        
        e=xmin:xdelta:xmax2;
        h=hist(thisGFP(thisGFP>xmin & thisGFP<xmax2), e);
        
        %hnorm=h;
        hnorm=h./max(h);
        %semilogx(e,hnorm,'LineWidth',1,'Color',[.8 .8 .8]); hold on;
        
        repH=[repH; hnorm];
    
    end
    %semilogx(e,mean(repH),'k','LineWidth',2); hold on;
    ypos=(length(As)-iA)*.66;
    %b=bar(e,ypos+mean(repH), .9, 'FaceColor',CT(iA,:)); hold on;
    %b(1).BaseValue = ypos;
    %axis([xmin .9*xmax 0 1]);
    
    axis tight
    thisColor=CT2(floor(mean(1.5*mean(repH))*100),:);
    
    jbfill(e,ones(1,length(repH)).*ypos+mean(repH),ones(1,length(repH)).*ypos,thisColor,'k',1,1);
    
    %area(e,ypos+mean(repH),ypos); hold on;
    %plot([0 xmax],[ypos ypos],'k-','LineWidth',1); hold on;
    set(gca,'fontsize',18);
    %xticks(10.^(1:5));
    
    %if iA==ceil(length(As)/2)+1
    %    if  day==2
    %        ylabel('          Frequency','FontSize',24);
    %    end
    %end
    %if iA<length(As)
    yticks([]);
    %end
    
    if iA==1
        text(xmax2*.99, ypos*1+.15, ['No drug'],'FontSize',16,'HorizontalAlignment','right','VerticalAlignment','bottom');
    else
        text(xmax2*.99, ypos*1+.15, [strAs{iA} ''],'FontSize',16,'HorizontalAlignment','right','VerticalAlignment','bottom');
    end
    
    if day==1
        title([strAs{iA}]);
    end
    
    if iA==length(As)
        xlabel('GFP intensity (a.u.)','FontSize',20);
        
    else
        %xticks([]);
    end
    
end

%print('-dpng','-r300','figures/AMNIS_ridgelineP.png');
export_fig '../../figures/Fig-1E.pdf';



