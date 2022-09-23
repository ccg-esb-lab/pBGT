clc
clear all
close all

run('../matlab_code/lib/addpath_recurse');
addpath_recurse('../matlab_code/lib/');


%% LOAD DATA

dataPath='../../data/Figure_1/';
fileName='all_strains';

filePath=[dataPath,fileName];
data = readtable(filePath,'Delimiter', '\t');

xmin=1e1;
xmax=10e5;
xdelta=5000;
nreps=4;

lbl_strains={'MG','MG:GT','MG/pBGT','MG/G55U','MG/G54U'};
mic_strains=[4,512,8192,32768,32768];
rel_mic_strains=mic_strains./512;

fitness_strains=[1.01,1,.943,.793,.557];

cost_strains=100.*(1-fitness_strains);
pcn_strains=[0,1,19.12,44.50,88.93];
std_pcn_strains=[0, 0, 1.56, 3.81, 15.65];

vmax_strains=[0.48, 0.44, 0.42, .36, .53];  %FIND DATA - CHARLY



%% GFP vs PCN
labels={'MG:GT','MG/pBGT','MG/pBGT G54U','MG/pBGT G55U','MG1655'};
istrains=[1,2,3,4,5]; %plasmid-bearing strains
%dfstrains={'G54U','G55U-','G55U','MGGT','MG','pBGT'};
dfstrains={'MG','MGGT','pBGT','G54U','G55U'};

figure(); clf('reset');set(gcf,'DefaultLineLineWidth',3); set(gcf, 'color', 'white'); hold all


CT=[.75 .75 .75; cbrewer('qual', 'Pastel2', length(labels)-1) ];


cvGFP=zeros(1,length(istrains));
meanGFP=zeros(1,length(istrains));
stdGFP=zeros(1,length(istrains));
steGFP=zeros(1,length(istrains));

for i=(1:length(istrains))
    istrain=istrains(i);
    this_strain=string(dfstrains(i));
    disp(this_strain)
    datastrain =data(data.strain==this_strain,:);
    disp(head(datastrain))
    thisGFP=datastrain.GFP;
    
    meanGFP(i)=mean(thisGFP);
    stdGFP(i)=std(thisGFP);
    cvGFP(i)=stdGFP(i)/meanGFP(i);
    steGFP(i)=std(thisGFP)/sqrt(nreps);
  
    
end
 
A=pcn_strains(istrains);
B=meanGFP;

b = polyfit(A, B, 1);
f = polyval(b, A);
Bbar = mean(B);
SStot = sum((B - Bbar).^2);
SSreg = sum((f - Bbar).^2);
SSres = sum((B - f).^2);
R2 = 1 - SSres/SStot;
R = corrcoef(A,B);
Rsq = R(1,2).^2;

plot(f, A,'k:','LineWidth',1);


rmax=50;
rmin=6;
for i=1:length(istrains)
    
    plot([meanGFP(istrains(i)), meanGFP(istrains(i))], [pcn_strains(istrains(istrains(i)))-std_pcn_strains(istrains(i)), pcn_strains(istrains(i))+std_pcn_strains(istrains(i))],  '-k', 'MarkerFaceColor', 'k', 'MarkerSize',10,'LineWidth',1);
    plot([meanGFP(istrains(i))-steGFP(istrains(i)), meanGFP(istrains(i))+steGFP(istrains(i))], [pcn_strains(istrains(i)), pcn_strains(istrains(i))],  '-k', 'MarkerFaceColor', 'k', 'MarkerSize',10,'LineWidth',1);
    plot(meanGFP(istrains(i)), pcn_strains(istrains(i)),  'o','MarkerFaceColor',CT(i,:),'MarkerEdgeColor',[.5 .5 .5],'MarkerSize', rmin+rmax*mic_strains(i)/mic_strains(4),'LineWidth',1);

end

set(gca,'fontsize',20);
ylabel('PCN','FontSize',24);
xlabel('Mean GFP intensity','FontSize',24);

text(2e3,100,['R^2=',num2str(Rsq)],'FontSize',20,'VerticalAlignment','top');
ylim([-5 100]);
xlim([-.5e4 1.5e5])

export_fig('../../figures/Fig-1D.pdf');


