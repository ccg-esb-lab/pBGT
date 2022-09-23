clc
clear all
close all

run('../matlab_code/lib/addpath_recurse');
addpath_recurse(['../matla' ...
    'b_code/lib/']);


%% IMPORT DATA

dataDir='../../data/Figure_3/';  %tmp
fileName='strains_survival_prob.csv';

M = csvread([dataDir,fileName],1);

deltaT=60.*M(:,1);
MG_surv=100.*M(:,2);
pBGT_surv=100.*M(:,3);
G54U_surv=100.*M(:,4);
G55U_surv=100.*M(:,5);


CT=cbrewer('qual', 'Pastel2', 4);
color_MG=[0 0 0];
color_pBGT=CT(1,:);
color_G54U=CT(2,:);
color_G55U=CT(3,:);
    
    
%%

figure(); clf('reset');set(gcf,'DefaultLineLineWidth',3); set(gcf, 'color', 'white'); hold all
%set(gcf,'Position',[1000         918         560         240])

p_pBGT=plot(deltaT, pBGT_surv, '-', 'Color',color_pBGT); hold on;
p_G54U=plot(deltaT, G54U_surv, '-', 'Color',color_G54U); hold on;
p_G55U=plot(deltaT, G55U_surv, '-', 'Color',color_G55U); hold on;
p_MG=plot(deltaT, MG_surv, '-', 'Color',color_MG); hold on;

plot([90,90],[0,100],'r:')

set(gca,'fontsize',18);

%xticks(deltaT);
xlim([0 deltaT(end)]);
xticks(0:60:240);
xlabel('Duration of antibiotic exposure (minutes)','FontSize',20);
ylabel('Survival probability (%)','FontSize',20);

legend([p_pBGT,p_G54U, p_G55U, p_MG],{'MG:pBGT','MG:G54U','MG:G55U','MG:GT'},'FontSize',18,'Location','NorthEast');
legend boxoff

export_fig '../../figures/Fig-3B.pdf'


