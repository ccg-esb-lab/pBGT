clc
clear all
close all

run('../matlab_code/lib/addpath_recurse');
addpath_recurse('../matlab_code/lib/');


%%

dataDir='../../data/uJ_data/pBGT-Ramp/data/';  %tmp
frame2min=2;
death=[152, 84, 136, 128];

%%
figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'fontsize',20);

set(gcf, 'Units','normalized','Position',[0. 0. 0.8 0.18]);

str={};
p=[];

for i=1:4

    fileName=[dataDir,'xy3_segmented_cell',num2str(i),'.csv'];
    M = readtable(fileName);
    idead=find(M.Slice==death(i));
    lineH=plot((frame2min.*M.Slice(1:idead)),(M.Mean(1:idead)),'-','LineWidth',3); hold all;
    color = get(lineH, 'Color');
    p(i)=plot(frame2min.*M.Slice(idead), M.Mean(idead), 'o', 'MarkerFaceColor',color, 'MarkerEdgeColor',color);
    str{i}=['',num2str(i)];
end
ylim([0 300]);
set(gca,'FontSize',16);
xlabel('Time (minutes)');
ylabel('GFP intensity (a.u.)');
xticks(0:60:360)
xlim([0, 320]);
legend(p, str,'Location','EastOutside')

export_fig('../../figures/Fig-2B.pdf');