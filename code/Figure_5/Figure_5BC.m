clc
clear all
close all

run('../matlab_code/lib/addpath_recurse');
addpath_recurse('../matlab_code/lib/');

%%

dataPath_MGGT='../../data/uJ_data/MGGT-Pulse/data/lineages_status/';
frames_MGGT=10:34;

dataPath_pBGT='../../data/uJ_data/pBGT-Pulse/data/lineages_status/';
frames_pBGT=20:44;

frames=1:frames_pBGT(end)-frames_pBGT(1)+1;
num_frames=length(frames_pBGT);
frame2min=10;

color_light_normal=[0.38431     0.61176           1]; %'#629CFF'
color_light_elongated=[0.98824     0.97647     0.58431]; %'#FCF995'
color_light_dead=[0.97255     0.46275     0.42745]; %'#F8766D'

%% MGGT

xy_files_MGGT=dir([dataPath_MGGT,'*.csv']);

num_normal_MGGT_reps=zeros(length(xy_files_MGGT),num_frames);
num_elongated_MGGT_reps=zeros(length(xy_files_MGGT),num_frames);
num_dead_MGGT_reps=zeros(length(xy_files_MGGT),num_frames);

for ifile=1:length(xy_files_MGGT)
    fileName_MGGT=xy_files_MGGT(ifile).name;
    disp(['Loading ',fileName_MGGT]);
    filePath_MGGT=strcat(dataPath_MGGT,'',fileName_MGGT);
    
    DATA_MGGT = readtable(filePath_MGGT,'ReadVariableNames', true);

    for iframe=1:num_frames

        %CHROMOSOMAL
        this_frame_MGGT=frames_MGGT(iframe);
        DATA_frame_MGGT=DATA_MGGT(ismember(DATA_MGGT.frame,this_frame_MGGT),:) ;

        num_elongated_MGGT_reps(ifile,iframe)=height(DATA_frame_MGGT(ismember(DATA_frame_MGGT.state, 2),:));
        num_normal_MGGT_reps(ifile,iframe)=height(DATA_frame_MGGT(ismember(DATA_frame_MGGT.state, 1),:));
        num_dead_MGGT_reps(ifile,iframe)=height(DATA_frame_MGGT(ismember(DATA_frame_MGGT.state, 3),:));

    end
end


num_normal_MGGT=sum(num_normal_MGGT_reps);
num_dead_MGGT=sum(num_dead_MGGT_reps);
num_elongated_MGGT=sum(num_elongated_MGGT_reps);

%% pBGT

xy_files=dir([dataPath_pBGT,'*.csv']);

num_normal_pBGT_reps=zeros(length(xy_files),num_frames);
num_elongated_pBGT_reps=zeros(length(xy_files),num_frames);
num_dead_pBGT_reps=zeros(length(xy_files),num_frames);

for ifile=1:length(xy_files)
    fileName_pBGT=xy_files(ifile).name;
    disp(['Loading ',fileName_pBGT]);
    filePath_pBGT=strcat(dataPath_pBGT,'',fileName_pBGT);
    DATA_pBGT = readtable(filePath_pBGT,'ReadVariableNames', true);

    for iframe=1:num_frames

        %PLASMID
        this_frame_pBGT=frames_pBGT(iframe);
        DATA_frame_pBGT=DATA_pBGT(ismember(DATA_pBGT.frame,this_frame_pBGT),:) ;

        num_elongated_pBGT_reps(ifile,iframe)=height(DATA_frame_pBGT(ismember(DATA_frame_pBGT.state, 2),:));
        num_normal_pBGT_reps(ifile,iframe)=height(DATA_frame_pBGT(ismember(DATA_frame_pBGT.state, 1),:));
        num_dead_pBGT_reps(ifile,iframe)=height(DATA_frame_pBGT(ismember(DATA_frame_pBGT.state, 3),:));

    end
end

num_normal_pBGT=sum(num_normal_pBGT_reps);
num_dead_pBGT=sum(num_dead_pBGT_reps);
num_elongated_pBGT=sum(num_elongated_pBGT_reps);


%% PLOT POPULATION FRACTIONS

figure(); clf('reset'); set(gcf, 'color', 'white'); hold all
set(gcf,'Units','normalized','Position',[0.1 .1 .5 .4])

num_total_pBGT=num_normal_pBGT+num_elongated_pBGT+num_dead_pBGT;

times=(frames-1)*frame2min;


subaxis(2,1,1,'SpacingVert',0.005,'PaddingTop',0.1); 
num_total_MGGT=num_normal_MGGT+num_elongated_MGGT+num_dead_MGGT;
normal_upper=100*num_normal_MGGT./num_total_MGGT;
normal_lower=zeros(1,num_frames);
elongated_upper=100*(num_normal_MGGT+num_elongated_MGGT)./num_total_MGGT;
elongated_lower=100*num_normal_MGGT./num_total_MGGT;
dead_upper=100*(num_normal_MGGT+num_elongated_MGGT+num_dead_MGGT)./num_total_MGGT;
dead_lower=100*(num_normal_MGGT+num_elongated_MGGT)./num_total_MGGT;

plot(times, normal_upper, '-k', 'LineWidth',1); hold on;
plot(times, elongated_upper, '-k', 'LineWidth',1); hold on;

jbfill(times,normal_lower,normal_upper,color_light_normal,color_light_normal,1,0.6);
jbfill(times,elongated_lower,elongated_upper,color_light_elongated,color_light_elongated,1,0.6);
jbfill(times,dead_lower,dead_upper,color_light_dead,color_light_dead,1,0.6);

xticks([]);
ylim([0, 100]);
xlim([0,240]);
set(gca,'fontsize',20);


subaxis(2,1,2,'SpacingVert',0.005,'PaddingTop',0.1);  

normal_upper=100*num_normal_pBGT./num_total_pBGT;
normal_lower=zeros(1,num_frames);
elongated_upper=100*(num_normal_pBGT+num_elongated_pBGT)./num_total_pBGT;
elongated_lower=100*num_normal_pBGT./num_total_pBGT;
dead_upper=100*(num_normal_pBGT+num_elongated_pBGT+num_dead_pBGT)./num_total_pBGT;
dead_lower=100*(num_normal_pBGT+num_elongated_pBGT)./num_total_pBGT;

plot(times, normal_upper, '-k', 'LineWidth',1); hold on;
plot(times, elongated_upper, '-k', 'LineWidth',1); hold on;

jbfill(times,normal_lower,normal_upper,color_light_normal,color_light_normal,1,0.6);
jbfill(times,elongated_lower,elongated_upper,color_light_elongated,color_light_elongated,1,0.6);
jbfill(times,dead_lower,dead_upper,color_light_dead,color_light_dead,1,0.6);

xticks(0:60:240);
xlim([0,240]);
ylim([0, 100]);
xlabel('Time (minutes)','FontSize',20);
set(gca,'fontsize',18);
%ylabel('Relative abundance (%)','FontSize',20);

%print(gcf, '-dpdf', 'figures/Fig-5BD.pdf');
export_fig('../../figures/Fig-5BD.png');
export_fig('../../figures/Fig-5BD.pdf');


