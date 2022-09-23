%% Demonstration of boxplotx
% Copyright (c) 2016, W J Whiten  BSD License
%
% The program boxplotx is for those who do not have the statistics
% toolbox 
%
% For details see doc boxplotx
%%
%

% % Fixes for latex
% % \usepackage[a4paper, total={6in, 9in}]{geometry} 
% % \pagebreak
% % \LARGE 
% % \Huge 
% % \huge 
% % \textit{\textbf{  }}
% % 
% % \setlength{\itemsep}{2ex}
% % \item[]
% % 
% % \includegraphics [width=6in]{example1_02.eps}

%% Simple usage 
% Boxplots are generated from the columns of a matrix
%
%
rng(0)  % reset random generator
x1=randn(200,10);
boxplotx(x1)
%%
%
%% Additional box plots can be added
% Set hold on for the current graph and call boxplotx again
%
% Plotting options for lines and points can be changed
%
%
hold on
boxplotx(randn(100,2),'lines','g-','points','r*')

%%
%
%% Box plots can have different numbers of points
% Also position and width can be set using second and third default
% arguments
%
% Each boxplot can be labelled
%
%
x2={randn(50,1),randn(1,300),randn(100,1)};
figure
boxplotx(x2,[1,4,6],[0.5,2,1])

%%
%
%% Bar percentages and boxplot labels
% Position and width can can be given as name & value options
%
% Labels can be put under the boxplots
%
figure
boxplotx(x2,'qpct',30,'bpct',5,'lines',{'r-','b-','g-'},  ...
    'xpos',[1,4,6],'width',[0.5,2,1],  ...
    'labels',{'50 Points','200 Points','100 Points'})

%%
%
%
% For more detail on options see doc boxplotx

