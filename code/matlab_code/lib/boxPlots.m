function boxPlots(f,datasets,featureLabels,labelAngle,colors,axisType,showGrid,fontSize,outliers,Legend,curvature,xLb,yLb)

% boxPlots: plots boxplots.
% 
% Arguments
%   datasets:         cell array containing datasets.
%   FeatureLabels:    feature labels.
%   labelAngle:       angle of xlabel.
%   Colors:           array containing color triplets for each dataset.
%   axisType:         'linear' to set y axis linear.
%                     'log' to set y axis logarithmic.
%   showGrid:         'on' to show grid.
%                     'off' to do not show grid.
%   fontSize:         scalar to set font size.
%   outliers:         up_20_d to show only outliers lying above the third
%                     quartile. 20 is the size of outlier dots. d is the
%                     type of outlier dots (eg. d for diamond).
%                     down_20_d to show only outliers lying below the first
%                     quartile.
%                     both_20_d to show all outliers.
%   Legend:           'Setname1_Setname2_...SetnameN' is the name of each
%                     Dataset separated by underline symbol.
%   Curvature:        a scalar value range from 0 to 1 that adds curved
%                     edges to rectangle of boxplots.
%   xLb:              string of xlabel axis property.
%   yLb:              string of ylabel axis property.
% 
% Example:
% 
% Copy the code below without '%' and paste it into command window and press Enter
% -------------------------------------------------------------------------
% OBS. If the code below is uncomented the function will fall into an
% infinte recursion bug
% -------------------------------------------------------------------------
% a=40*randn([1000,3]);
% b=50*randn([1000,3]);
% c=60*randn([1000,3]);
% features={a,b,c};
% 
% boxPlots(features,{'Feature 01','Feature 02','Feature 03'},45,[89 117 164; 95 158 110; 45 87 96],'linear','on',14,'both_50_d','DATASET 01_DATASET 02_DATASET 03',0.1,'Features','Values')
% 
% To export a high resolution image of the boxplot type the code below
% print('nameOfTheImage','-dpng','-r300')
% '-dpng' for png format
% '-r300' for 300 dpi

[nColors,~]=size(colors);
Sets=length(datasets);

if(nColors~=Sets)
    disp('The number of colors must be the same as the number of sets.')
    return
end

[~,nFeaturesL]=size(featureLabels);
[~,nFeatures]=size(cell2mat(datasets(1)));

if(nFeaturesL~=nFeatures)
    disp('The number of feature labels must be the same as the number of features.')
    return
end

f=figure;
ax=axes;

set(ax,'XGrid','off')
set(ax,'YGrid','off')

if(strcmp(axisType,'linear'))
%    set('YScale','linear')
elseif(strcmp(axisType,'log'))
 %   set('YScale','log')
else
    disp('The number os colors must be the same as the number of sets.')
    delete(f)
    return
end

XAxisSub = Sets*nFeatures+1;

xl = xlim;
xLimMin=xl(1);
xLimMax=xl(2);
xWidth = xLimMax-xLimMin;
xWidthSub = xWidth/XAxisSub;
pos=[];
fLabel=[];

for j=1:nFeatures
    for i=1:Sets
        %sort the feature matrix in ascending order
        FeaturesSorted = sort(cell2mat(datasets(i)),1);
        %median of the entire set
        secondQuartile = median(FeaturesSorted(:,j));
        subsetFeaturesDOWN=FeaturesSorted(:,j);
        subsetFeaturesUP=FeaturesSorted(:,j);
        subsetFeaturesDOWN(FeaturesSorted(:,j)>secondQuartile)=[];
        subsetFeaturesUP(FeaturesSorted(:,j)<secondQuartile)=[];
        firstQuartile = median(subsetFeaturesDOWN);
        thirdQuartile = median(subsetFeaturesUP);
        max_ = thirdQuartile+(1.5*(thirdQuartile-firstQuartile));
        min_ = firstQuartile-(1.5*(thirdQuartile-firstQuartile));
        
        SET=cell2mat(datasets(i));
        Max_(j,i) = max(SET(:,j));
        Min_(j,i) = min(SET(:,j));
        
        hold on
        %plot the box
        rectangle('Position',[((xWidthSub*(i+Sets*(j-1)))-((xWidthSub*0.8)/2)),firstQuartile,xWidthSub*0.8,(thirdQuartile-firstQuartile)],...
            'FaceColor',convertColorRange(colors(i,:)),...
            'EdgeColor','k',...
            'LineWidth',1.5,...
            'curvature',[curvature curvature])
        
        %plot the vertical lines
        plot([xWidthSub*(i+Sets*(j-1)),xWidthSub*(i+Sets*(j-1))],[min_,firstQuartile],'k-','LineWidth',1.5)
        plot([xWidthSub*(i+Sets*(j-1)),xWidthSub*(i+Sets*(j-1))],[thirdQuartile,max_],'k-','LineWidth',1.5)
        
        %plot the horizontal max and min lines
        plot([((xWidthSub*(i+Sets*(j-1)))-((xWidthSub*0.5)/2)),((xWidthSub*(i+Sets*(j-1)))-((xWidthSub*0.5)/2))+(xWidthSub*0.5)],[min_,min_],'k-','LineWidth',1.5)
        plot([((xWidthSub*(i+Sets*(j-1)))-((xWidthSub*0.5)/2)),((xWidthSub*(i+Sets*(j-1)))-((xWidthSub*0.5)/2))+(xWidthSub*0.5)],[max_,max_],'k-','LineWidth',1.5)
        
        %plot the second quartile horizontal line
        plot([((xWidthSub*(i+Sets*(j-1)))-((xWidthSub*0.8)/2)),((xWidthSub*(i+Sets*(j-1)))-((xWidthSub*0.8)/2))+(xWidthSub*0.8)],[secondQuartile,secondQuartile],'k-','LineWidth',1.5)
        
        %define the outlier limits
        
        upperOutlierLim=thirdQuartile+1.5*(thirdQuartile-firstQuartile);
        lowerOutlierLim=firstQuartile-1.5*(thirdQuartile-firstQuartile);
        
        Out_=strsplit(outliers,'_');
        
        if(strcmp(Out_(1),'up'))
            outLUP=SET(:,j); outLUP(outLUP<upperOutlierLim)=[];
            x=((xWidthSub*(i+Sets*(j-1)))*ones([length(outLUP),1]));
            y=outLUP;
            st_up=scatter(x,y,str2double(Out_(2)),char(Out_(3)));
            set(st_up,'MarkerEdgeColor',[0 0 0])
            set(st_up,'MarkerFaceColor',convertColorRange(colors(i,:)))
        elseif(strcmp(Out_(1),'down'))
            outLDOWN=SET(:,j); outLDOWN(outLDOWN>lowerOutlierLim)=[];
            x=((xWidthSub*(i+Sets*(j-1)))*ones([length(outLDOWN),1]));
            y=outLDOWN;
            st_down=scatter(x',y',str2double(Out_(2)),char(Out_(3)));
            set(st_down,'MarkerEdgeColor',[0 0 0])
            set(st_down,'MarkerFaceColor',convertColorRange(colors(i,:)))
        elseif(strcmp(Out_(1),'both'))
            outLUP=SET(:,j); outLUP(outLUP<upperOutlierLim)=[];
            x=((xWidthSub*(i+Sets*(j-1)))*ones([length(outLUP),1]));
            y=outLUP;
            st_up=scatter(x,y,str2double(Out_(2)),char(Out_(3)));
            set(st_up,'MarkerEdgeColor',[0 0 0])
            set(st_up,'MarkerFaceColor',convertColorRange(colors(i,:)))
            outLDOWN=SET(:,j); outLDOWN(outLDOWN>lowerOutlierLim)=[];
            x=((xWidthSub*(i+Sets*(j-1)))*ones([length(outLDOWN),1]));
            y=outLDOWN;
            st_down=scatter(x',y',str2double(Out_(2)),char(Out_(3)));
            set(st_down,'MarkerEdgeColor',[0 0 0])
            set(st_down,'MarkerFaceColor',convertColorRange(colors(i,:)))
        end
        pos=[pos,(xWidthSub*(i+Sets*(j-1)))];
    end
end

if(min(Min_(:))<0)
    min_global=min(Min_(:));
    max_global=max(Max_(:));
    ylim([min_global-abs(min_global*0.2),max_global+max_global*0.6])
else
    min_global=min(Min_(:));
    max_global=max(Max_(:));
    ylim([min_global-abs(min_global*0.2),max_global+max_global*0.6])
end

xlim([0,1])
%set(ax,'Box','on')

if(strcmp(showGrid,'on'))
    grid on
else
    grid off
end

xticks(pos)
for i=1:Sets
    fLabel=[fLabel,[featureLabels(i),featureLabels(i),featureLabels(i)]];
end
xticklabels(fLabel)
xtickangle(labelAngle)
set(gca,'FontSize',fontSize)
xlabel(xLb)
ylabel(yLb)

%Plot the legend
xl = xlim; width_axis = xl(2)-xl(1); width_legend_box = width_axis/Sets;
yl = ylim; height_axis = yl(2)-yl(1); height_legend_box = (height_axis*0.05);
%L_=strsplit(Legend,'_');

%rectangle('Position',[xl(1)+(width_axis*0.05),yl(2)-(height_axis*0.15),(width_axis*0.9),(height_axis*0.1)],'FaceColor','w','EdgeColor','k','LineWidth',1.5)
%for i=1:Sets
%    rectangle('Position',[xl(1)+(width_axis*0.1)+(i-1)*(width_legend_box*0.9),yl(2)-(height_axis*0.1)-(height_legend_box/2),width_legend_box*0.1,height_legend_box],'FaceColor',convertColorRange(colors(i,:)),'EdgeColor','k','LineWidth',0.5)
%    text(xl(1)+(width_axis*0.12)+(width_legend_box*0.1)+(i-1)*(width_legend_box*0.9),yl(2)-(height_axis*0.075)-(height_legend_box/2),L_(i),'FontUnits','normalized')
%end

end

function color_ = convertColorRange(colors)

R = colors(1); R_ = R/255;
G = colors(2); G_ = G/255;
B = colors(3); B_ = B/255;

color_ = cat(2,R_,G_,B_);

end