close all; clear all; clear classes; clc

cd('/home/kpalmer/AnacondaProjects/Localisation/Scripts')
load('PMRF_localizations_parm_1250_1600_18Jan17_010650_1159_2317.mat')
load('PMRF_array_struct_20180108.mat')
master = array_struct(1).master;
% Pre-allocate the space


% Pre-allocate the cross correlation structure
Call_space =[];

%% Cluster PMRF using ladders
close all
fs =6000;
array_idx=1
LSQ_plot_thres =.0008; % Plotting only
cor_thresh = .2;
time_cut = 20*60; %30 minutes
NChidlHyd = 2; % Number of pairs of hydrophones


% Set up data parameters
PRMRFCluster  = clusterDataGPLclass();
PRMRFCluster.array_struct = array_struct;
PRMRFCluster.localize_struct =localize_struct;
PRMRFCluster.fs =fs;
PRMRFCluster.time_cut =time_cut
PRMRFCluster.CrossScore = 100000
PRMRFCluster.cutoff =.006

% Run the clustering algorithim up to the the similarity matrix
simMatLowMemory(PRMRFCluster)
% Create the Cluster ID's
updateClusterID(PRMRFCluster)

%% Plot clusters and the un-clustered data

[k1 k2]= find(PRMRFCluster.localize_struct.hyd(19).score <LSQ_plot_thres);


title1 =[ 'Localized Data Clusterd Using Ladder Structure and Mean'...
    newline 'Correlation Between ', num2str(NChidlHyd),' Hydrophone Pairs'];
title2 = 'Localized Data Plotted Against Time of Detection';


x = PRMRFCluster.arrivalTable.Loc(k2,1);
y = PRMRFCluster.arrivalTable.Loc(k2,2);

chain = PRMRFCluster.chains;
Cluster_id = PRMRFCluster.Cluster_id;

color_var_drift= zeros(size(x,1));
ColorVals = jet(length(k2));
colorscale2 = (PRMRFCluster.localize_struct.hyd(19).rtimes(k2)- ...
    min(PRMRFCluster.localize_struct.hyd(19).rtimes(k2)))/fs;



% Create color variable
cluster_n = 1; % initial cluster color (other than 0)

% Color all localizations in the cluster the same
for ii=1:length(chain)
    
    if chain(ii).n >1
        color_var_drift([ii, chain(ii).index]) = (log10(chain(ii).n))+.001;
        cluster_n =cluster_n+1;
    end
end


figure
fig = get(groot,'CurrentFigure');


x_trim_drift=x;
y_trim_drift=y;
color_var_drift = color_var_drift(k2);


cbar1label = 'log10(Number of Detections in Cluster)';
cbar2label = 'Elapsed Time Since First Detection (sec)';

% text labels, none for unknown data
txt_str2 =[''];
x_txt = max(x)-.1
y_txt = min(y) +.01


txt_str1 =[num2str(length(unique(Cluster_id))), ' Groups Predicted'];
y_txt = min(y);



%Plot data
fig
subplot(2,1, 1)
scatter(x_trim_drift, y_trim_drift, 40, color_var_drift, 'filled');
title(title1)
ylabel('Lat'); xlabel('Lon')
if ~isempty(LSQ_plot_thres) 
    colormap((jet))
else
    colormap(lines(ncolors))
end
text(x_txt, y_txt, txt_str1)
h = colorbar;
ylabel(h, cbar1label)
hold off

fig
subplot(2,1, 2)
scatter(x_trim_drift, y_trim_drift, 40, colorscale2, 'filled');
title(title2)
text(x_txt, y_txt, txt_str2)
ylabel('Lat'); xlabel('Lon')
h = colorbar;
ylabel(h, cbar2label)



