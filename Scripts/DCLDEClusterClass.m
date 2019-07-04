1+1
%%
clear all; close all; clc

close all; clear all
clear classes; clc

whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2})
load(whereAmI{3})
clear whereAmI


hydrophone_struct= struct();
for ii=1:size(dclde_2013_meta,1)
    hydrophone_struct(ii).name = num2str(dclde_2013_meta(ii,1));
    hydrophone_struct(ii).location = dclde_2013_meta(ii,[11:12]);
    hydrophone_struct(ii).depth= abs(dclde_2013_meta(ii, 13));
    hydrophone_struct(ii).channel=ii;
end
hyd_arr = struct2cell(hydrophone_struct);
hyd_arr =vertcat(hyd_arr{2,:,:});

parent =5;
fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;

% Set up data parameters
DCLDE2013Cluster  = clusterDataGPLclass();

DCLDE2013Cluster.array_struct = array_struct;

DCLDE2013Cluster.localize_struct =localize_struct;

DCLDE2013Cluster.fs =fs;
DCLDE2013Cluster.Parent_hyd = 5;
DCLDE2013Cluster.time_cut =time_cut
DCLDE2013Cluster.CrossScore = 100000
DCLDE2013Cluster.cutoff =.006
% Run the clustering algorithim up to the the similarity matrix
simMatLowMemory(DCLDE2013Cluster)
% Create the Cluster ID's
updateClusterID(DCLDE2013Cluster)

%% Hum... not keen with the results, going to lower the cuttoff

DCLDE2013Cluster.cutoff =.05
DCLDE2013Cluster.chains =[];
updateClusterID(DCLDE2013Cluster)

%% Plotting

[~, k1]= find(DCLDE2013Cluster.localize_struct.hyd(5).score(7,:) < -49.5745e-6);
[~, k3]= find(sum(~isnan(DCLDE2013Cluster.localize_struct.hyd(5).score)) > 1);
k2 = intersect(k1, k3);


title1 =[ 'Localized Data Clusterd Using Ladder Structure and Mean'...
    newline 'Correlation Between ', num2str(NChidlHyd),' Hydrophone Pairs'];
title2 = 'Localized Data Plotted Against Time of Detection';


x = DCLDE2013Cluster.arrivalTable.Loc(k2,1);
y = DCLDE2013Cluster.arrivalTable.Loc(k2,2);

na_idx =find(~isnan(x))


chain = DCLDE2013Cluster.chains;
Cluster_id = DCLDE2013Cluster.Cluster_id;

color_var_drift= zeros(size(x,1),1);
color_var_drift_clustid = ones(size(x,1),1);
ColorVals = jet(length(k2));
colorscale2 = (DCLDE2013Cluster.localize_struct.hyd(5).rtimes(k2)- ...
    min(DCLDE2013Cluster.localize_struct.hyd(5).rtimes(k2)))/fs;



% Create color variable
cluster_n = 1; % initial cluster color (other than 0)

% Color all localizations in the cluster the same
for ii=1:length(chain)
    
    if chain(ii).n >1
        color_var_drift([ii, chain(ii).index]) = (log10(chain(ii).n))+.001;
        color_var_drift_clustid([ii, chain(ii).index]) = cluster_n;
        cluster_n =cluster_n+1;
    end
end


figure
fig = get(groot,'CurrentFigure');


x_trim_drift=x;
y_trim_drift=y;
color_var_drift = color_var_drift(k2);
color_var_drift_clustid =color_var_drift_clustid(k2)


cbar1label = 'log10(Number of Detections in Cluster)';
cbar2label = 'Elapsed Time Since First Detection (sec)';

% text labels, none for unknown data
txt_str2 =[''];
x_txt = max(x)-.1
y_txt = min(y) +.01


txt_str1 =[num2str(length(unique(Cluster_id))), ' Groups Predicted'];
y_txt = min(y);

% Randomly change values of the color ids for easer plotting
clustids = unique(color_var_drift_clustid);
new_ids = datasample([1:length(unique(color_var_drift_clustid))], ...
    length(unique(color_var_drift_clustid)), 'Replace', false); 
color_var_drift_clustid_shuff =color_var_drift_clustid;

for ii = 1:length(unique(color_var_drift_clustid))

color_var_drift_clustid_shuff(color_var_drift_clustid == clustids(ii)) = new_ids(ii);
end


%Plot data
figure
colormap((plasma))
scatter(x_trim_drift, y_trim_drift, 40, color_var_drift_clustid, 'filled');
title(title1)
ylabel('Latitude'); xlabel('Longitude')
h = colorbar;
ylabel(h, cbar1label)
hold on

scatter( dclde_2013_meta(:,12), dclde_2013_meta(:,11), 140, 'k', 'filled', 'd')
scatter(dclde_2013_meta(5,12), dclde_2013_meta(5,11), 140, 'r', 'filled', 'd')
hold off

figure
hold on
colormap((viridis))
scatter(x_trim_drift, y_trim_drift, 40, colorscale2, 'filled');
title(title2)
text(x_txt, y_txt, txt_str2)
ylabel('Latitude'); xlabel('Longitude')
scatter( dclde_2013_meta(:,12), dclde_2013_meta(:,11), 140, 'k', 'filled', 'd')
scatter(dclde_2013_meta(5,12), dclde_2013_meta(5,11), 140, 'r', 'filled', 'd')
h = colorbar;
ylabel(h, cbar2label)





