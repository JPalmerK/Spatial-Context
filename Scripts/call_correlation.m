close all; clear all; clc

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
LSQ_plot_thres =.001;
cor_thresh = -.4;
time_cut =12; % minutes
NChidlHyd = 1; % Number of pairs of hydrophones


[chain, Cluster_id] = MultiplHydrophonCluster(array_struct,localize_struct,...
    array_idx,fs, NChidlHyd, time_cut, cor_thresh, LSQ_plot_thres);






%% Cluster PMRF using ladders
n_corr = 4;
array_id=19;
Corr_coef_map = make_LSQ_corr_map(array_struct,...
    localize_struct, array_id, n_corr);

% correlate out to 15 min (dive time of a whale)
time_cut = 15*60;

% Set up initial times
x=localize_struct.hyd(array_struct(1).master).coordinates(2,2,:);
y=localize_struct.hyd(array_struct(1).master).coordinates(2,1,:);
times=localize_struct.hyd(array_struct(1).master).rtimes;
score=localize_struct.hyd(array_struct(1).master).score(2,:);
% calltype=localize_struct.hyd(array_struct(j).master).calltype;

x=squeeze(x);
y=squeeze(y);
score=squeeze(score);


for jj=1:size(Corr_coef_map,1)
    
    % Creat the chain links
    chain = make_Ladder_linkeges(squeeze(Corr_coef_map(jj,:,:)),...
        localize_struct, 19, 6000, 0.01, time_cut);

    
    %     % Create a histogram for the number of chains for the number of
    %     % hydrophone pairs
    %     n_per_chain = zeros(1, length(chain));
    %     for ii = 1:length(chain)
    %         n_per_chain(ii) = chain(ii).n;
    %     end
    %
    %     [~,edges] = histcounts(log10(n_per_chain),20);
    %     figure(2)
    %     subplot(2,2, jj)
    %     histogram(n_per_chain,10.^edges)
    %     set(gca,'xscale','log')
    %     title([ 'Clustered using ' num2str(jj) ' hydrophone Pairs'])
    %
    
    % For localized calls, plot according to group colo
    % Create color variable
    color_var= zeros(size(x));
    cluster_n = 1; % initial cluster color (other than 0)
    
    % Color all localizations in the cluster the same
    for ii=1:length(chain)
        
        if chain(ii).n >1
            color_var([ii, chain(ii).index]) = cluster_n;
            cluster_n =cluster_n+1;
        end
    end
    
    % trim the localizations
    
    [k1 k2]= find(score <.01);
    x_trim=x(k2);
    y_trim=y(k2);
    color_var = color_var(k2);
    
    figure(4)
    subplot(2,2, jj)
    scatter(x_trim, y_trim, 20, color_var, 'filled');
    title([ 'Clustered Using ' num2str(jj) ' Hydrophone Pairs'])
    
    
end

%% Add clock drift to one of the instruments

localize_struct_drift = localize_struct;

% Add one second of clock drift to the first slave instrument
localize_struct_drift.hyd(19).delays(:,1) = localize_struct_drift.hyd(19).delays(:,1) + 5;

% Run the analysis for correlation etc

n_corr = 2;
array_id=19;
Corr_coef_map_drift = make_LSQ_corr_map(array_struct,...
    localize_struct_drift, array_id, n_corr);




% Set up initial times
x=localize_struct.hyd(array_struct(1).master).coordinates(2,2,:);
y=localize_struct.hyd(array_struct(1).master).coordinates(2,1,:);
times=localize_struct.hyd(array_struct(1).master).rtimes;
score=localize_struct.hyd(array_struct(1).master).score(2,:);
% calltype=localize_struct.hyd(array_struct(j).master).calltype;

x=squeeze(x);
y=squeeze(y);
score=squeeze(score);


for jj=1:size(Corr_coef_map,1)
    
    % Creat the chain links
    chain_drift = make_Ladder_linkeges(squeeze(Corr_coef_map_drift(jj,:,:)),...
        localize_struct, 19, 0.9, time_cut);
    
    
    %     % Create a histogram for the number of chains for the number of
    %     % hydrophone pairs
    %     n_per_chain = zeros(1, length(chain_drift));
    %     for ii = 1:length(chain_drift)
    %         n_per_chain(ii) = chain_drift(ii).n;
    %     end
    %
    %     [~,edges] = histcounts(log10(n_per_chain),20);
    %     figure(2)
    %     subplot(2,2, jj)
    %     histogram(n_per_chain,10.^edges)
    %     set(gca,'xscale','log')
    %     title([ 'Clustered using ' num2str(jj) ' hydrophone Pairs'])
    
    
    % For localized calls, plot according to group colo
    % Create color variable
    color_var_drift= zeros(size(x));
    cluster_n = 1; % initial cluster color (other than 0)
    
    % Color all localizations in the cluster the same
    for ii=1:length(chain_drift)
        
        if chain_drift(ii).n >1
            color_var_drift([ii, chain_drift(ii).index]) = cluster_n;
            cluster_n =cluster_n+1;
        end
    end
    
    % trim the localizations
    
    [k1 k2]= find(score <.01);
    x_trim_drift=x(k2);
    y_trim_drift=y(k2);
    color_var_drift = color_var_drift(k2);
    
    figure(4)
    subplot(2,2, jj+2)
    scatter(x_trim_drift, y_trim_drift, 20, color_var_drift, 'filled');
    title([ '1s Drift Clustered Using ' num2str(jj) ' Hydrophone Pairs'])
    
    
end



%% Add 1s RANDOM clock drift to one of the instruments

localize_struct_drift = localize_struct;

% Add random clock drift (+-1 sec) to the slave instruments
for ii=1:4
    random_jitter = rand(size(localize_struct_drift.hyd(19).delays(:,1)));
    localize_struct_drift.hyd(19).delays(:,1) = ...
        localize_struct_drift.hyd(19).delays(:,ii) + random_jitter;
end
% Run the analysis for correlation etc

n_corr = 2;
array_id=19;
Corr_coef_map_drift = make_LSQ_corr_map(array_struct,...
    localize_struct_drift, array_id, n_corr);

time_cut = 500;



% Set up initial times
x=localize_struct.hyd(array_struct(1).master).coordinates(2,2,:);
y=localize_struct.hyd(array_struct(1).master).coordinates(2,1,:);
times=localize_struct.hyd(array_struct(1).master).rtimes;
score=localize_struct.hyd(array_struct(1).master).score(2,:);
% calltype=localize_struct.hyd(array_struct(j).master).calltype;

x=squeeze(x);
y=squeeze(y);
score=squeeze(score);


for jj=1:size(Corr_coef_map,1)
    
    % Creat the chain links
    chain_drift = make_Ladder_linkeges(squeeze(Corr_coef_map_drift(jj,:,:)),...
        localize_struct, 19, 0.9, time_cut);
    
    %
    %     Create a histogram for the number of chains for the number of
    %     hydrophone pairs
    %     n_per_chain = zeros(1, length(chain_drift));
    %     for ii = 1:length(chain_drift)
    %         n_per_chain(ii) = chain_drift(ii).n;
    %     end
    %
    %     [~,edges] = histcounts(log10(n_per_chain),20);
    %     figure(2)
    %     subplot(2,2, jj)
    %     histogram(n_per_chain,10.^edges)
    %     set(gca,'xscale','log')
    %     title([ 'Clustered using ' num2str(jj) ' hydrophone Pairs'])
    
    
    % For localized calls, plot according to group colo
    % Create color variable
    color_var_drift= zeros(size(x));
    cluster_n = 1; % initial cluster color (other than 0)
    
    % Color all localizations in the cluster the same
    for ii=1:length(chain_drift)
        
        if chain_drift(ii).n >1
            color_var_drift([ii, chain_drift(ii).index]) = cluster_n;
            cluster_n =cluster_n+1;
        end
    end
    
    % trim the localizations
    
    [k1 k2]= find(score <.1);
    x_trim_drift=x(k2);
    y_trim_drift=y(k2);
    color_var_drift = color_var_drift(k2);
    
    figure(4)
    subplot(2,2, jj+2)
    scatter(x_trim_drift, y_trim_drift, 20, color_var_drift, 'filled');
    title([ 'Random 1s Drift, Clustered Using ' num2str(jj) ' Hydrophone Pairs'])
    
    
end



%% Add linear clock drift to one of the instruments

% Linear clock drift up to 5min per 6 months (average from DPS
% conversation)

localize_struct_linear_drift = localize_struct;

% Model coefficients (minutes per 6 months, convert to seconds per day)
lm_slope = [0 ,1 ,3, 5]*60/(6*30*24*60*60);


% assume clock starts 1 hr before first cals
fs =6000;
time = localize_struct.hyd(19).rtimes/fs;% seconds now yes?
time = time - time(1) + (60*60);

time_add = [];

for ii=1:4
    % amount of clock drift for each units and time to add to the delays
    time_add(:,ii) = time*lm_slope(ii);
    
    % Add the clock drift to the arrival times
    localize_struct_linear_drift.hyd(19).delays(:,ii) = ...
        localize_struct_linear_drift.hyd(19).delays(:,ii) + time_add(:,ii);
    
end




n_corr = 2;
array_id=19;
Corr_coef_map_drift = make_LSQ_corr_map(array_struct,...
    localize_struct_linear_drift, array_id, n_corr);

time_cut = 500;


% Make another figure with the correlations but trim the lousy LSQ values
score=localize_struct.hyd(array_struct(1).master).score(2,:);
[k1 k2]= find(score <.01);
for ii = 1:size(Corr_coef_map_drift,1)
    
    % Create the correlation heat map
    figure(1)
    subplot(2,ceil(n_corr/2), ii)
    imagesc(squeeze(Corr_coef_map_drift(ii,k2,k2)));
    colormap jet
    
    if ii == 1
        title([num2str(ii) ' pair of hydropones'])
    else
        title([num2str(ii) ' pairs of hydropones'])
    end
    
end




% Set up initial times
x=localize_struct.hyd(array_struct(1).master).coordinates(2,2,:);
y=localize_struct.hyd(array_struct(1).master).coordinates(2,1,:);
times=localize_struct.hyd(array_struct(1).master).rtimes;


x=squeeze(x);
y=squeeze(y);
score=squeeze(score);


for jj=1:size(Corr_coef_map_drift,1)
    
    % Creat the chain links
    chain_drift = make_Ladder_linkeges(squeeze(Corr_coef_map_drift(jj,:,:)),...
        localize_struct, 19, 0.01, time_cut);
    
    %
    %     Create a histogram for the number of chains for the number of
    %     hydrophone pairs
    %     n_per_chain = zeros(1, length(chain_drift));
    %     for ii = 1:length(chain_drift)
    %         n_per_chain(ii) = chain_drift(ii).n;
    %     end
    %
    %     [~,edges] = histcounts(log10(n_per_chain),20);
    %     figure(2)
    %     subplot(2,2, jj)
    %     histogram(n_per_chain,10.^edges)
    %     set(gca,'xscale','log')
    %     title([ 'Clustered using ' num2str(jj) ' hydrophone Pairs'])
    
    
    % For localized calls, plot according to group colo
    % Create color variable
    color_var_drift= zeros(size(x));
    cluster_n = 1; % initial cluster color (other than 0)
    
    % Color all localizations in the cluster the same
    for ii=1:length(chain_drift)
        
        if chain_drift(ii).n >1
            color_var_drift([ii, chain_drift(ii).index]) = cluster_n;
            cluster_n =cluster_n+1;
        end
    end
    
    % trim the localizations for plotting
    x_trim_drift=x(k2);
    y_trim_drift=y(k2);
    color_var_drift = color_var_drift(k2);
    
    figure(4)
    subplot(2,2, jj+2)
    scatter(x_trim_drift, y_trim_drift, 20, color_var_drift, 'filled');
    title([ 'Clock Drift 0-5 min per 6months ' num2str(jj) ' Hydrophone Pairs'])
    colormap 'jet'
    
    
end



% Make a fgiure with time as coloring agent
figure(5)
scatter(x_trim_drift, y_trim_drift, 20, log10(time(k2)), 'filled');
title([ 'Clock Drift 0-5 min per 6months ' num2str(jj) ' Hydrophone Pairs'])
colormap 'jet'


%% Load the DCLDE2013 5day  dataset and give it a go

clear all; close all; clc

cd('/home/kpalmer/AnacondaProjects/Localisation/Scripts')

% Load the meta data for array plotting
dclde_2013_meta = xlsread('/cache/kpalmer/quick_ssd/data/DCLDE_2013_10Channel/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx');

% Create labels for plotting locations (undocumented function)
%labels=sprintfc('%d',1:10)
%
% plot(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 'o')
% text(dclde_2013_meta(:,12), dclde_2013_meta(:,11),labels,'VerticalAlignment','top','HorizontalAlignment','left')

load('DCLDE2013_RW_localizationsDCLDE_2013_10_Chan_chan_5_1_671.mat')
array_struct = localize_struct.hyd(5).array_struct;

% Trim the lat/lon space for efficeincy
array_struct.latgrid = array_struct.latgrid(1:15:length(array_struct.latgrid));
array_struct.longrid = array_struct.longrid(1:10:length(array_struct.longrid));

for ii=1:(length(array_struct.toa_diff)-1)
    
    mm = cell2mat(array_struct.toa_diff(ii+1));
    
    mm = mm(1:15:(size(mm,1)),1:10:(size(mm,2)));
    
    array_struct.toa_diff(ii+1) = {mm};
    clear mm
    
end

array_id=5; % center hydrophone

% trim the scores
score=localize_struct.hyd(array_id).score(5,:);
[~, k2]= find(score <.02);

% Trim calls
localize_struct.hyd(array_id).score = localize_struct.hyd(array_id).score(:,k2);

% Trim times
localize_struct.hyd(array_id).rtimes = localize_struct.hyd(array_id).rtimes(:,k2);

% Trim corrdinates
localize_struct.hyd(array_id).coordinates = localize_struct.hyd(array_id).coordinates(:,:,k2);

% trim dex
localize_struct.hyd(array_id).dex = localize_struct.hyd(array_id).dex(k2);

% trim coord time
localize_struct.hyd(array_id).coord_time = localize_struct.hyd(array_id).coord_time(k2,:);

% trim cross correlation score
localize_struct.hyd(array_id).cross_score = localize_struct.hyd(array_id).cross_score(k2,:);

% and delays (not dealing with CC matrix atm)
localize_struct.hyd(array_id).delays = localize_struct.hyd(array_id).delays(k2,:);


% Trim lat/lon

% Pick center array and number of hydrophone pairs to include
n_corr = 3;

% how close can calls be in time before they get new clusters in seconds
time_cut = 60*5;

% which hydrophones to include in the sub array, must be same lenght as
% n_corr

hyd_saves = [2,3,8];
%%

array_idx = 1
fs=2000;
time_cut = 20;
NChidlHyd=2;
cor_thresh=-1;
LSQ_plot_thres =1;

[chain, Cluster_id] = MultiplHydrophonCluster(array_struct,localize_struct,...
    array_idx,fs, NChidlHyd, time_cut, cor_thresh, LSQ_plot_thres);




%%

fs =2000;

% Set up initial times
x=localize_struct.hyd(array_struct(1).master).coordinates(5,2,:);
y=localize_struct.hyd(array_struct(1).master).coordinates(5,1,:);
times=localize_struct.hyd(array_struct(1).master).rtimes;


x=squeeze(x);
y=squeeze(y);
score=squeeze(localize_struct.hyd(array_id).score(5,:));

% more trimming for plot locations
[~, k2]= find(score <.01);

chain=struct();

for jj=1:size(Corr_coef_map,1)
    
    % Creat the chain links
    temp_chain = make_Ladder_linkeges(squeeze(Corr_coef_map(jj,:,:)),...
        localize_struct, array_id, fs, 0.2, time_cut);
    chain(jj).cc_dat = temp_chain;
    
    % Get the sort order of the chains by number of calls in each chain
    [~, chain_order] = sort([chain(jj).cc_dat(1:end).n],'descend');
    
    
    % For localized calls, plot according to group colo
    % Create color variable
    color_var_drift= zeros(size(x))/0;
    cluster_n = 1; % initial cluster color (other than 0)
    
    % Color all localizations in the cluster the same
    for ii=1:length(chain(jj).cc_dat)
        idx = chain_order(ii);
        if chain(jj).cc_dat(idx).n >=1
            color_var_drift([idx, chain(jj).cc_dat(idx).index]) = cluster_n;
            cluster_n =cluster_n+1;
        end
    end
    
    % trim the localizations for plotting
    x_trim_drift=x(k2);
    y_trim_drift=y(k2);
    color_var_drift = color_var_drift(k2);
    color_var_drift = color_var_drift - min(color_var_drift);
    
    figure(2)
    subplot(2,2, jj+1)
    scatter(x_trim_drift, y_trim_drift,...
        20, color_var_drift, 'filled');
    title(['Clusterd using '  num2str(jj) ' Hydrophone Pair'])
    colormap 'lines'
    hold on
    scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80,...
        'k', 'filled', 'd')
    
    clear temp_chain
    
end
%%

subplot(2,1, 1)
scatter(x_trim_drift, y_trim_drift, 20, times, 'filled');
title([ 'Elapsed Time' num2str(jj) ' (Reality)'])
colormap 'jet'
hold on
scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80, 'k', ...
    'filled', 'd')


% Make a fgiure with time as coloring agent
figure(6)
scatter(x_trim_drift, y_trim_drift, 20, times, 'filled');
title([ 'Clock Drift 0-5 min per 6months ' num2str(jj) ' Hydrophone Pairs'])
colormap 'jet'









%% Make a big datafram with all the relavent information for plotting


rng(1)

x=squeeze(localize_struct.hyd(array_struct(1).master).coordinates(5,2,:));
y=squeeze(localize_struct.hyd(array_struct(1).master).coordinates(5,1,:));
times=   localize_struct.hyd(array_struct(1).master).rtimes;
scores = localize_struct.hyd(array_struct(1).master).score(end,:);

% Create an empty structure for the clusters
all_dat_cluster = struct();


for nn =1:length(chain)

% create matrix with event id,x, y, time, score, and cluster id 
all_dat_cluster(nn).chain = zeros(sum([chain(1).cc_dat.n])+1, 4);
all_dat_cluster(nn).chain(:,1) = 0/0;
all_dat_cluster(nn).chain(:,2) = x;
all_dat_cluster(nn).chain(:,3) = y;
all_dat_cluster(nn).chain(:,4) = times'/fs; % in seconds
all_dat_cluster(nn).chain(:,5) = score;

% Assign each call to an ID and a color
for ii =1:length(chain(nn).cc_dat)

     idx = chain(nn).cc_dat(ii).index;
     all_dat_cluster(nn).chain(idx, 1) =ii;

end
all_dat_cluster(nn).chain(idx, end) = ii+1;

% Sort the data by time
[~,idx] = sort(all_dat_cluster(nn).chain(:,4));
all_dat_cluster(nn).chain = all_dat_cluster(nn).chain(idx,:);

% Create a colorscheme to plot with that will maintain the same colors
group_colors = jet(length(chain(nn).cc_dat));
all_dat_cluster(nn).group_colors_shuffle = group_colors(randperm(size(group_colors,1)),:);

% Now make an animated map that steps through time and displays the groups
% together

end






% define timescale in seconds (7 days continuous in DCLDE 32013)
tt = [0:(10*60):(7*24*60*60)];

% Durtion (in sampes) time calls will be allowed to stay on the plot (1hr)
expired_dur = 15*60;

group_colors_shuffle = group_colors(randperm(size(group_colors,1)),:);

% Extra threshold for plotting
thresh =.01;
fid = figure(4);
% writerObj = VideoWriter('AllLocs_lsqlessthan0p2.avi'); % Name it.
% writerObj.FrameRate = 60; % How many frames per second.
% open(writerObj);

for ii=1:length(tt)
    
    
    % get the calls within the timestamp
    call_idx = find(all_dat_cluster(1).chain(:,4)>=(tt(ii)-expired_dur)...
        & (all_dat_cluster(1).chain(:,4)<= tt(ii)));
    
    % If there are calls within the time window plot them
    if ~isempty(call_idx)
        cluster_ids=struct();
        % Get the names of the clusters included in the first 10 min
        for jj=1:length(all_dat_cluster)
            
            % get the unique cluster id's
            cluster_ids(jj).group = unique(all_dat_cluster(jj).chain(call_idx,1));
            
            % add all the indexes
            idxs =[];
            for kk=1:length(cluster_ids(jj).group)
                
                idxs = [idxs; find(all_dat_cluster(jj).chain(:,1) == cluster_ids(jj).group(kk))];
                
            end
            
            % IDXs is row index for the x, y, time and color
            cluster_ids(jj).idx = idxs;
        end
        
        
  
        
        
        % Find the indexes of all of the calls with the ids
        figure(4)
        
        % Plot the clusters
        for ll = 1:length(all_dat_cluster)
            
            % Raw index of the calls to plot
            plotted_call_idx = cluster_ids(ll).idx;
            
            % TRimmed based on threshold
            trimmed_idx = find(all_dat_cluster(ll).chain(plotted_call_idx,5)<=thresh)
            
            subplot(3,1,ll)
            % Index of the plotted calls
            scatter(all_dat_cluster(ll).chain(plotted_call_idx,2),...
                all_dat_cluster(ll).chain(plotted_call_idx,3), 20,...
                all_dat_cluster(ll).group_colors_shuffle(all_dat_cluster(ll).chain(plotted_call_idx,1),:),...
                'filled')
            hold on
            scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80, 'k', ...
                'filled', 'd')
            % add master array
            scatter(dclde_2013_meta(5,12), dclde_2013_meta(5,11), 80, 'r', ...
                'filled', 'd')
            % child arrays
%             scatter(dclde_2013_meta(child_hyds,12), dclde_2013_meta(child_hyds,11), 80, 'b', ...
%                 'filled', 'd')
            xlim([-70.6, -70])
            ylim([41.8, 42.5])  
            title(['Raw: Clustered Using ' num2str(ll) 'Hydrophone Pairs'])
            

            txt = {'Elapsed Time (min):', num2str(tt(ii)/60)};
            text(-70.4, 41.85,txt)

            hold off 
            pause(0.01)

%             % Index of the plotted calls
%             figure(5)
%             subplot(3,1,ll)
%             scatter(all_dat_cluster(ll).chain(trimmed_idx,2),...
%                 all_dat_cluster(ll).chain(trimmed_idx,3), 20,...
%                 all_dat_cluster(ll).group_colors_shuffle(all_dat_cluster(ll).chain(trimmed_idx,1),:),...
%                 'filled')
%             hold on
%             scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80, 'k', ...
%                 'filled', 'd')
%             % add master array
%             scatter(dclde_2013_meta(5,12), dclde_2013_meta(5,11), 80, 'r', ...
%                 'filled', 'd')
%             % child arrays
% %             scatter(dclde_2013_meta(child_hyds,12), dclde_2013_meta(child_hyds,11), 80, 'b', ...
% %                 'filled', 'd')
%             xlim([-70.6, -70])
%             ylim([41.8, 42.5])  
%             title(['Trimmed: Clustered Using ' num2str(ll) 'Hydrophone Pairs'])
%             hold off 
%            
           
%         frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%         writeVideo(writerObj, frame);
           
        end
        
    end
    
end
hold off
close(writerObj); % Saves the movie.





%% Make an animated map

% define timescale
tt = [0:(10*60):(7*24*60*60)] *2000;

% Durtion (in sampes) time calls will be allowed to stay on the plot (1hr)
expired_dur = 30*60*fs;


% just use the trimmed times
x=squeeze(localize_struct.hyd(array_struct(1).master).coordinates(5,2,:));
y=squeeze(localize_struct.hyd(array_struct(1).master).coordinates(5,1,:));
times=localize_struct.hyd(array_struct(1).master).rtimes;
x=x(k2);
y=y(k2);
times = times(k2);


color_map_time = jet(length(x));

color_map_group = struct();
color_map_group(1).colors = jet(length(chain(1).cc_dat));
color_map_group(1).color_map_group_shuffle = color_map_group(1).colors(randperm(size(color_map_group(1).colors,1)),:);
color_map_group(2).colors = jet(length(chain(2).cc_dat));
color_map_group(2).color_map_group_shuffle = color_map_group(2).colors(randperm(size(color_map_group(2).colors,1)),:);
color_map_group(3).colors = jet(length(chain(3).cc_dat));
color_map_group(3).color_map_group_shuffle = color_map_group(3).colors(randperm(size(color_map_group(3).colors,1)),:);



% which channels were used as children
child_hyds = localize_struct.hyd(array_id).array_struct.slave(1:n_corr);

for ii = 1:length(tt)
    
    % get the calls within the timestamp
    call_idx = find(times>=(tt(ii)-expired_dur) & (times<= tt(ii+1)));
    
    if ~isempty(call_idx)
        figure(1)
%         % Plot the values colored by time
%         subplot(2,1,1)
%         scatter(x(call_idx), y(call_idx), 20, color_map_time(call_idx,:), 'filled')
%         hold on
%         scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80, 'k', ...
%             'filled', 'd')
%         % add master array
%         scatter(dclde_2013_meta(5,12), dclde_2013_meta(5,11), 80, 'r', ...
%             'filled', 'd')
%         % child arrays
%         scatter(dclde_2013_meta(child_hyds,12), dclde_2013_meta(child_hyds,11), 80, 'y', ...
%             'filled', 'd')
%         title('Unclustered Calls')
%         xlim([-70.6, -70])
%         ylim([41.8, 42.5])
%         hold off
%         
        % Color by cluster
        figure(4)
        subplot(3,2,1)
        cluster_color = color_map_group(1).color_map_group_shuffle(call_idx)+1;
        scatter(x(call_idx), y(call_idx), 20, cluster_color, 'filled')
        hold on
        scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80, 'k', ...
            'filled', 'd')
        % add master array
        scatter(dclde_2013_meta(5,12), dclde_2013_meta(5,11), 80, 'r', ...
            'filled', 'd')
        % child arrays
        scatter(dclde_2013_meta(child_hyds,12), dclde_2013_meta(child_hyds,11), 80, 'b', ...
            'filled', 'd')
        xlim([-70.6, -70])
        ylim([41.8, 42.5])
        hold off
        
        
        
        
        subplot(3,2,3)
        cluster_color = color_var_drift(call_idx)+1;
        scatter(x(call_idx), y(call_idx), 20, color_map_group(cluster_color,:), 'filled')
        hold on
        scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80, 'k', ...
            'filled', 'd')
        % add master array
        scatter(dclde_2013_meta(5,12), dclde_2013_meta(5,11), 80, 'r', ...
            'filled', 'd')
        % child arrays
        scatter(dclde_2013_meta(child_hyds,12), dclde_2013_meta(child_hyds,11), 80, 'b', ...
            'filled', 'd')
        xlim([-70.6, -70])
        ylim([41.8, 42.5])
        hold off
        
        
        subplot(3,2,3)
        cluster_color = color_var_drift(call_idx)+1;
        scatter(x(call_idx), y(call_idx), 20, color_map_group(cluster_color,:), 'filled')
        hold on
        scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80, 'k', ...
            'filled', 'd')
        % add master array
        scatter(dclde_2013_meta(5,12), dclde_2013_meta(5,11), 80, 'r', ...
            'filled', 'd')
        % child arrays
        scatter(dclde_2013_meta(child_hyds,12), dclde_2013_meta(child_hyds,11), 80, 'b', ...
            'filled', 'd')
        xlim([-70.6, -70])
        ylim([41.8, 42.5])
        hold off       
        
        
        pause(.05)
    end
    
    
end









%% Plot by time but by cluster

group_ids = unique(color_var_drift);


for ii=1:(length(tt)-1)
    
    
    % get the calls within the timestamp
    call_idx = find(times>=(tt(ii)-expired_dur) & (times<= tt(ii+1)));
    
    % Get the cluster ids
    cluster_id = unique(color_var_drift(call_idx));
    
    % add the indexes of all the clusters to call_idx
    
    for jj=1:length(unique(cluster_id))
        
        % Find all values of the cluster
        new_idx = find(color_var_drift == cluster_id(jj));
        
        if length(new_idx)>1
            call_idx = [call_idx, new_idx(2:end)'];
            %disp(['added one ' num2str(ii)])
        end
        
    end
    
    % Color by cluster
    cluster_color = color_var_drift(call_idx)+1;
    
    scatter(x(call_idx), y(call_idx), 20, color_map_group_shuffle(cluster_color,:), 'filled')
    
    hold on
    scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80, 'k', ...
        'filled', 'd')
    % add master array
    scatter(dclde_2013_meta(5,12), dclde_2013_meta(5,11), 80, 'r', ...
        'filled', 'd')
    % child arrays
    scatter(dclde_2013_meta(child_hyds,12), dclde_2013_meta(child_hyds,11), 80, 'p', ...
        'filled', 'd')
    xlim([-70.6, -70])
    ylim([41.8, 42.5])
    hold off
    
    
end







%% Different animated plot that plots the groups as a function of time
group_id = unique(color_var_drift);

% Sort group ID's by time
[~, idx] = sort(times);



% Create a dummy plotted variable
been_plotted = zeros(size(times));

while sum(been_plotted)<length(been_plotted)
    
    % get the group of the first zero event and the ID of that event
    group_val = find(been_plotted == 0,1);
    group_id = color_var_drift(group_val);
    
    % get the gorup indexes
    group_idx = find(color_var_drift == group_id);
    
    % Get the times of the group
    time_vals = times(group_idx);
    
    % Get lat and lon
    x_lon = x(group_idx);
    y_lat = y(group_idx);
    
    
    figure(2)
    scatter(x_lon, y_lat, 20, (time_vals-min(time_vals))/fs/60/60, 'filled')
    title(['Cluster ' num2str(ii)])
    xlim([-70.6, -70])
    ylim([41.8, 42.5])
    
    c = colorbar;
    c.Label.String = 'Cluster Duration (min)';
    
    hold on
    scatter(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 80, 'k', ...
        'filled', 'd')
    hold off
    
    pause(max(.1, length(group_idx)/500))
    
    been_plotted(group_idx) =1;
    disp(num2str(sum(been_plotted)))
    
    
    
    
    
    
end













