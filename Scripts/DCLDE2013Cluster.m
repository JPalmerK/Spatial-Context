
%% Load the DCLDE2013 5day  dataset and give it a go

clear all; close all; clc

cd('/home/kpalmer/AnacondaProjects/Localisation/Scripts')

% Load the meta data for array plotting
dclde_2013_meta = xlsread(strcat('/cache/kpalmer/quick_ssd/data/',...
    'DCLDE_2013_10Channel/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx'));

% Convert the meta data to a structure in the format that the
% GPL/localisation code expects
hydrophone_struct= struct();
for ii=1:size(dclde_2013_meta,1)
    hydrophone_struct(ii).name = num2str(dclde_2013_meta(ii,1));
    hydrophone_struct(ii).location = dclde_2013_meta(ii,[11:12]);
    hydrophone_struct(ii).depth= abs(dclde_2013_meta(ii, 13));
    hydrophone_struct(ii).channel=ii;
end



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


%% Run the cluster analysis
% Pick center array and number of hydrophone pairs to include

array_idx = 1;
fs=2000;
time_cut = 40;
NChidlHyd= 1;
cor_thresh=.01;
LSQ_plot_thres =1;

[chain, Cluster_id] = MultiplHydrophonCluster(array_struct,localize_struct,...
    array_idx,fs, NChidlHyd, time_cut, cor_thresh, LSQ_plot_thres, hydrophone_struct);

%% Create figures for workshop


% 1) Hydropone TOA

% Create figure
fig = figure;
colormap(bone);
plot_val =[array_struct.toa_diff{1,2}];

axes1 = axes('Parent',fig);
hold(axes1,'on');
imagesc(array_struct.longrid, (array_struct.latgrid),plot_val)
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 40, 'k', 'filled', 'd')
scatter(hyd_arr([5,1],2), hyd_arr([5,1],1), 40, 'r', 'filled', 'd')
% Create ylabel
ylabel({'Lat'},'FontWeight','bold');
% Create xlabel
xlabel({'Lon'},'FontWeight','bold');
% Create title
title({'Arrival Time Differences (DCLDE 2013)'},'FontSize',14);
xlim(axes1,[-70.7074671052631 -69.9178242481203]);
ylim(axes1,[41.6642345183486 42.6807597477064]);
box(axes1,'on');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'FontWeight','bold','Layer','top');
% Create colorbar
colorbar('peer',axes1);

h = colorbar;
ylabel(h, 'Expected Time Difference of Arrival (s)')



for ii =1:4
    kk=37
    % Pick a set of call delays (number of calls)
    delays = localize_struct.hyd(ParentHyd).delays(kk,:);
    
    % Get the tdoa space
    toa_space = (cell2mat(array_struct(array_idx).toa_diff(ii+1)));
    
    % potential toa space
    LSQ_VAL = abs(toa_space - delays(ii));
    
    LSQ_space = LSQ_VAL;
    
    hold off
    % Create figure
    fig = figure;
    colormap(copper);
    
    axes1 = axes('Parent',fig);
    hold(axes1,'on');
    imagesc(array_struct.longrid, (array_struct.latgrid), (LSQ_space))
    hold on
    scatter(hyd_arr(:,2), hyd_arr(:,1), 60, 'w', 'filled', 'd')
    scatter(hyd_arr([5,ii],2), hyd_arr([5,ii],1), 60, 'r', 'filled', 'd')
    % Create ylabel
    ylabel({'Lat'},'FontWeight','bold');
    % Create xlabel
    xlabel({'Lon'},'FontWeight','bold');
    % Create title
    title(['Call 1 Hydrophones 5 and ', num2str(ChildHyds(ii))],'FontSize',14);
    xlim([-70.7074671052631 -69.9178242481203]);
    ylim([41.6642345183486 42.6807597477064]);
    box(axes1,'on');
    axis(axes1,'ij');
    % Set the remaining axes properties
    set(axes1,'FontWeight','bold','Layer','top');
    % Create colorbar
    colorbar('peer',axes1);
    
    h = colorbar;
    ylabel(h, 'TOA Space - Measured Delay')
    
end



for call_id =1:3
    
    
    
    tdoa_diffs = Px(1,call_id,:);
    arrival1 = reshape((tdoa_diffs), 201,201).';
    
    hold off
    % Create figure
    fig = figure;
    colormap(copper);
    
    axes1 = axes('Parent',fig);
    hold(axes1,'on');
    imagesc(array_struct.longrid, array_struct.latgrid, arrival1)
    hold on
    scatter(hyd_arr(:,2), hyd_arr(:,1), 60, 'w', 'filled', 'd')
    scatter(hyd_arr([5,1],2), hyd_arr([5,1],1), 60, 'r', 'filled', 'd')
    %  scatter(localize_struct.hyd(5).coordinates(end, 2, call_id),...
    %      localize_struct.hyd(5).coordinates(end, 1, call_id),40, 'b', 'filled')
    % Create ylabel
    ylabel({'Lat'},'FontWeight','bold');
    % Create xlabel
    xlabel({'Lon'},'FontWeight','bold');
    % Create title
    title({['Call ',num2str(call_id),' TDOA Space Hydrophones 1 and 5']},'FontSize',14);
    xlim([-70.7074671052631 -69.9178242481203]);
    ylim([41.6642345183486 42.6807597477064]);
    box(axes1,'on');
    axis(axes1,'ij');
    % Set the remaining axes properties
    set(axes1,'FontWeight','bold','Layer','top');
    % Create colorbar
    colorbar('peer',axes1);
    
    h = colorbar;
    ylabel(h, 'TOA Space - Measured Delay')
end





% Mean for call 1
sub = Px(:,37,:);
sub = 1./-(1-(1+ exp(reshape(squeeze(mean(sub,1)), 201,201).')))



hold off
% Create figure
fig = figure;
colormap(bone);

axes1 = axes('Parent',fig);
hold(axes1,'on');
imagesc(array_struct.longrid, array_struct.latgrid, sub)
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 60, 'w', 'filled', 'd')
%  scatter(localize_struct.hyd(5).coordinates(end, 2, call_id),...
%      localize_struct.hyd(5).coordinates(end, 1, call_id),40, 'b', 'filled')
% Create ylabel
ylabel({'Lat'},'FontWeight','bold');
% Create xlabel
xlabel({'Lon'},'FontWeight','bold');
% Create title
title({['Call ',num2str(1),' TDOA Space Hydrophones 1 and 5']},'FontSize',14);
xlim([-70.7074671052631 -69.9178242481203]);
ylim([41.6642345183486 42.6807597477064]);
box(axes1,'on');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'FontWeight','bold','Layer','top');
% Create colorbar
colorbar('peer',axes1);

h = colorbar;
ylabel(h, 'TOA Space - Measured Delay')




















