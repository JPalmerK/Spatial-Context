
close all; clear all; clc

cd('/home/kpalmer/AnacondaProjects/Localisation/Scripts')


% Load metadata to create Hydrophone Structur
dclde_2013_meta = xlsread(strcat('/cache/kpalmer/quick_ssd/data/',...
    'DCLDE_2013_10Channel/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx'));

% Load the DCLDE array structure (just for the array structure, depth, and SSP)
load('DCLDE2013_RW_localizationsDCLDE_2013_10_Chan_chan_5_1_671.mat')



% Parent hydrophone for the analysis
parent =5;
fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;
array_struct = localize_struct.hyd(parent).array_struct;

% Trim the lat/lon space for efficeincy (oops!)
array_struct.latgrid = array_struct.latgrid(1:15:length(array_struct.latgrid));
array_struct.longrid = array_struct.longrid(1:10:length(array_struct.longrid));

% Dummy for cleaning grids
for ii=1:(length(array_struct.toa_diff)-1)
    
    mm = cell2mat(array_struct.toa_diff(ii+1));
    
    mm = mm(1:15:(size(mm,1)),1:10:(size(mm,2)));
    
    array_struct.toa_diff(ii+1) = {mm};
    clear mm
    
end



% Convert the meta data to a structure in the format that the
% GPL/localisation code expects
hydrophone_struct= struct();
for ii=1:size(dclde_2013_meta,1)
    hydrophone_struct(ii).name = num2str(dclde_2013_meta(ii,1));
    hydrophone_struct(ii).location = dclde_2013_meta(ii,[11:12]);
    hydrophone_struct(ii).depth= abs(dclde_2013_meta(ii, 13));
    hydrophone_struct(ii).channel=ii;
end
hyd_arr = struct2cell(hydrophone_struct);
hyd_arr =vertcat(hyd_arr{2,:,:});




%% Create the moving agents in the enviornment

% Each simulation has an overarching simulation parameter set and 
% each agent in the simulation has two parameter sets: calling and movement. 
% The calling parameter set describes the agent's calling pattern. 
% Calling options include consistant and random. Options include constant, random, and
% bout*. Also provide the call depth (MUST MATCH LSQ simulation depth).
% Soundspeed (here a single value)

% Movement parameters dictate how the agents move in the envornment and
% include, model, max_speed, initial lat/long, bearing (for directed travel and
% constant). Model options included random walk/search, directed travel, and
% linear

clear spaceWhale

% Simulation parameters that apply to the whole system
n_hrs = 5;
param_sim.dur= 60*60*n_hrs; % duration of the entire simulation
param_sim.fs =2000; % not currently used but will  be
param_sim.ssp = ssp(1,2);
param_sim.depth = grid_depth;

% add the simulation parameters to the spaceWhale
spaceWhale.param_sim = param_sim;

% Set whale start point in meters (ll corner of array)
parm_movement.max_speed = 8.4;
parm_movement.duration = 60*60*n_hrs; % currently needs to be the same for all

% Make call parameters
parm_calling.model = 'bout';


% Create 10 agents with different movement behaviours
for ii=1:7
    parm_calling.frequency = 400*rand; % calling every 200 seconds on average
    
    % Movement parameters for each  agent(s)
    parm_movement.lat_init =  min(hyd_arr(:,1)) + (max(hyd_arr(:,1)) - min(hyd_arr(:,1)))*rand(1);
    parm_movement.lon_init =  min(hyd_arr(:,2)) + (max(hyd_arr(:,2)) - min(hyd_arr(:,2)))*rand(1);
    parm_movement.bearing = 180 -180*rand(1); %degrees
    
    
    % Duration - random number of hours (uniform) times some portion of
    % that (rand)
    parm_movement.duration = floor(60*60*(randi(n_hrs)*rand(1))); % how long is the agent around
    
    
    parm_movement.start_time = randi([0, param_sim.dur-parm_movement.duration],1)
    
    
%     % Set last animal to random walk
     if mod(ii,2) ==1
% %         parm_movement.model = 'linear';
%     else
        parm_movement.model = 'directed travel';
%         
%     elseif ii==3
     else
         parm_movement.model = 'search'
     end
    % Assume max swimming speed of 10mph http://www.listenforwhales.org/page.aspx?pid=451
    
    
    parm_calling.start_time = 1;
    spaceWhale.agent(ii).parm_movement = parm_movement;
    spaceWhale.agent(ii).parm_calling = parm_calling;

end


% add the movement
[spaceWhale] = AgentMovement(spaceWhale, array_struct, hydrophone_struct);


%% Plot the agents and the hydrophones

figure(1) 

examp.arrivalTable.ArrivalSec(:,1)

% Make color scale
ColorVals = lines(length(spaceWhale.agent));
TimeColorVals = parula(spaceWhale.param_sim.dur);

for ii=1:length(spaceWhale.agent)
    
    subplot(1,2,2)
    hold on 
    plot(spaceWhale.agent(ii).location(:,2),...
        spaceWhale.agent(ii).location(:,1),...
        'Color', [ColorVals(ii,:)])

    % Add points indicating when the agent produced a call
    call_times = spaceWhale.agent(ii).call_times;
    scatter(spaceWhale.agent(ii).location(call_times,2),...
        spaceWhale.agent(ii).location(call_times,1),40,...
       [ColorVals(ii,:)], 'filled')
   
   
   subplot(1,2,1)
    hold on
    % Add points indicating when the agent produced a call
    call_times = spaceWhale.agent(ii).call_times;
    scatter(spaceWhale.agent(ii).location(call_times,2),...
        spaceWhale.agent(ii).location(call_times,1),40,...
        [TimeColorVals(call_times,:)], 'filled')

   
   

end
subplot(1,2,1)
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')

subplot(1,2,2)
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')


% Create ylabel
ylabel({'Lat'},'FontWeight','bold');
% Create xlabel
xlabel({'Lon'},'FontWeight','bold');
% Create title
title({'Simulated Movement and Call Times'},'FontSize',14);

%% Transform spaceWhales into localization structure 
% for handing of to MultiplHydrophonCluster


% Clean out locations/times that are outside the array(ish) boundaries
lat_lim = [42 42.35];
lon_lim = [-70.55 -70.1];

% Combine TOA (r times on the parent) and TDOA delays from space whales
rtimes = [];
agentID =[]; 
TDOAs = [];
coords =[];

for ii  = 1: length(spaceWhale.agent)
    
    %Call Time Locs
    call_timeloc = spaceWhale.agent(ii).location(spaceWhale.agent(ii).call_times,:);
    
    in_array_idx = find(call_timeloc(:,1)>=lat_lim(1) & call_timeloc(:,1)<=lat_lim(2) & ...
        call_timeloc(:,2)>=lon_lim(1) & call_timeloc(:,2)<=lon_lim(2));
    
    % lat long coordinates at times where there was a call
    coords = [coords; call_timeloc(in_array_idx,:)];
    

    % Get the arrival times of the calls at the parent
    rtimes = [rtimes spaceWhale.agent(ii).Arrival_times(in_array_idx,parent)'];
    
    % get the TDO and out the parent hydrophone column
    tdoa = struct2array(spaceWhale.agent(ii).tdoa(parent));  
    tdoa(:,parent) = [];
    TDOAs = [TDOAs; tdoa(in_array_idx,:)];
    
    % Retain agent ID for validation
    agentID = [agentID repmat(ii,1, length(tdoa))];
    

    
end
clear tdoa

% Get hte arrival time order 
[~, aTimeOrder] =  sort(rtimes);

localize_struct_sim = struct();
localize_struct_sim.hyd(parent).cc_matrix = [];
localize_struct_sim.hyd(parent).delays = TDOAs(aTimeOrder,:);
localize_struct_sim.hyd(parent).rtimes = rtimes(aTimeOrder)*fs; % arrival times in samples
localize_struct_sim.hyd(parent).coordinates= coords(aTimeOrder,:);
localize_struct_sim.hyd(parent).agent_ids = agentID(aTimeOrder);
localize_struct_sim.hyd(parent).score = ones(size(agentID));

%% Create the clusters 

close all
array_idx = 1;
threshs = [.01 0.03 0.05];
a = zeros(length(threshs),3);
b = zeros(length(threshs),3);
c = zeros(length(threshs),3);
d = zeros(length(threshs),3);

for jj = 1:1
    for ii = 1:1

    NChidlHyd = ii; % number of pairs of hydrphones
    time_cut = 10*60; % don't correlate anything greater than 5 min apart
    cor_thresh = threshs(jj);

    [chain, Cluster_id] = MultiplHydrophonCluster(array_struct,...
        localize_struct_sim, array_idx,fs, NChidlHyd, time_cut,...
        cor_thresh, [], hydrophone_struct);

    Cluster_id = Cluster_id+1;
   

    newClusterId = allignclusters(localize_struct_sim.hyd(parent).agent_ids, Cluster_id);
    [e,f,g,h] = RandIndex(newClusterId, localize_struct_sim.hyd(parent).agent_ids');

    a(jj, ii) = e;
    b(jj, ii) = f;
    c(jj, ii) = g;
    d(jj, ii) = h;

    end
end

 %% Create Confusion matrix 
 
 
% Stick the real id's next to the predicted clusters
truthAndpreds =[ localize_struct_sim.hyd(parent).agent_ids' newClusterId];
truthAndpreds(:,3)=0; % True class
truthAndpreds(:,4)=0; % Score
truthAndpreds(:,5)=0; % Prediction after applying likelihood ratio

% Simulate a binary upcall detector assume approximately 1/5th of the
% humpback calls are upcalls

% Likelyhood of the detector tagging a right whale upcall as such is
% approximately 0.9 so pulled from a normal distribution centered on 2.5
% and inv.logit transformed

% Likelihood of the detector tagging a humpback call as a right whale is
% 1/5  so one fifth of the simulated classifier thresholds should fall over
% 0.5.


% Assume any call with a classification threshold above 0.5 is tagged as RW
% and any detection classification less than 0.5 is tagged as a humpback


for ii =1:length(spaceWhale.agent)


    % If it's odd it's definitely a humpbacks
    if mod(ii,2)
    
        disp('Right Whale')
        
        % Get the indicies of the calls of that individual
        rw_idx = find(truthAndpreds(:,1) == ii);
        
        % Pull from distribution to estimate classificaiton probability 
        rw_mean = 1.5;
        rw_sd = 1;
        score =  rw_sd.*randn(length(rw_idx),1) + rw_mean;
        score = 1./(1+exp(-score));
        
        truthAndpreds(rw_idx,3) =1;
        truthAndpreds(rw_idx,4) = score'
        
        % Likelihood ratio
        LR = log(prod(score./(1-score)));
        
        truthAndpreds(rw_idx,5) = LR;
    
    
    else
        disp('Humpback')
        
        % Get the indicies of the calls of that individual
        mn_idx = find(truthAndpreds(:,1) == ii);
        
        % Pull from distribution to estimate classificaiton probability 
        mw_mean = -0.75;
        mn_sd = 0.5;
        score =  mn_sd.*randn(length(mn_idx),1) + mw_mean;
        score = 1./(1+exp(-score));
        
        truthAndpreds(mn_idx,3) = 0;
        truthAndpreds(mn_idx,4) = score';
    
        
        
        % Likelihood ratio
        LR = log(prod(score./(1-score)));
        
        truthAndpreds(mn_idx,5) = LR;
    
    end
    
    
    
end


% Create the confusion matrix for the unclustered data
tot_mn = sum(truthAndpreds(:,3) ==0);
tot_eg = sum(truthAndpreds(:,3) ==1);


% Proportion of right whale calls correctly classified as humpack
Tp_eg = sum(truthAndpreds(:,4) > 0.5 & truthAndpreds(:,3)==1);

% Proportion of right whale calls incorrectly tagged as humpbacks
Fn_eg = sum(truthAndpreds(:,4) < 0.5 & truthAndpreds(:,3)==1);


% Proportion of humpback calls correctly classified as humpack
Tp_mn = sum(truthAndpreds(:,4) < 0.5 & truthAndpreds(:,3)==0);

% Proportion of humpback callsincorrectly tagged as right whale
Fp_mn = sum(truthAndpreds(:,4) > 0.5 & truthAndpreds(:,3)==0);



% Create the confusion matrix for the clustered data

% Proportion of right whale calls correctly classified as humpack
Tp_egc = sum(truthAndpreds(:,5) > 5 & truthAndpreds(:,3)==1);

% Proportion of right whale calls incorrectly tagged as humpbacks
Fn_egc = sum(truthAndpreds(:,5) < (1/5) & truthAndpreds(:,3)==1);


% Proportion of humpback calls correctly classified as humpack
Tp_mnc = sum(truthAndpreds(:,5) < (1/5) & truthAndpreds(:,3)==0);

% Proportion of humpback callsincorrectly tagged as right whale
Fp_mnc = sum(truthAndpreds(:,5) > 5 & truthAndpreds(:,3)==0);


  
 
 x = -4:.02:4;
 y = 1./(1+exp(-x));
 
 plot(x,y)
 
 
 

 
 
 
 

 %%
% close all
% a1 = zeros(4,3);
% b1 = zeros(4,3);
% c1 = zeros(4,3);
% d1 = zeros(4,3);
% 
% for jj = 1:4
%     for ii = 1:3
% 
%     NChidlHyd = ii % number of pairs of hydrphones
%     time_cut = 5*60; % don't correlate anything greater than 5 min apart
%     cor_thresh = threshs(jj);
% 
%     [chain, Cluster_id] = MultiplHydrophonCluster(array_struct,...
%         localize_struct_sim, array_idx,fs, NChidlHyd, time_cut,...
%         cor_thresh, [], hydrophone_struct);
% 
%     disp([num2str(length(unique(Cluster_id))), ' Clusters'])
% 
% 
%     [e,f,g,h] = RandIndex(localize_struct_sim.hyd(parent).agent_ids',...
%         Cluster_id-1);
% 
%     a1(jj, ii) = e;
%     b1(jj, ii) = f;
%     c1(jj, ii) = g;
%     d1(jj, ii) = h;
% 
%     end
% end
% 






