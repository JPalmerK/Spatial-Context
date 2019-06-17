
close all; clear all; clc
%
% cd('/home/kpalmer/AnacondaProjects/Localisation/Scripts')
%
%
% % Load metadata to create Hydrophone Structur
% dclde_2013_meta = xlsread(strcat('/cache/kpalmer/quick_ssd/data/',...
%     'DCLDE_2013_10Channel/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx'));
% % Load the DCLDE array structure (just for the array structure, depth, and SSP)
% load('/home/kpalmer/AnacondaProjects/Localisation/Scripts/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_array_struct1_671.mat')
% load('/home/kpalmer/AnacondaProjects/Localisation/Scripts/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_localize_struct_671.mat')


% Load metadata to create Hydrophone Structur
dclde_2013_meta = xlsread('C:\Users\Kaitlin\Desktop\DCL2013_NEFSC_SBNMS_200903_metadata.xlsx');
load('D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_array_struct1_671.mat')
load('D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_localize_struct_671.mat')


% Parent hydrophone for the analysis
parent =5;
fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;



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


%% Create the agents in the enviornment

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
n_hrs = 8;
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
parm_calling.model = 'random';


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
    parm_movement.start_time = randi([0, param_sim.dur-parm_movement.duration],1);
    
    if mod(ii,2) ==1
        parm_movement.model = 'directed travel';

    else
        parm_movement.model = 'search';
    end
    % Assume max swimming speed of 10mph http://www.listenforwhales.org/page.aspx?pid=451
    
    
    parm_calling.start_time = 1;
    spaceWhale.agent(ii).parm_movement = parm_movement;
    spaceWhale.agent(ii).parm_calling = parm_calling;
    
end


% add the movement
[spaceWhale] = AgentMovement(spaceWhale, array_struct, hydrophone_struct);








%% Plot the agents 

% Make color scale
ColorVals = lines(length(spaceWhale.agent));

for ii=1:length(spaceWhale.agent)
    
    call_times = spaceWhale.agent(ii).call_times;
    
    %subplot(1,3,1)
    
    plot(spaceWhale.agent(ii).location(:,2),...
        spaceWhale.agent(ii).location(:,1),...
        'Color', [ColorVals(ii,:)])
    hold on
    scatter(spaceWhale.agent(ii).location(call_times,2),...
        spaceWhale.agent(ii).location(call_times,1),40,...
        [ColorVals(ii,:)], 'filled')
    
%     subplot(1,3,2)
%     scatter3(spaceWhale.agent(ii).location(call_times,2),...
%         spaceWhale.agent(ii).location(call_times,1),...
%         [0 diff(call_times)],...
%         40,[ColorVals(ii,:)], 'filled')
%     hold on
%     
%     
%     subplot(1,3,3)
%     hold on
%     plot(call_times/60, ones(size(call_times))*ii*.5,...
%         'Color', [ColorVals(ii,:)])
%     xlabel('Time in Simulation (min)')
end

subplot(1,3,1)
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')
xlim([min(hyd_arr(:,2))-.05  max(hyd_arr(:,2))+.05])
ylim([min(hyd_arr(:,1))-.05  max(hyd_arr(:,1))+.05])


subplot(1,3,2)
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')
xlim([min(hyd_arr(:,2))-.05  max(hyd_arr(:,2))+.05])
ylim([min(hyd_arr(:,1))-.05  max(hyd_arr(:,1))+.05])


% Create ylabel
ylabel({'Lat'},'FontWeight','bold');
% Create xlabel
xlabel({'Lon'},'FontWeight','bold');
% Create title
title({'Simulated Movement and Call Times'},'FontSize',14);
xlim([min(hyd_arr(:,2))-.05  max(hyd_arr(:,2))+.05])
ylim([min(hyd_arr(:,1))-.05  max(hyd_arr(:,1))+.05])


%% Create a table containing arrival times, ranges, and true ids

% Step through all agents and grab the call times, arrivals, range and
% call id

T = table(spaceWhale.agent(1).call_times', ...
    spaceWhale.agent(1).Arrival_times,...
    spaceWhale.agent(1).RangeKm,...
    ones(length(spaceWhale.agent(1).call_times),1),...
    'VariableNames',{'callTimes','ArrivalSec', 'RangeKm', 'ID'});

for ii = 2:length(spaceWhale.agent)
    Tnew =table(spaceWhale.agent(ii).call_times', ...
        spaceWhale.agent(ii).Arrival_times,...
        spaceWhale.agent(ii).RangeKm,...
        ones(length(spaceWhale.agent(ii).call_times),1)*ii,...
        'VariableNames',{'callTimes','ArrivalSec', 'RangeKm', 'ID'});
    T = [T ; Tnew];
end

% Clear out any detections greater than 8km from the receiver
distance_bool = ones(size(T.RangeKm));
distance_bool(find(T.RangeKm>8)) = nan;

T.RangeKm = T.RangeKm.*distance_bool;
T.ArrivalSec = T.ArrivalSec.*distance_bool;

NChidlHyd = 2;
array_idx = 1;
cor_thresh =.03;
time_cut = 60*15;

% Array containing the parent hydrophone and two selected child hydrophones
Arrivals2Chan = [T.ArrivalSec(:,[5,1,2]) T.ID];

% Remove arrivals where there are two or more Nans
Arrivals2Chan = Arrivals2Chan(sum(~isnan(Arrivals2Chan(:,1:3)),2)>1,:);

% Remove calls that weren't detected on the parent hydrophone
Arrivals2Chan = Arrivals2Chan(~isnan(Arrivals2Chan(:,1)),:);

% Sort by start time
[~,sortidx] = sort(Arrivals2Chan(:,1));
Arrivals2Chan = Arrivals2Chan(sortidx,:);

% Create a matrix that incorporates miss association
ArrivalsMiss = zeros(size(Arrivals2Chan));
ArrivalsMiss(:,1)  = Arrivals2Chan(:,1);
ArrivalsMiss(:,end)  = Arrivals2Chan(:,end);


for channum =1:NChidlHyd
    
    % Distance between two sensors - maximum exptected TDOA
    % depth between calling whale and array
    depth_range = hydrophone_struct(5).depth...
        - hydrophone_struct(channum).depth ;
    
    % Calculate the horizontal distance between the calling
    % locations and the hydrophone
    horizontal_distance = arrayfun(@(lats, lons)...
        vdist(lats, lons,...
        hydrophone_struct(5).location(1),...
        hydrophone_struct(5).location(2)),...
        hydrophone_struct(channum).location(1),...
        hydrophone_struct(channum).location(2));
    
    
    % Maximum expected tdoa between channel 1 and 2
    MaxTOA_1 = sqrt(depth_range^2 + horizontal_distance.^2)/1500;
    
    
    % Step through the arrivals and select a call at random from all calls
    % within the MaxTOA_1
    for ii = 1:length(ArrivalsMiss)
        
        corr_idxs = find(abs(Arrivals2Chan(:,channum+1)...
            -Arrivals2Chan(ii,1))<MaxTOA_1+4);
        
        % If more than one index falls in the correlation zone pick first one
        % difference
        if length(corr_idxs)>1
            [~, corr_idxs_win] = min(abs(Arrivals2Chan(corr_idxs,2)...
                -Arrivals2Chan(ii,1)));
            
            corr_idxs = corr_idxs(corr_idxs_win);
            
        elseif isempty(corr_idxs)
            continue;
        end
        ArrivalsMiss(ii, channum+1) = Arrivals2Chan(corr_idxs,...
            channum+1);    
    end
    
end

% Change the 0's to nans to make easier below
ArrivalsMiss(ArrivalsMiss == 0) = nan;

% Create the TDOA values
TDOA_vals = [ArrivalsMiss(:, 1)- ArrivalsMiss(:, 2),...
    ArrivalsMiss(:, 1)-ArrivalsMiss(:, 3)];


%% Create probability space and expanded probability spaces
%ArrivalsMiss = Arrivals2Chan;
% Swim speed (s) and speed of sound (c)
s =12;
c=1500;

% Create structure for each call and it's projected space
Proj_space = struct();


for ii =1:length(ArrivalsMiss)
    
    
    % Pick a set of call delays (number of calls)
    delays = (TDOA_vals(ii,:));
    
    % Remove nans, but keep the indexes for hydrophone pairs
    hyd_idx = find(~isnan(delays));
    delays = delays(~isnan(delays));
    
    averageLSQ_space = [];
    
    % Iterate through the hydrophone pairs and get the combined TDOA
    for jj = 1:length(delays)
        
        % Get the tdoa space
        toa_space = (cell2mat(array_struct(1).toa_diff(jj+1)));
        
        % potential toa space
        LSQ_VAL = (toa_space - delays(jj)).^2;
        
        % Convert tdoa space to probability space
        averageLSQ_space(:,:,jj) = LSQ_VAL;
        
    end
    
    % If there were multiple arrivals of the same call on each hydrophone,
    % take the average of the normalized PDF space
    if length(delays)>1
        
        % average along the third axis
        averageLSQ_space = mean(averageLSQ_space,3);
        
    end
    
    % Do the time projection up to the maximum correlation time
    % figure out which calls are within the alotted fimeframe
    time_gaps = ArrivalsMiss(ii:end, 1)-ArrivalsMiss(ii, 1);
    
    % Only interested in time gaps less than the time cut off (5min)
    time_gaps = time_gaps(time_gaps<time_cut);
    
    % sigma (error) valuse from the normal distribution
    sigma = (time_gaps * (2*s)/c);
    
    % Iterate through  the time gaps and create the projections
    for jj= 1:length(time_gaps)
        
        % Grow the LSQ space
        LSQ_space =  normpdf(averageLSQ_space, 0, (1+sigma(jj)));
        
        Proj_space(ii).projection(:,:,jj) = LSQ_space;
    end
    
end

%% Plots Demsonstrating Time projection of the Probility space

figure
subplot(1,3,1)
imagesc(array_struct.latgrid, array_struct.longrid,...
    (squeeze(Proj_space(1).projection(:,:,1))'))
hold on
scatter(hyd_arr(:,1), hyd_arr(:,2), 80, 'k', 'filled', 'd')
title_str= ['Call 1'];
xlabel('Lat')
ylabel('Lon')
title(title_str)

subplot(1,3,2)
imagesc(array_struct.latgrid, array_struct.longrid,...
    (squeeze(Proj_space(1).projection(:,:,3))'))
hold on
scatter(hyd_arr(:,1), hyd_arr(:,2), 80, 'k', 'filled', 'd')
title_str= ['Projected ',...
    num2str(round(Arrivals2Chan(3,1)-Arrivals2Chan(1,1))), ' sec'];
xlabel('Lat')
ylabel('Lon')
title(title_str)

subplot(1,3,3)
imagesc(array_struct.latgrid, array_struct.longrid,...
    (squeeze(Proj_space(1).projection(:,:,end))'))
hold on
scatter(hyd_arr(:,1), hyd_arr(:,2), 80, 'k', 'filled', 'd')
title_str= ['Projected ',...
    num2str(...
    round((Arrivals2Chan(size(Proj_space(1).projection,3),1)-Arrivals2Chan(1,1))))...
    , ' sec'];
xlabel('Lat')
ylabel('Lon')
title(title_str)
colorbar


%% Create the similarity matrix

%Previously, in coding, we created the TDOA projections for each call and
% projected it up to 10 minutes, stopping where each potentially
% confounding call could be. Now we step through each call space and see if
% it and it's projections match up with the next calls in the sequence

Sim_mat = zeros(length(TDOA_vals))/0;
for ii =1:(length(Proj_space)-1)
    
    % The call we are currently looking at
    Current_call = Proj_space(ii).projection;
    
    % Compare the call and it's temporal projections to the next calls in
    % the sequence. If only one call in the sequence move to the next.
    coef_vals = 1;
    if length(size(Current_call))>2
        
        % Figure out the max number of calls to look ahead, runs into
        % problems near the end otherwise
        
        max_lookahead = min([size(Current_call,3)-1,...
                            (length(Proj_space)-1)-size(Current_call,3)-1+ii]);
        
        % Don't go beyoned the end of the matrix
        if max_lookahead >0
            for jj =1: max_lookahead
                
                % Call a projected at time of call b
                call_a =Current_call(:,:,jj);
                
                % Next call in the sequence, if more than one dimension only
                % take first grid space
                call_b = Proj_space(ii+jj-1).projection;
                
                if length(size(call_b))>2
                    call_b = call_b(:,:,1);
                end
                
                
                % Normalize both call a and call b such that we are looking at
                % the probabiliyt that eacg grid space COULD have the whale in
                % it
                call_b = (call_b - min(min(call_b)))./...
                    ( max(max(call_b)) - min(min(call_b)));
                
                call_a = (call_a - min(min(call_a)))./...
                    ( max(max(call_a)) - min(min(call_a)));
                
                
                % Size of area represented by the next call
                b_hiprob = find(call_b>.99);
                a_hiprob = find(call_a>.99);
                size_b_highprob = length(b_hiprob);
                
                % Determine the percentage of the next call that is covered by
                % the time expanded current call (mm)
                coef = length(intersect(a_hiprob, b_hiprob))/size_b_highprob;
                coef_vals = [coef_vals, coef];
                
            end
        end
    end
    
    % Insert the coefficients
    Sim_mat(ii,ii:(ii+length(coef_vals)-1)) = coef_vals;
    Sim_mat(ii:(ii+length(coef_vals)-1),ii) = coef_vals;
    
    
end

% Plot of the similarity matrix 
figure; 
mycolormap = [ ones(1,3); jet(30)];
imagesc(Sim_mat)
colormap(mycolormap);
colorbar
title('Call Space Simmilarity')
xlabel('Call ID')
ylabel('Call ID')

%% Ladder linkages for new reduced plot
cutoff = 0.9999; % Correlation threshold cut off
time = ArrivalsMiss;

% Elapsed time since first call
time=(time-time(1));  % seconds 
time0=time; % keep original copy for later indexing

% Replicate Corr_coef_map for shrinking matrix
Corr_coef_map = Sim_mat;
[s1,s2]=size(Corr_coef_map);
nc=1; % cluster counter
clear chain dex ss

while s1>1
    flag=1;
    
    % call element index, start at 1 and examine columns
    elt=1; % first element of cluster is always upper left hand element
    index=[];
    cluster=[];
    
    
    % Find all calls that are sufficiently correlated with this one.
    while flag == 1
        % Add call to current cluster, retaining the time and index
        % of the call element
        cluster = [cluster, time(elt)]; 
        index   = [index,elt];
        
        % entire column for each element
        vec=Corr_coef_map(:,elt);
        
        % Only look at the elements below
        vec=vec(elt+1:end);
        
        % Look for correlations that exceed the threshold
        k=find(vec>cutoff);
        
        if length(k)==0
            flag=0;  % Nobody met the criterion
        else
            elt=k(1)+elt;  % Set next call to the first that meets criterion
        end
    end % flag test
    
    % now cut cluster at first large jump:
    k=find(diff(cluster)>time_cut);
    if length(k)>0
        k=k(1);
        % Remove everything after the first large gap
        cluster=cluster(1:k);
        index=index(1:k);
    end
    
    % aggregate clusters
    chain(nc).links=cluster;
    chain(nc).n = length(cluster);
    for j=1:length(cluster)
        chain(nc).index(j)=find(time0==cluster(j));
    end
    nc=nc+1;
    
    % filter out this cluster to reduce matrix
    [s1,s2]=size(Corr_coef_map);
    rs1=setdiff([1:s1],index);
    Corr_coef_map=Corr_coef_map(rs1,:);
    Corr_coef_map=Corr_coef_map(:,rs1);
    time=time(rs1);
    [s1,s2]=size(Corr_coef_map);
    
end % while s1 loop

chain = struct(chain);


%

%% Transform spaceWhales into localization structure

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

for jj = 1:2
    for ii = 1:2
        
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
% 
% 
% % Stick the real id's next to the predicted clusters
% truthAndpreds =[ localize_struct_sim.hyd(parent).agent_ids' newClusterId];
% truthAndpreds(:,3)=0; % True class
% truthAndpreds(:,4)=0; % Score
% truthAndpreds(:,5)=0; % Prediction after applying likelihood ratio

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






