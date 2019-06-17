function [spaceWhale] = AgentMovement(spaceWhale, array_struct, hydrophone_struct)

% Function for creating animal movement in the enviornment.
% Returnse a structure containing the XY coordinates for each agent
% in space whale as well as the times at which a call is produced

% Input:
% spaceWhale- structure containing agents with associated movmeent
%               and calling parameters
% Movement model options- 'linear', 'directed travel', 'random walk /
% searching'
% array_struct - from GPL
% hydrophone_stuct - structure containing hydrophone locations from
%   gpl setup

% Returns:
%   Adds location, calling times, call arrival times for each
%   hydrophone and each agent (whale)
%   tdoa for each parent/child pair

% This probably wants to be turned into an object...


%% Grid Setup and received time setup
% Meters per decimal degree lat and lon

% Get distance in meters between the lower and upper right
grid_v = vdist(min(array_struct.latgrid), min(array_struct.longrid),...
    max(array_struct.latgrid), min(array_struct.longrid));


% Get distance in meters between the lower left and lower right
grid_h = vdist(min(array_struct.latgrid), max(array_struct.longrid),...
    min(array_struct.latgrid), min(array_struct.longrid));


% Grid X/Y space
deltalat_space = grid_v/ (length(array_struct.latgrid)-1);
deltalon_space = grid_h/ (length(array_struct.longrid)-1);

% Create new lat/lon grids
lat_grid_m = [0:deltalat_space:grid_v];
lon_grid_m = [0:deltalon_space:grid_h];

% Decimal degrees per m lat and lon
deg_per_m_lat = (max(array_struct.latgrid)-...
    min(array_struct.latgrid))/max(lat_grid_m);

deg_per_m_lon = (max(array_struct.longrid)...
    -min(array_struct.longrid))/max(lon_grid_m);


% All simulation times (in seconds)
all_tt = [0:spaceWhale.param_sim.dur];


%% Create movement model with Calling Behaviour

% Get a list of agents where locations (required variable) have not been
% filled out

% First check if location is a field, if not add it
if ~ isfield(spaceWhale.agent, 'location')
    for ii =1:length(spaceWhale.agent)
        spaceWhale.agent(ii).location = [];
    end
    
end


% Determine which agents have not had localizations added. This step allows
% for agents to be added after the fact without updating the existing agent
% locations
if length(spaceWhale.agent) >1
    empty_locs = find(cellfun(@isempty,{spaceWhale.agent.location}));
else
    empty_locs = 1;
end



% For each agent where there are not locations presently existing, make
% the movement model as directed by the calling and movement parameters

for zz= 1:length(empty_locs)
    
    % Go to the agent that has not had it's location created
    ii = empty_locs(zz);
    
    % Extract parameters for the agents moving and calling behaviour
    calling_parms = spaceWhale.agent(ii).parm_calling;
    movement_parms = spaceWhale.agent(ii).parm_movement;
    
    
    % Time in seconds, starts when animal appears ends after duration
    tt = [movement_parms.start_time:...
        movement_parms.start_time+ movement_parms.duration];
    
    % Initialize bearing
    bearing=zeros(1, length(tt));
    
    
    % Create matrix for arrival times
    Arrival_times = [];
    
    % determine linear or random walk model
    switch movement_parms.model
        case 'linear'
            % Travel in a straight line - out of date don't use
            % bearing and speed
%             bearing = (pi/180)*bearing+movement_parms.bearing;
%             speed= zeros(1, length(tt))+ movement_parms.max_speed;
            
            
        case 'directed travel'
            % how often to sample (every n seconds- 45 sec to 5 min)
            n = max(5, randi([floor(spaceWhale.agent(zz).parm_movement.duration/300),...
                floor(spaceWhale.agent(zz).parm_movement.duration/45)]));
            disp([num2str(n) ' movements directed travel'])
            
            % directional travel but vearing
            speed_dist = makedist('Normal', movement_parms.max_speed-2,1);
            speed_dist = truncate(speed_dist, 0,  movement_parms.max_speed);
            
            % useful catch for Alaxander days (Some days are like that,
            % even in Australia)
            try
                sample_time = [1, sort(datasample(tt((2:end-1)), n-2,...
                    'Replace', false)), tt(end)];
            catch ME
                print('blarg')
            end
            
            % Agent speed
            speed = random(speed_dist, length(sample_time), 1)';
            
            % vear plus or minus 20 deg while walking
            bearing = (pi/180)*(movement_parms.bearing +...
                40-80*rand(length(sample_time),1))';
            
            % integrate over all time
            try
                speed = interp1(sample_time, speed, tt, 'nearest');
            catch ME
                disp('blarg')
            end
            
            % Smooth between the bearings
            bearing = interp1(sample_time, bearing, tt, 'nearest');
            
        case {'random walk', 'search'}
            % how often to sample (every n seconds 15 seconds to 1.5 minutes)
            n = max(5,(randi([...
                floor(spaceWhale.agent(zz).parm_movement.duration/120),...
                floor(spaceWhale.agent(zz).parm_movement.duration/30)])));
            
            disp([num2str(n) ' movements random walk'])
            % Random walk movement pattern
            speed_dist = makedist('Normal', movement_parms.max_speed-2,1);
            speed_dist = truncate(speed_dist, 0,  movement_parms.max_speed);
            
            % Times at which movement is observed
            sample_time = [1, sort(datasample(tt(2:end-1), n-2, 'Replace', false)), tt(end)];
            
            % Speed and bearing of the agent
            speed = random(speed_dist, length(sample_time), 1)';
            bearing = (pi/180)*(360*rand(1,length(sample_time)));
            
            % interpolate the location throughout time for smoother
            % movement
            speed = interp1(sample_time, speed, tt, 'nearest');
            bearing = interp1(sample_time, bearing, tt, 'nearest');
            
        case {'brownian'}
            % Brownian movement pattern- not implemented
            %
%             xdis = cumsum(normrnd(length(tt), 0,1 ))
%             ydis = cumsum(normrnd(length(tt), 0,1 ))
            
        case {'fill space'}
            % Fill the grid area with sources - not implemented
%             whale_lat = sort(repmat(array_struct.latgrid, 1,...
%                 length(array_struct.longrid)))';
%             whale_lon =repmat(array_struct.longrid, 1,...
%                 length(array_struct.longrid))';
%             
%             spaceWhale.agent(ii).location = [whale_lat, whale_lon];
%             
%             % Set/override callin parameters
%             calling_parms.model= 'fill space';
            
        case{'stationary'}
%             % One stationary agent in the array - not implemented
%             whale_lat = movement_parms.lat_init ;
%             whale_lon = movement_parms.lon_init ;
%             
%             spaceWhale.agent(ii).location = [whale_lat, whale_lon];
%             calling_parms.model= 'fill space';
%             
            
        otherwise
            disp('Movement model not present');
            
            
    end
    
    
    
    % Check if lat/lon need to be filled in from
    if isfield(spaceWhale.agent(ii), 'location') ==0 ||...
            isempty(spaceWhale.agent(ii).location)
        
        % Convert bearing to lat and lon
        [lon_dist, lat_dist] = pol2cart(bearing, 1);
        
        % Whale latitude and logitude referenced from the pram file (WGS)
        whale_lat=  movement_parms.lat_init +...
            cumsum(speed.*lat_dist*deg_per_m_lat);
        whale_lon=  movement_parms.lon_init +...
            cumsum(speed.*lon_dist*deg_per_m_lon);
        
        % add location information to agent
        spaceWhale.agent(ii).location = [whale_lat', whale_lon'];
    end
    
    
    % Determine when in the track calls are produced
    switch calling_parms.model
        
        % Case one, animal makes calls at a constant rate
        case 'constant'
            % Constant calling rate, e.g. pseudo tracking animal
            
            % Check if calling frequency (i.e. i-call-i is present, break
            % if not)
            
            if ~isfield(calling_parms, 'frequency')
                
                disp(['No call rate (frequency) provided', newline, ...
                    'add call rate to calling_parms.frequency'])
                return
                
            end
            
            % times (indexes) when calls were produced
            call_produced = floor([calling_parms.start_time:...
                calling_parms.frequency:...
                length(tt)]);
        case 'random'
            % Produce calls at random times
            % Number of times to produce a call
            % how often to sample (every n seconds- 30 sec to 10 min)
            n_samples = randi([floor(...
                spaceWhale.agent(zz).parm_movement.duration/600),...
                floor(spaceWhale.agent(zz).parm_movement.duration/30)]);
            
            % Times at which calls were produced
            call_produced =   sort(datasample(1:length(tt),...
                n_samples, 'Replace', false));
            
        case 'Bout'
            % Calls produced in bouts, nested poisson distribution not
            % implemented
            
        case 'fill space'
            % For spatial analysis using non-agent based method. Fill all
            % locations with calls. Not implemented
            call_produced = 1:size(spaceWhale.agent(ii).location,1);
            
        otherwise
            disp('Go away or we shall taunt you a secONd time.')
    end
    
    % Update the agent with calling times
    spaceWhale.agent(ii).call_times = call_produced;
    
    % Update based on agent start time
    
    % For each hydrophone get the range between the calling animal
    % and the hydrophone, then calculate the arrival time
    TL = zeros(length(call_produced),length(hydrophone_struct));
    Range = zeros(length(call_produced),length(hydrophone_struct));
    
    
    for jj=1:length(hydrophone_struct)
        
        % depth between calling whale and array
        depth_range = hydrophone_struct(jj).depth -...
            spaceWhale.param_sim.depth;
        
        
        % Calculate the horizontal distance between the calling
        % locations and the hydrophone
        horizontal_distance = arrayfun(@(lats, lons)...
            vdist(lats, lons,...
            hydrophone_struct(jj).location(1),...
            hydrophone_struct(jj).location(2)),...
            spaceWhale.agent(ii).location(call_produced,1),...
            spaceWhale.agent(ii).location(call_produced,2));
        
        
        %arrival time between call and hydrophone jj
        dist = (depth_range^2 + horizontal_distance.^2).^(0.5);
        
        % Time the call arrives at the hydrophone after it's produced
        Arrival_times(:,jj) = spaceWhale.agent(ii).call_times' + ...
            (dist / spaceWhale.param_sim.ssp);
        
        % Transmission loss cylindirical spreading
        TL(:,jj) = 10*log10(dist);
        Range(:,jj) = dist;
        
    end
    
    
    % Update the agent call arrival times adjust for start time
    spaceWhale.agent(ii).Arrival_times = Arrival_times +...
        movement_parms.start_time;
    spaceWhale.agent(ii).TL = TL;
    spaceWhale.agent(ii).RangeKm = Range/1000;
    
    % Create TDOA space for a combination of every hydrophone
    for jj=1:length(hydrophone_struct)
        master= jj;
        
        for kk=1:10
            % TDOA between hydrophone 1 and remaining hydrophones
            TDOA_struct.hyd(master).td(:,kk)= Arrival_times(:,kk)...
                -Arrival_times(:,master);
            
            
        end
    end
    
    
    % Add tdoa structure to agent
    spaceWhale.agent(ii).tdoa = TDOA_struct.hyd;
    clear TDOA_struct
    
    
end





end