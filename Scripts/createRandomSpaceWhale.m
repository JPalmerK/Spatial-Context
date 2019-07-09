function spaceWhale = createRandomSpaceWhale(n_hrs, nAgents, hyd_arr,...
    array_struct,hydrophone_struct, ssp, grid_depth, subArray)

clear spaceWhale

% Simulation parameters that apply to the whole system
param_sim.dur= 60*60*n_hrs; % duration of the entire simulation
param_sim.fs =2000; % not currently used but will  be
param_sim.ssp = ssp(1,2); % Sound speed
param_sim.depth = grid_depth; % Agent calling depth

% add the simulation parameters to the spaceWhale
spaceWhale.param_sim = param_sim;


% Assume max swimming speed of 10mph
% http://www.listenforwhales.org/page.aspx?pid=451
parm_movement.max_speed = 8.4;

% currently needs to be the same for all
parm_movement.duration = 60*60*n_hrs; 

% Make call parameters
parm_calling.model = 'Bout';
%parm_calling.model = 'random';


% Create 10 agents with different movement behaviours
for ii=1:nAgents
    % calling every 200 seconds on average
    %parm_calling.frequency = 400*rand; 
    
    % Movement parameters for each  agent(s)
    parm_movement.lat_init =  min(hyd_arr(subArray,1)) + ...
        (max(hyd_arr(subArray,1)) - min(hyd_arr(subArray,1)))*rand(1);
    parm_movement.lon_init =  min(hyd_arr(subArray,2)) +...
        (max(hyd_arr(subArray,2)) - min(hyd_arr(subArray,2)))*rand(1);
    parm_movement.bearing = 180 -360*rand(1); %degrees
    
    % Duration - random number of hours (uniform) times some portion of
    % that (rand)
    
%     parm_movement.duration = floor(...
%         60*60*(randi(n_hrs)*rand(1))); % how long is the agent around
    
% how long is the agent around
    parm_movement.duration = param_sim.dur-2; 
    
    
    parm_movement.start_time = randi([1,...
        param_sim.dur-parm_movement.duration],1);
    
    if mod(ii,2) ==1
        parm_movement.model = 'search';
        
    else
        parm_movement.model = 'directed travel';
    end
    
    % Second in agent's life when it starts calling
    parm_calling.start_time = 1;
    
    % Add the parameter to each agent
    spaceWhale.agent(ii).parm_movement = parm_movement;
    spaceWhale.agent(ii).parm_calling = parm_calling;
    
end



% add the movement
[spaceWhale] = AgentMovement(spaceWhale, array_struct, hydrophone_struct);

end






