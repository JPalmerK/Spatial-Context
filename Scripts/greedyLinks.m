%% Experiment-Clock drift,
% Determin the performance of the clustering algorithims under various
% clock drift conditions

% Run experiments
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

% Number of hours the experiment runs and the number of agents included
n_hrs = 0.75;
n_agents = 6;


clear spaceWhale examp
close all
%
spaceWhale=[];
% Create new agents
[spaceWhale] =  createRandomSpaceWhale(n_hrs,n_agents, hyd_arr,...
    array_struct,hydrophone_struct, ssp, grid_depth, array_struct.slave([1,2,3]));

% Populate data and parameters
examp = simulationClass();
examp.spaceWhale = spaceWhale;
examp.array_struct = array_struct;
examp.hydrophone_struct = hydrophone_struct;
examp.spaceWhale= spaceWhale;
examp.cutoff = .95;
examp.time_cut = 10*60;
examp.randomMiss =0;
examp.child_idx = [1,2,3];


% Second method, ideal
examp.clearCalcValues();
simMatIdeal(examp);
examp.drawSimMat

simmat = examp.Sim_mat;
exampChains = examp.chains;

%% Create the chains going backwards




times0 = examp.arrivalArray(:,1);

% Step through looking for chains
simmatTemp = simmat;
cutoff = .85;
% val =[];
% time =[];
% 
% for ii=1:size(simmatTemp,1)-1
%     val = [val; simmatTemp(ii+1:end,ii)];
%     time = [time; times0(ii+1:end)];
%     
% end
% scatter(time, val)


chains =struct;
ids =1:length(times0);
times = times0;
clusterId =1;


chainID =1;
curr_id =1;


while length(times0)>2
links = times0(curr_id)
index = ids(curr_id)
vals = simmatTemp(curr_id, curr_id+1:end)';
timevals = times0(curr_id+1:end)-times0(curr_id); % elapsed time since call
idcalls = ids(curr_id+1:end);

while nanmax(vals)>cutoff
    curr_id
    normvals = vals./timevals;
    [~,id] = max(normvals);
    links = [links timevals(id)];
    index = [index idcalls(id)];
    curr_id = idcalls(id);
    
    vals = simmatTemp(curr_id, curr_id+1:end)';
    timevals = times0(curr_id+1:end)-times0(curr_id); % elapsed time since call
    idcalls = ids(curr_id+1:end);
end

    % end of cluster remove chosen values

    chains(chainID).links =links;
    chains(chainID).index = index;
    chains(chainID).n = length(index);
        times0(index)=[];
    simmatTemp(:,index)=[];
    ids(index)=[];
    chainID =chainID+1
    curr_id =1;
    index=[];
    links =[];
end  
    
    






















