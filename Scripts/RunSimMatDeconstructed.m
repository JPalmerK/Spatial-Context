%% Experiment- Number of Agents
% Setup sensitivity to number of agents

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

child_idx = [1 2 3];
%% Create the agents



[spaceWhale] =   createRandomSpaceWhale(0.75, 9, hyd_arr,...
    array_struct,hydrophone_struct, ssp, grid_depth,...
    [array_struct.master, array_struct.slave(child_idx)]);

simStruct=struct();
simStruct.spaceWhale=spaceWhale;
simStruct.array_struct=array_struct;
simStruct.truncateKm=10;
simStruct.s = 8
simStruct.PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
% Position/soundspeed uncertainty, see Eva's paper
% receiver position = 5 m -> 5/1500 = 0.0004 sec
% height LSQ peak (check with Tyler) = .1 sec
% sound speed profile = 50 m/s -> 50/1500 =0.3 sec
simStruct.drift=0;
simStruct.maxEltTime =60;
simStruct.cutoff = .25;
simStruct.child_idx= child_idx;
simStruct.randomMiss=0;
simStruct.arrivalTable = UpdateArrTable(simStruct);
simStruct.arrivalArray= UpdateArrArray(simStruct);
simStruct.TDOA_vals = UpdateTDOA(simStruct);
simStruct.Sim_mat= simMatIdealXcorrDist(simStruct);
simStruct.chains =updateChainsEncounterFirst(simStruct);
simStruct.Cluster_id= updateClusterID(simStruct)

AdjRand= getRand(simStruct)


