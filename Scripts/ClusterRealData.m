%% Experiment-Clock drift, 
% Determin the performance of the clustering algorithims under various
% clock drift conditions

% Run experiments
close all; clear all
clear classes; clc

whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2});
load(whereAmI{3});
load(whereAmI{4});
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


if ~isfield(localize_struct.hyd(5), 'calls')
    localize_struct.hyd(5).calls = hyd(5).detection.calls;
end
parent =5;
fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;
%% Clear some of the chaff from the localize struct
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

%%

close all; 
    % Populate data and parameters
    examp = clusterClassGPL();
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.cutoff = .9;
    examp.time_cut = 70*60;
    examp.randomMiss =0;
    examp.s =8;
    examp.child_idx = [1:8];
    examp.localize_struct =localize_struct;
    examp.limitTime =36*60*60;
    examp.maxEltTime =60*10;
    examp.calls =localize_struct.hyd(5).calls;
    examp.callParm =  hyd(5).detection.parm;       
    
    examp.clearCalcValues
    UpdateArrTable(examp) % Run this first!
    simMatTDOAonly(examp)
    examp.updateClusterID
    
    examp.exportRavenTxt()
    examp.drawSimMat
    examp.drawAgents
    
    %%
    examp.clearCalcValues
    UpdateArrTable(examp) % Run this first!
    toaOnlyCluster(examp);
    examp.drawAgents
    %%
    examp.clearCalcValues;
    examp.simMatIdeal;
    examp.updateChains;
    examp.updateClusterID;
    examp.drawAgents
%% Compare GPL clustered detections with verified calls


% Load the calls
 load('/home/kpalmer/AnacondaProjects/bin/SpatialContext/Scripts/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_calls1_671.mat')


























