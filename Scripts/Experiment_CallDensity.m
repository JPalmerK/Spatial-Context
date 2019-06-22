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


%% Set up experiment thresholds and results structre

% Similarity thresholds from the sensitivity analysis
% Method 1 - TDOA only
% Method 2 - Spatial Grids
% Method 3 - Spatial GRids, ad hoc
% Method 4/Baseline - TOA only


simThresh = [0.91, 0.45, 0.55];
timeThresh = [41 13*60 14*60 7*60];
nRuns = 100;
nAgents = round(linspace(3,10,6));

% Thresholds fo rthe snesitivty analysis
% Run a default example first
TimeThresh =linspace(0,2000,50);
SimThresh = linspace(0,1,50);


perf_methbaseline = struct('RandMat', cell(1, length(nAgents)));
perf_meth1 = perf_methbaseline;
perf_meth2 = perf_methbaseline;
perf_meth3 = perf_methbaseline;



%% Run the loop


% Outter Loop - numberof agents in the simulation
for ii =1:length(nAgents)
    disp([num2str(nAgents), ' Agents'])
    agentNum = nAgents(ii);
    
    for iter=1:nRuns
        
        disp(['Run ', num2str(iter)])
        
        close all
        clear spaceWhale examp
        %disp([num2str(iter) ' of ' num2str(nRuns) ' runs'])
        % Create new agents
        [spaceWhale] =  createRandomSpaceWhale(0.75,agentNum, hyd_arr,...
            array_struct,hydrophone_struct, ssp, grid_depth);
        
        % Populate data and parameters
        examp = simulationClass();
        examp.spaceWhale = spaceWhale;
        examp.array_struct = array_struct;
        examp.hydrophone_struct = hydrophone_struct;
        examp.spaceWhale= spaceWhale;
        examp.time_cut = 10*60;
        examp.randomMiss =0;
        examp.UpdateArrArray()
        
        %% First method, baseline
        examp.clearCalcValues();
        senMat = runSensitivtyLp(examp,TimeThresh);
        perf_methbaseline(ii).RandMat = cat(3, perf_methbaseline(ii).RandMat, senMat);
        
        
        
        %% Second method, TDOA only        

        % First method, TDOA only
        examp.clearCalcValues();
        examp.time_cut = max(TimeThresh);
        simMatTDOAonly(examp);
        senMat = runSensitivtyLp(examp,TimeThresh,SimThresh);
        figure; imagesc(senMat); colorbar
        perf_meth1(ii).RandMat = cat(3, perf_meth1(ii).RandMat, senMat);
        
        
        
        %% Thrid method, ideal localization
        examp.clearCalcValues();
        examp.time_cut = max(TimeThresh);
        simMatIdeal(examp);
        senMat = runSensitivtyLp(examp,TimeThresh,SimThresh);
        figure; imagesc(senMat); colorbar
        perf_meth2(ii).RandMat = cat(3, perf_meth2(ii).RandMat, senMat);
        
        %% Fourth method, ad hoc
        examp.clearCalcValues();
        examp.time_cut = max(TimeThresh);
        simMatadHoc(examp);
        senMat = runSensitivtyLp(examp,TimeThresh,SimThresh);
        figure; imagesc(senMat); colorbar
        perf_meth3(ii).RandMat = cat(3, perf_meth3(ii).RandMat, senMat);
        

    end
    
    
    
    
end

%% Make Plots


figure
subplot(2,2,2)
imagesc(squeeze(nanmedian(perf_meth1(1).RandMat,3)))
colorbar; caxis([0 .4])

subplot(2,2,3)
imagesc(squeeze(nanmedian(perf_meth2(1).RandMat,3)))
colorbar; caxis([0 .4])

subplot(2,2,4)
imagesc(squeeze(nanmedian(perf_meth3(1).RandMat,3)))
colorbar; caxis([0 .4])





