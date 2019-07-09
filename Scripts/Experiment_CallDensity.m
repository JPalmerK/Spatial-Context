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


nRuns = 100;
nAgents = round(linspace(3,9,4));
% Thresholds fo rthe snesitivty analysis
% Run a default example first
TimeThresh =linspace(0,120,50);
SimThresh = linspace(0,1,50);


perf_methbaseline = struct('RandMat', cell(1, length(nAgents)), 'predAgents', cell(1, length(nAgents)));
perf_meth1 = perf_methbaseline;
perf_meth2 = perf_methbaseline;
perf_meth3 = perf_methbaseline;



%% Run the loop

child_idx = [1 2 3];

% Outter Loop - numberof agents in the simulation
for ii =1:length(nAgents)
    disp([num2str(nAgents(ii)), ' Agents'])
    agentNum = nAgents(ii);
    
    for iter=1:nRuns
        
        disp(['Run ', num2str(iter)])
        
        %close all
        clear spaceWhale examp
        %disp([num2str(iter) ' of ' num2str(nRuns) ' runs'])
        % Create new agents- id values of the hydrophone array
        [spaceWhale] =   createRandomSpaceWhale(0.75,8, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth,...
        [array_struct.master, array_struct.slave(child_idx)]);
    
        % Populate data and parameters
        examp = simulationClass();
        examp.spaceWhale = spaceWhale;
        examp.array_struct = array_struct;
        examp.hydrophone_struct = hydrophone_struct;
        examp.spaceWhale= spaceWhale;
        examp.randomMiss =0;
        examp.UpdateArrArray();
        examp.child_idx = child_idx;
        examp.time_cut = 5*60;
        UpdateArrArray(examp)
        
        
        %% First method, baseline
        examp.clearCalcValues();
        toaOnlyCluster(examp);
        [senMat nAgePreds] = runSensitivtyLp(examp,TimeThresh);
        perf_methbaseline(ii).RandMat = cat(3, perf_methbaseline(ii).RandMat, senMat);
        perf_methbaseline(ii).predAgents = cat(3, perf_methbaseline(ii).predAgents, nAgePreds);
        %figure; scatter(nAgents, senMat);
        %title('baseline')
        
        
        %% Second method, TDOA only        
        
        % First method, TDOA only
        examp.clearCalcValues();
        simMatTDOAonly(examp);
        [senMat, nAgePreds] = runSensitivtyLp(examp,TimeThresh,SimThresh);
        perf_meth1(ii).RandMat = cat(3, perf_meth1(ii).RandMat, senMat);
        perf_meth1(ii).predAgents = cat(3, perf_meth1(ii).predAgents, nAgePreds);

        %% Thrid method, ideal localization
        examp.clearCalcValues();
        simMatIdeal(examp);
        [senMat, nAgePreds] = runSensitivtyLp(examp,TimeThresh,SimThresh);
        %figure; imagesc(senMat); colorbar
        perf_meth2(ii).RandMat = cat(3, perf_meth2(ii).RandMat, senMat);
        perf_meth2(ii).predAgents = cat(3, perf_meth2(ii).predAgents, nAgePreds);
       % figure; scatter(reshape(nAgents,[],1), reshape(senMat,[],1))
        %title('Method 2')
        %% Fourth method, ad hoc
        examp.clearCalcValues();
        simMatadHoc(examp);
        [senMat, nAgePreds] = runSensitivtyLp(examp,TimeThresh,SimThresh);
        %figure; imagesc(senMat); colorbar
        perf_meth3(ii).RandMat = cat(3, perf_meth3(ii).RandMat, senMat);
        perf_meth3(ii).predAgents = cat(3, perf_meth3(ii).predAgents, nAgePreds);
        %figure; scatter(reshape(nAgents,[],1), reshape(senMat,[],1))
        %title('Method 3')

    end
    
    
    
    
end



%% Try to link number of agents and the thresholds

vals = [];
methLabel =[];
agentsSim =[];
agentgroupid=[];
quantLvl1 =2;
quantLvl2 =100;
agentQuant=[2 100]
for jj =1:length(nAgents)
    for ii =1:nRuns
        
        aa = perf_meth1(jj).RandMat(:,:,ii);
        bb = perf_meth1(jj).predAgents(:,:,ii);
        %agentQuant = quantile(reshape(bb,[],1), [quantLvl1, quantLvl2]);
        valsl =aa(find(bb>=agentQuant(1) & bb<= agentQuant(2) & bb>1));
        vals = [vals; valsl];
        methLabel = [methLabel; ones(size(valsl))];
        agentsSim = [agentsSim; ones(size(valsl))*nAgents(jj)];
        
        aa = perf_meth2(jj).RandMat(:,:,ii);
        bb = perf_meth2(jj).predAgents(:,:,ii);
        %agentQuant = quantile(reshape(bb,[],1), [quantLvl1, quantLvl2]);
        valsl =aa(find(bb>=agentQuant(1) & bb<= agentQuant(2) & bb>1));
        vals = [vals; valsl];
        methLabel = [methLabel; ones(size(valsl))+1];
        agentsSim = [agentsSim; ones(size(valsl))*nAgents(jj)];

        
        aa = perf_meth3(jj).RandMat(:,:,ii);
        bb = perf_meth3(jj).predAgents(:,:,ii);
        %agentQuant = quantile(reshape(bb,[],1), [quantLvl1, quantLvl2]);
        valsl =aa(find(bb>=agentQuant(1) & bb<= agentQuant(2) & bb>1));
        vals = [vals; valsl];
        methLabel = [methLabel; ones(size(valsl))+2];
        agentsSim = [agentsSim; ones(size(valsl))*nAgents(jj)];

        
        aa = perf_methbaseline(jj).RandMat(:,:,ii);
        bb = perf_methbaseline(jj).predAgents(:,:,ii);
        %agentQuant = quantile(reshape(bb,[],1), [quantLvl1, quantLvl2]);
        valsl =aa(find( bb>1));
        vals = [vals; valsl];
        methLabel = [methLabel; ones(size(valsl))+3];
        agentsSim = [agentsSim; ones(size(valsl))*nAgents(jj)];
        
    end
end

figure; 
for ii =1:length(perf_meth1)
    idx = find(agentsSim == nAgents(ii));
    subplot(2,3,ii);
    boxplot(vals(idx), methLabel(idx))
    title([num2str(nAgents(ii)), ' agents in simulation' ])
    xlabel('Method Number')
    ylabel('Adjusted Rand Index')
    ylim([-.15,1.05])
end




%% Make Plots of The Sensitivity Space
close all;
tics = fliplr(TimeThresh);

for ii =1:length(perf_meth1)
figure

yval = max([ max(max(squeeze(nanmedian(perf_meth1(ii).RandMat,3)))),...
     max(max(squeeze(nanmedian(perf_meth2(ii).RandMat,3)))),...
      max(max(squeeze(nanmedian(perf_meth3(ii).RandMat,3))))]);

subplot(2,2,1)
plot(TimeThresh/60,(squeeze(nanmedian(perf_methbaseline(ii).RandMat,3))))
xlabel('Time Threshold (s)')
ylabel('Adjusted Ran Index')
title('Baseline')



subplot(2,2,2)
imagesc(SimThresh, tics, flipud(squeeze(nanmedian(perf_meth1(ii).RandMat,3)))')
set(gca,'YDir','normal')
colorbar; caxis([0 yval])
xlabel('Time Threshold (s)')
ylabel('Similarity Threshold')
title('Method 1- TDOA')



subplot(2,2,3)
imagesc(SimThresh, tics, flipud(squeeze(nanmedian(perf_meth2(ii).RandMat,3)))')
set(gca,'YDir','normal')
colorbar; caxis([0 yval])
xlabel('Time Threshold (s)')
ylabel('Adjusted Ran Index')
title('Method 2- Spatial Ideal')


subplot(2,2,4)
imagesc(SimThresh, tics, flipud(squeeze(nanmedian(perf_meth3(ii).RandMat,3)))')
set(gca,'YDir','normal')
colorbar; caxis([0 yval])
xlabel('Time Threshold (s)')
ylabel('Adjusted Ran Index')
title('Method 3- Spatial ad hoc')

mtit([num2str(nAgents(ii)), ' agents in Simulation'])

end
































