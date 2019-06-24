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


nRuns = 10;
nAgents = round(linspace(3,11,5));

% Thresholds fo rthe snesitivty analysis
% Run a default example first
TimeThresh =linspace(0,20*60,50);
SimThresh = linspace(0,1,50);


perf_methbaseline = struct('RandMat', cell(1, length(nAgents)), 'predAgents', cell(1, length(nAgents)));
perf_meth1 = perf_methbaseline;
perf_meth2 = perf_methbaseline;
perf_meth3 = perf_methbaseline;



%% Run the loop


% Outter Loop - numberof agents in the simulation
for ii =1:length(nAgents)
    disp([num2str(nAgents(ii)), ' Agents'])
    agentNum = nAgents(ii);
    
    for iter=1:nRuns
        
        disp(['Run ', num2str(iter)])
        
        %close all
        clear spaceWhale examp
        %disp([num2str(iter) ' of ' num2str(nRuns) ' runs'])
        % Create new agents
        [spaceWhale] =  createRandomSpaceWhale(0.75, agentNum, hyd_arr,...
            array_struct,hydrophone_struct, ssp, grid_depth, [5, 1,2]); 
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
        [senMat nAgePreds] = runSensitivtyLp(examp,TimeThresh);
        perf_methbaseline(ii).RandMat = cat(3, perf_methbaseline(ii).RandMat, senMat);
        perf_methbaseline(ii).predAgents = cat(3, perf_methbaseline(ii).predAgents, nAgePreds);
        %figure; scatter(nAgents, senMat);
        title('baseline')
        
        
        %% Second method, TDOA only        
        
        % First method, TDOA only
        examp.clearCalcValues();
        examp.time_cut = max(TimeThresh);
        simMatTDOAonly(examp);
        [senMat, nAgePreds] = runSensitivtyLp(examp,TimeThresh,SimThresh);
        
        perf_meth1(ii).RandMat = cat(3, perf_meth1(ii).RandMat, senMat);
        perf_meth1(ii).predAgents = cat(3, perf_meth1(ii).predAgents, nAgePreds);
        %         figure; imagesc(senMat); colorbar
        %         figure; imagesc(nAgents); colorbar
        %figure; scatter(reshape(nAgents,[],1), reshape(senMat,[],1))
        %title('Method 1')
        %% Thrid method, ideal localization
        examp.clearCalcValues();
        examp.time_cut = max(TimeThresh);
        simMatIdeal(examp);
        [senMat, nAgePreds] = runSensitivtyLp(examp,TimeThresh,SimThresh);
        %figure; imagesc(senMat); colorbar
        perf_meth2(ii).RandMat = cat(3, perf_meth2(ii).RandMat, senMat);
        perf_meth2(ii).predAgents = cat(3, perf_meth2(ii).predAgents, nAgePreds);
       % figure; scatter(reshape(nAgents,[],1), reshape(senMat,[],1))
        %title('Method 2')
        %% Fourth method, ad hoc
        examp.clearCalcValues();
        examp.time_cut = max(TimeThresh);
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
quantLvl1 =0;
quantLvl2 =.1;

for jj =1:length(perf_meth1)
    for ii =1:184
        
        aa = perf_meth1(jj).RandMat(:,:,ii);
        bb = perf_meth1(jj).predAgents(:,:,ii);
        agentQuant = quantile(reshape(bb,[],1), [quantLvl1, quantLvl2]);
        valsl =aa(find(bb>=agentQuant(1) & bb<= agentQuant(2)));
        vals = [vals; valsl];
        methLabel = [methLabel; ones(size(valsl))];
        agentsSim = [agentsSim; ones(size(valsl))*nAgents(jj)];
        
        aa = perf_meth2(jj).RandMat(:,:,ii);
        bb = perf_meth2(jj).predAgents(:,:,ii);
        agentQuant = quantile(reshape(bb,[],1), [quantLvl1, quantLvl2]);
       valsl =aa(find(bb>=agentQuant(1) & bb<= agentQuant(2)));
        vals = [vals; valsl];
        methLabel = [methLabel; ones(size(valsl))+1];
        agentsSim = [agentsSim; ones(size(valsl))*nAgents(jj)];

        
        aa = perf_meth3(jj).RandMat(:,:,ii);
        bb = perf_meth3(jj).predAgents(:,:,ii);
        agentQuant = quantile(reshape(bb,[],1), [quantLvl1, quantLvl2]);
        valsl =aa(find(bb>=agentQuant(1) & bb<= agentQuant(2)));
        vals = [vals; valsl];
        methLabel = [methLabel; ones(size(valsl))+2];
        agentsSim = [agentsSim; ones(size(valsl))*nAgents(jj)];

        
        aa = perf_methbaseline(jj).RandMat(:,:,ii);
        bb = perf_methbaseline(jj).predAgents(:,:,ii);
        agentQuant = quantile(reshape(bb,[],1), [quantLvl1, quantLvl2]);
       valsl =aa(find(bb>=agentQuant(1) & bb<= agentQuant(2)));
        vals = [vals; valsl];
        methLabel = [methLabel; ones(size(valsl))+3];
        agentsSim = [agentsSim; ones(size(valsl))*nAgents(jj)];
        
    end
end

figure; 
for ii =1:length(perf_meth1)
    idx = find(agentsSim == nAgents(ii));
    subplot(2,3,ii);boxplot(vals(idx), methLabel(idx))
    title([num2str(nAgents(ii)), ' agents in simulation' ])
    xlabel('Method Number')
    ylabel('Adjusted Rand Index')
    ylim([-.15,1.05])
end





%%

% 
% figure; scatter(reshape(bb,[],1), reshape(aa,[],1))
% randQuant =quantile(reshape(aa,[],1), .9);
% 
% highQuantAgents = bb(find(aa>randQuant));
% hist(highQuantAgents,4)
% [randIDx, randIDy] = find(aa>randQuant);
% hist(aa(sub2ind(size(aa), randIDx, randIDy)));
% hist(bb(sub2ind(size(bb), randIDx, randIDy)));
% 
% figure; scatter(aa(sub2ind(size(aa), randIDx, randIDy)), bb(sub2ind(size(bb), randIDx, randIDy)));
% 
% 
% n_agents = reshape(perf_meth2(1).predAgents(randIDx,randIDx,1);
% 
% timeVal = TimeThresh(randIDy)
% simVal = SimThresh(randIDy)
% scatter(timeVal,simVal)

%% Make Plots of The Sensitivity Space
close all;
for ii =1:length(perf_meth1)
figure

subplot(2,2,1)
plot(TimeThresh/60,(squeeze(nanmedian(perf_methbaseline(ii).RandMat,3))))
xlabel('Time Threshold (min)')
ylabel('Adjusted Ran Index')
title('Baseline')
ylim([-.05 .1])


subplot(2,2,2)
imagesc(SimThresh, TimeThresh/60,squeeze(nanmedian(perf_meth1(ii).RandMat,3)))
%colorbar; caxis([0 .4])
ylabel('Time Threshold (min)')
xlabel('Similarity Threshold')
title('Method 1- TDOA')

subplot(2,2,3)
imagesc(SimThresh, TimeThresh/60,squeeze(nanmedian(perf_meth2(ii).RandMat,3)))
%colorbar; caxis([0 .4])
ylabel('Time Threshold (min)')
xlabel('Similarity Threshold')
title('Method 2- Spatial Ideal')


subplot(2,2,4)
imagesc(SimThresh, TimeThresh/60,squeeze(nanmedian(perf_meth3(ii).RandMat,3)))
%colorbar; caxis([0 .4])
ylabel('Time Threshold (min)')
xlabel('Similarity Threshold')
title('Method 3- Spatial ad hoc')

mtit([num2str(nAgents(ii)), ' agents in Simulation'])

end
































