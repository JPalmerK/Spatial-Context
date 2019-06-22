% Experiment 1 - Sensitivity Analysis

% Run experiments
close all; clear all
clear classes; clc

whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2})
load(whereAmI{3})


%%
% Parent hydrophone for the analysis
parent =5;
fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;

% Convert the meta data to a structure (same format as GPL expects for
% conistancy with application to real data)

hydrophone_struct= struct();
for ii=1:size(dclde_2013_meta,1)
    hydrophone_struct(ii).name = num2str(dclde_2013_meta(ii,1));
    hydrophone_struct(ii).location = dclde_2013_meta(ii,[11:12]);
    hydrophone_struct(ii).depth= abs(dclde_2013_meta(ii, 13));
    hydrophone_struct(ii).channel=ii;
end
hyd_arr = struct2cell(hydrophone_struct);
hyd_arr =vertcat(hyd_arr{2,:,:});


%% Run the Experiment
nIters = 100;
n_agents = 7;

% Run a default example first
TimeThresh =linspace(0,2000,50);
SimThresh = linspace(0,1,50);

ExpScoresMeth1 = zeros(nIters, length(TimeThresh), length(SimThresh))/0;
ExpScoresMeth2 = ExpScoresMeth1;
ExpScoresMeth3 = ExpScoresMeth1;
ExpScoresBaseline = ExpScoresMeth1;


for iter =1:nIters
    clear spaceWhale examp
    close all
    disp(num2str(iter))
    
    % Create new agents
    [spaceWhale] =  createRandomSpaceWhale(.75, n_agents, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth);
    
    % Populate data and parameters
    examp = simulationClass();
    examp.spaceWhale = spaceWhale;
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.spaceWhale= spaceWhale;
    examp.time_cut = 1500;
    examp.randomMiss =0;
    

    % First method, TDOA only
    examp.clearCalcValues();
    simMatTDOAonly(examp);
   
    for jj = 1:length(SimThresh)
        examp.Cluster_id =[];
        examp.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=(TimeThresh(ii));
            examp.updateChains;
            examp.getRand();
            ExpScoresMeth1(iter, ii,jj) = examp.AdjRand;
            
        end
        
    end
    
    figure
    imagesc(squeeze(ExpScoresMeth1(iter,:,:))); colorbar
    title('Method 1')
    
    
    % Second method, Ideal
    examp.clearCalcValues();
    examp.time_cut = 1500;
    simMatIdeal(examp);
    
    

        
    for ii = 1:length(TimeThresh)
        
        examp.Cluster_id =[];
        examp.time_cut=(TimeThresh(ii));
        examp.updateChains;
        examp.getRand();
        ExpScoresMeth2(iter, 1,ii) = examp.AdjRand;
        
    end
        
    
    figure
    imagesc(squeeze(ExpScoresMeth2(iter,:,:))); colorbar;
    colorbar
    title('Method 2')
    
    % Third method, ad hoc
    examp.clearCalcValues();
    examp.time_cut = 1500;
    simMatadHoc(examp);
    
    
    for jj = 1:length(SimThresh)
        examp.Cluster_id =[];
        examp.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=(TimeThresh(ii));
            examp.updateChains;
            examp.getRand();
            ExpScoresMeth3(iter, ii,jj) = examp.AdjRand;
            
        end
        
    end
    
    figure
    imagesc(squeeze(ExpScoresMeth3(iter,:,:))); colorbar;
    title('Method 3')
    
    
    % Fourth method, baseline
    examp.clearCalcValues();
    examp.toaOnlyCluster();
    examp.getRand();
    
    for jj = 1:length(SimThresh)
        examp.Cluster_id =[];
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=(TimeThresh(ii));
            examp.toaOnlyCluster();
            examp.getRand();
            ExpScoresBaseline(iter, ii,1) = examp.AdjRand;
            
        end
        
    end
    
    figure
    imagesc(squeeze(ExpScoresBaseline(iter,:,:)))
    title('Baseline')
    colorbar
    
    
end
%%

% Mean Values
figure
subplot(2,2,1)
plot(TimeThresh/60, nanmedian(ExpScoresBaseline(:,:,1),1))
xticks((1:9:length(TimeThresh)))
xticklabels({round(TimeThresh(1:9:end),2)})
ylabel('Median Adjusted Rand')
title('Baseline')

subplot(2,2,2)
imagesc(squeeze(nanmedian(ExpScoresMeth1)))
colorbar;
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 1')

subplot(2,2,3)
imagesc(squeeze(nanmedian(ExpScoresMeth2)))
colorbar; 
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 2')

subplot(2,2,4)
imagesc(squeeze(nanmedian(ExpScoresMeth3)))
colorbar; 
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 3')


%% Figure of Upper 95th percentile


figure
subplot(2,2,1)
plot(TimeThresh/60, squeeze(prctile(ExpScoresBaseline(:,:,1),95)))
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
ylabel('95th Percentile Adjusted Rand')
title('Baseline')

subplot(2,2,2)
imagesc(squeeze(prctile(ExpScoresMeth1,95)))
colorbar;
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 1')

subplot(2,2,3)
imagesc(squeeze(prctile(ExpScoresMeth2,95)))
colorbar; 
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 2')

subplot(2,2,4)
imagesc(squeeze(prctile(ExpScoresMeth3,95)))
colorbar; 
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 3')

%% Figure of lower 95th percentile


figure
subplot(2,2,1)
plot(TimeThresh/60, squeeze(prctile(ExpScoresBaseline(:,:,1),5)))
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
ylabel('95th Percentile Adjusted Rand')
title('Baseline')

subplot(2,2,2)
imagesc(squeeze(prctile(ExpScoresMeth1,5)))
colorbar;
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 1')

subplot(2,2,3)
imagesc(squeeze(prctile(ExpScoresMeth2,5)))
colorbar; 
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 2')

subplot(2,2,4)
imagesc(squeeze(prctile(ExpScoresMeth3,5)))
colorbar; 
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 3')