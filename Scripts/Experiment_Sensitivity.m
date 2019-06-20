% Experiment 1 - Sensitivity Analysis



% Run experiments
close all; clear all
clear classes; clc

cd('/home/kpalmer/AnacondaProjects/Spatial-Context/Scripts')
% Load metadata to create Hydrophone Structur
dclde_2013_meta = xlsread(strcat('/cache/kpalmer/quick_ssd/data/',...
    'DCLDE_2013_10Channel/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx'));
% Load the DCLDE array structure (just for the array structure, depth, and SSP)
load('DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_array_struct1_671.mat')
load('DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_localize_struct_671.mat')

%%
close all
clear all
cd('D:\Anaconda Projects\Spatial-Context\Scripts')
% Load metadata to create Hydrophone Structur
dclde_2013_meta = xlsread('C:\Users\Kaitlin\Desktop\DCL2013_NEFSC_SBNMS_200903_metadata.xlsx');
load('D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_array_struct1_671.mat')
load('D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_localize_struct_671.mat')

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
nIters = 5;
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
    [spaceWhale] =  createRandomSpaceWhale(.75, 7, hyd_arr,...
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
    
    
    for jj = 1:length(SimThresh)
        examp.Cluster_id =[];
        examp.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=(TimeThresh(ii));
            examp.updateChains;
            examp.getRand();
            ExpScoresMeth2(iter, ii,jj) = examp.AdjRand;
            
        end
        
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
            ExpScoresBaseline(iter, ii,jj) = examp.AdjRand;
            
        end
        
    end
    
    figure
    imagesc(squeeze(ExpScoresBaseline(iter,:,:)))
    title('Baseline')
    colorbar
    
    
end
%%

% get range of values
figure
subplot(2,2,1)
imagesc(nanmean(ExpScoresBaseline,3))
colorbar;
yticks((1:2:length(TimeThresh)))
yticklabels({round(TimeThresh(1:2:end)/60,2)})
title('Baseline')

subplot(2,2,2)
imagesc(squeeze(nanmean(ExpScoresMeth1)))
colorbar;
yticks((1:6:length(TimeThresh)))
yticklabels({round(TimeThresh(1:6:end)/60,2)})
xticks((1:9:length(SimThresh)))
xticklabels({round(SimThresh(1:9:end),2)})
title('Method 1')

subplot(2,2,3)
imagesc(squeeze(nanmean(ExpScoresMeth2)))
colorbar; 
yticks((1:2:length(TimeThresh)))
yticklabels({round(TimeThresh(1:2:end)/60,2)})
xticks((1:4:length(SimThresh)))
xticklabels({round(SimThresh(1:4:end),2)})
title('Method 2')

subplot(2,2,4)
imagesc(squeeze(nanmean(ExpScoresMeth3)))
colorbar; 
yticks((1:2:length(TimeThresh)))
yticklabels({round(TimeThresh(1:2:end)/60,2)})
xticks((1:4:length(SimThresh)))
xticklabels({round(SimThresh(1:4:end),2)})
title('Method 3')


