% Accuracy comparison between the three methods using adjusted rand

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
% 
% clear all
% % Load metadata to create Hydrophone Structur
% dclde_2013_meta = xlsread('C:\Users\Kaitlin\Desktop\DCL2013_NEFSC_SBNMS_200903_metadata.xlsx');
% load('D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_array_struct1_671.mat')
% load('D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_localize_struct_671.mat')

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

%% Experiment 1 Ideal method assuming no drift

% Number of hours the experiment runs and the number of agents included
n_hrs = 0.5;
n_agents = 7;
nRuns = 200;

% Number of calls included in the analysis
n_calls = zeros(1, nRuns);

Meth1_aRand = zeros(1, nRuns);
Meth2_aRand = zeros(1, nRuns);
Meth3_aRand = zeros(1, nRuns);
Meth4_aRand = zeros(1, nRuns);



Meth1_aRandDrift_4 = zeros(1, nRuns);
Meth2_aRandDrift_4 = zeros(1, nRuns);
Meth3_aRandDrift_4 = zeros(1, nRuns);
Meth4_aRandDrift_4 = zeros(1, nRuns);


Meth1_aRandDrift_60 = zeros(1, nRuns);
Meth2_aRandDrift_60 = zeros(1, nRuns);
Meth3_aRandDrift_60 = zeros(1, nRuns);
Meth4_aRandDrift_60 = zeros(1, nRuns);


for ii =1:nRuns
    disp(num2str(ii))

    clear spaceWhale examp
    close all
%     
    spaceWhale=[];
    % Create new agents
    [spaceWhale] =  createRandomSpaceWhale(n_hrs,n_agents, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth);
    
    % Populate data and parameters
    examp = simulationClass();
    examp.spaceWhale = spaceWhale;
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.spaceWhale= spaceWhale;
    examp.time_cut = 15*60;
    examp.randomMiss =0;
    
    jj =1;
    
    disp(['Run ', num2str(ii), ' exp ',num2str(jj), ' of 12'])
    
    % First method, TDOA only
    simMatTDOAonly(examp);
    examp.getRand();
    Meth1_aRand(ii)=examp.AdjRand;
    

    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    
    % Second Method, ad hoc
    examp.clearCalcValues();
    simMatadHoc(examp);
    examp.getRand();
    Meth2_aRand(ii)=examp.AdjRand;
    
    
   jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    
    % Third method, ideal
    examp.clearCalcValues();
    simMatIdeal(examp);
    examp.getRand();
    Meth3_aRand(ii)= examp.AdjRand;
    
    
    % Fourth model Baseline
    examp.clearCalcValues();
    toaOnlyCluster(examp);
    examp.getRand();
    Meth4_aRand(ii)= examp.AdjRand;

    %%%%%%%%%%%% Slight Drift Present %%%%%%%%%%%%%%%
    

    % Shift by two seconds allow for 4 seconds of drift (default)
    examp.randomMiss =1;
    examp.assSec = 4;
    examp.drift=2;
    
    % Reset the parameters
    examp.clearCalcValues();    

    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % First method, TDOA only
    examp.clearCalcValues();
    simMatTDOAonly(examp);
    examp.getRand();
    Meth1_aRandDrift_4(ii)=examp.AdjRand;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Second Method, ad hoc
    examp.clearCalcValues();
    simMatadHoc(examp);
    examp.getRand();
    Meth2_aRandDrift_4(ii)=examp.AdjRand;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Third method, ideal
    examp.clearCalcValues();
    simMatIdeal(examp);
    examp.getRand();
    Meth3_aRandDrift_4(ii)= examp.AdjRand;  
    
    % Fourth model Baseline
    examp.clearCalcValues();
    toaOnlyCluster(examp);
    examp.getRand();
    Meth4_aRandDrift_4(ii)= examp.AdjRand; 
    
    
    %%%%%%%%%% Severe Drift Present %%%%%%%%%%%%%%%
     
    % Shift by 30 seconds allow for 60 seconds association time
    examp.assSec =60;
    examp.drift=30;
    
    % Reset the parameters
    examp.clearCalcValues();    
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % First method, TDOA only
    simMatTDOAonly(examp);
    examp.getRand();
    Meth1_aRandDrift_60(ii)=examp.AdjRand;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Second Method, ad hoc
    examp.clearCalcValues();
    simMatadHoc(examp);
    examp.getRand();
    Meth2_aRandDrift_60(ii)=examp.AdjRand;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Third method, ideal
    examp.clearCalcValues();
    simMatIdeal(examp);
    examp.getRand();
    Meth3_aRandDrift_60(ii)= examp.AdjRand;    
    
    % Fourth model Baseline
    examp.clearCalcValues();
    toaOnlyCluster(examp);
    examp.getRand();
    Meth4_aRandDrift_60(ii)= examp.AdjRand; 
    

    n_calls(ii) = length(examp.arrivalArray);
end

%%

edges = 0.5:.05:.9;


h1 = histcounts(Meth1_aRand,edges,  'Normalization', 'probability');
h2 = histcounts(Meth2_aRand,edges, 'Normalization', 'probability');
h3 = histcounts(Meth3_aRand,edges, 'Normalization', 'probability');
h10 = histcounts(Meth4_aRand,edges, 'Normalization', 'probability');

h4 = histcounts(Meth1_aRandDrift_4,edges, 'Normalization', 'probability');
h5 = histcounts(Meth2_aRandDrift_4,edges, 'Normalization', 'probability');
h6 = histcounts(Meth3_aRandDrift_4,edges, 'Normalization', 'probability');
h11 = histcounts(Meth4_aRandDrift_4,edges, 'Normalization', 'probability');


h7 = histcounts(Meth1_aRandDrift_60,edges, 'Normalization', 'probability');
h8 = histcounts(Meth2_aRandDrift_60,edges, 'Normalization', 'probability');
h9 = histcounts(Meth3_aRandDrift_60,edges, 'Normalization', 'probability');
h12 = histcounts(Meth4_aRandDrift_60,edges, 'Normalization', 'probability');


cmap = parula(4);


figure (1)
subplot(3,1,1)
b = bar(edges(2:end), [h1; h2; h3; h10]', 'FaceColor','flat');
for k = 1:4
    b(k).CData = cmap(k,:);
end
ylim([0 .6])
ylabel('Proportion of Runs')
title('Perfect Association')
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline', 'Location','best')

subplot(3,1,2)
b = bar(edges(2:end),[h4; h5; h6; h11]', 'FaceColor','flat')
for k = 1:4
    b(k).CData = k;
end
ylim([0 .6])
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline', 'Location','best')
ylabel('Proportion of Runs')
title('Random Association 4 Sec')


subplot(3,1,3)
b = bar(edges(2:end),[h7; h8; h9; h12]', 'FaceColor','flat')
for k = 1:4
    b(k).CData = k;
end
ylim([0 .6])
xlabel('Percent of Calls Correctly Classified')
ylabel('Proportion of Runs')
title('Random Association 60 Sec')
xlabel('Adjusted Rand Index')
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline','Location','best')

subplot(3,1,2)
hold on
scatter(n_calls, Meth1_aRand, '.') 
scatter(n_calls, Meth2_aRand, '.')
scatter(n_calls, Meth3_aRand, '.')
scatter(n_calls, Meth4_aRand, '.')
title('Perfect Association')
ylabel('Adjusted Rand Index')
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline', 'Location','best')


subplot(3,2,4)
hold on
scatter(n_calls, Meth1_aRandDrift_4, '.') 
scatter(n_calls, Meth2_aRandDrift_4, '.')
scatter(n_calls, Meth3_aRandDrift_4, '.')
scatter(n_calls, Meth4_aRandDrift_4, '.')
title('Random Association 4 Sec')
ylabel('Adjusted Rand Index')
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline')

subplot(3,2,6)
hold on
scatter(n_calls, Meth1_aRandDrift_60, '.') 
scatter(n_calls, Meth2_aRandDrift_60, '.')
scatter(n_calls, Meth3_aRandDrift_60, '.')
scatter(n_calls, Meth4_aRandDrift_60, '.')
ylabel('Adjusted Rand Index')
xlabel('Number of Calls in the Dataset')
title('Random Association 60 Sec')
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline')


%% Experiment 2- Classifier performance

% Median correct classification rate for species A is 0.8, 0.85, 0.95
score_mean = [1.39 1.74 2.95];
score_mean = [1 3 6];

score_sd = [4.75 4.75 4.75];
nRuns = 200;

% Number of calls included in the analysis
n_calls = zeros(1, nRuns);


perf_meth1 = struct;
perf_meth2 = struct;
perf_meth3 = struct;
perf_methbaseline = struct;


for ii=1:nRuns
    close all
    clear spaceWhale examp
    disp([num2str(ii) ' of ' num2str(nRuns) ' runs'])
    % Create new agents
    [spaceWhale] =  createRandomSpaceWhale(0.75,6, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth);
    
    % Populate data and parameters
    examp = simulationClass();
    examp.spaceWhale = spaceWhale;
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.spaceWhale= spaceWhale;
    examp.time_cut = 10*60;
    examp.randomMiss =0;


    
    % Method 1
    examp.clearCalcValues();
    simMatTDOAonly(examp);

    examp.SppCorrRate = score_mean(1);
    examp.SppCorrSd=score_sd(1);
    perf = examp.estClassifierPerf;
    perf_meth1.Low(ii) = perf;
    examp.truthAndpreds =[];
    
    examp.SppCorrRate = score_mean(2);
    examp.SppCorrSd=score_sd(2);
    perf = examp.estClassifierPerf;
    perf_meth1.Med(ii) = perf;
    examp.truthAndpreds=[];
    
    examp.SppCorrRate = score_mean(3);
    examp.SppCorrSd=score_sd(3);
    perf = examp.estClassifierPerf;
    perf_meth1.High(ii) = perf;
    examp.truthAndpreds=[];
    

    % Second Method, ideal
    examp.clearCalcValues();
    examp.simMatIdeal();
    examp.cutoff = quantile(reshape(examp.Sim_mat,[],1),.9);
    
    examp.SppCorrRate = score_mean(1)
    examp.SppCorrSd=score_sd(1);
    perf = examp.estClassifierPerf;
    perf_meth2.Low(ii) = perf;
    examp.truthAndpreds =[];
    
    examp.SppCorrRate = score_mean(2)
    examp.SppCorrSd=score_sd(2);
    perf = examp.estClassifierPerf;
    perf_meth2.Med(ii) = perf;
    examp.truthAndpreds=[];
    
    examp.SppCorrRate = score_mean(3)
    examp.SppCorrSd=score_sd(3);
    perf = examp.estClassifierPerf;
    perf_meth2.High(ii) = perf;
    examp.truthAndpreds=[];    
    
    
    % Third Method, ad hoc
    examp.clearCalcValues();
    simMatadHoc(examp);
    examp.cutoff = quantile(reshape(examp.Sim_mat,[],1),.9);

    examp.SppCorrRate = score_mean(1);
    examp.SppCorrSd=score_sd(1);
    updateChains(examp)
    perf = examp.estClassifierPerf;
    perf_meth3.Low(ii) = perf;
    examp.truthAndpreds =[];
    
    examp.SppCorrRate = score_mean(2);
    examp.SppCorrSd=score_sd(2);
    perf = examp.estClassifierPerf;
    perf_meth3.Med(ii) = perf;
    examp.truthAndpreds=[];
    
    examp.SppCorrRate = score_mean(3);
    examp.SppCorrSd=score_sd(3);
    perf = examp.estClassifierPerf;
    perf_meth3.High(ii) = perf;
    examp.truthAndpreds=[];

    % Fourth Method, baseline
    examp.clearCalcValues();
    examp.toaOnlyCluster();
    examp.getRand();

    examp.SppCorrRate = score_mean(1);
    examp.SppCorrSd=score_sd(1);
    perf = examp.estClassifierPerf;
    perf_methbaseline.Low(ii) = perf;
    examp.truthAndpreds =[];
    
    examp.SppCorrRate = score_mean(2);
    examp.SppCorrSd=score_sd(2);
    perf = examp.estClassifierPerf;
    perf_methbaseline.Med(ii) = perf;
    examp.truthAndpreds=[];
    
    examp.SppCorrRate = score_mean(3);
    examp.SppCorrSd=score_sd(3);
    perf = examp.estClassifierPerf;
    perf_methbaseline.High(ii) = perf;
    examp.truthAndpreds=[];
    

end

perf1 = cell2mat({perf_meth2.Low([1:end]).totCorrectClassMeth})'-cell2mat({perf_meth2.Low([1:end]).totCorrect}');
perf2 = cell2mat({perf_meth2.Med([1:end]).totCorrectClassMeth})'-cell2mat({perf_meth2.Med([1:end]).totCorrect}');
perf3 = cell2mat({perf_meth2.High([1:end]).totCorrectClassMeth})'-cell2mat({perf_meth2.High([1:end]).totCorrect}');
figure
hist([perf1 perf2 perf3],20)
legend('Low', 'Medium', 'High','Location','best')



%% Make Figures
low_perf_class = table(cell2mat({perf_meth1.Low([1:end]).totCorrectClassMeth}'),...
    cell2mat({perf_meth2.Low([1:end]).totCorrectClassMeth}'),...
    cell2mat({perf_meth3.Low([1:end]).totCorrectClassMeth}'),...
    cell2mat({perf_methbaseline.Low([1:end]).totCorrectClassMeth}'));

low_perf_unaided = table(cell2mat({perf_meth1.Low([1:end]).totCorrect}'),...
    cell2mat({perf_meth2.Low([1:end]).totCorrect}'),...
    cell2mat({perf_meth3.Low([1:end]).totCorrect}'),...
    cell2mat({perf_methbaseline.Low([1:end]).totCorrect}'));


aa = table2array(low_perf_class)-table2array(low_perf_unaided);
figure (2)
subplot(3,1,1)
hist(aa)
title('Change in Classificaiton Improvement- Low Detector Performance')
legend('TDOA Only', 'Idealized Spatial Model', 'AdHoc','Baseline','Location','best')


med_perf_class = table(cell2mat({perf_meth1.Med([1:end]).totCorrectClassMeth}'),...
    cell2mat({perf_meth2.Med([1:end]).totCorrectClassMeth}'),...
    cell2mat({perf_meth3.Med([1:end]).totCorrectClassMeth}'),...
    cell2mat({perf_methbaseline.Med([1:end]).totCorrectClassMeth}'));

med_perf_unaided = table(cell2mat({perf_meth1.Med([1:end]).totCorrect}'),...
    cell2mat({perf_meth2.Med([1:end]).totCorrect}'),...
    cell2mat({perf_meth3.Med([1:end]).totCorrect}'),...
    cell2mat({perf_methbaseline.Med([1:end]).totCorrect}'));


aa = table2array(med_perf_class)-table2array(med_perf_unaided);
subplot(3,1,2)
hist(aa)

title('Change in Classificaiton Improvement- Moderate Detector Performance')
legend('TDOA Only', 'Idealized Spatial Model', 'AdHoc','Baseline','Location','best')



high_perf_class = table(cell2mat({perf_meth1.High([1:end]).totCorrectClassMeth}'),...
    cell2mat({perf_meth2.High([1:end]).totCorrectClassMeth}'),...
    cell2mat({perf_meth3.High([1:end]).totCorrectClassMeth}'),...
    cell2mat({perf_methbaseline.High([1:end]).totCorrectClassMeth}'));

high_perf_unaided = table(cell2mat({perf_meth1.High([1:end]).totCorrect}'),...
    cell2mat({perf_meth2.High([1:end]).totCorrect}'),...
    cell2mat({perf_meth3.High([1:end]).totCorrect}'),...
    cell2mat({perf_methbaseline.High([1:end]).totCorrect}'));


aa = table2array(high_perf_class)-table2array(high_perf_unaided);
subplot(3,1,3)
hist(aa)


title('Change in Classificaiton Improvement- High Detector Performance')
legend('TDOA Only', 'Idealized Spatial Model', 'AdHoc','Baseline','Location','best')
%%

edges = -.1:.02:.1;
h1 = histcounts(aa(:,1),edges,  'Normalization', 'probability');
h2 = histcounts(aa(:,2),edges, 'Normalization', 'probability');
h3 = histcounts(aa(:,3),edges, 'Normalization', 'probability');
h10 = histcounts(aa(:,4),edges, 'Normalization', 'probability');

figure (1)
b = bar(edges(2:end), [h1; h2; h3; h10]', 'FaceColor','flat');
for k = 1:4
    b(k).CData = k;
end
xlabel('Percent of Calls Correctly Classified')
ylabel('Proportion of Runs')
title('Change in Classificaiton Improvement')
xlabel('Adjusted Rand Index')
legend('TDOA Only', 'Idealized Spatial Model', 'AdHoc','Baseline','Location','best')








%% Experiment 3- Threshold Sensitivity analysis


% Run a default example first
TimeThresh =[2 3 5 10 30 60 100 200 400 600 900 1200];

nIters =200;



% Get some reasonable thresholds
[spaceWhale] =  createRandomSpaceWhale(.75,n_agents, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth);
    
    % Populate data and parameters
    examp = simulationClass();
    examp.spaceWhale = spaceWhale;
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.spaceWhale= spaceWhale;
    examp.time_cut = 1500;
    examp.randomMiss =0;
    
% First method, tdoa only
simMatTDOAonly(examp);
SimThresh1 = quantile(reshape(examp.Sim_mat,[],1), [0:.02:1]);

% Second method, Ideal
examp.clearCalcValues();
examp.simMatIdeal();
SimThresh2 = quantile(reshape(examp.Sim_mat,[],1), [0:.02:1]);


% Third method, ad hoc
examp.clearCalcValues();
examp.simMatadHoc();
SimThresh3 = quantile(reshape(examp.Sim_mat,[],1), [0:.02:1]);


ExpScoresMeth1 = zeros(nIters, length(TimeThresh), length(SimThresh1));
ExpScoresMeth2 = ExpScoresMeth1;
ExpScoresMeth3 = ExpScoresMeth1;
ExpScoresBaseline = ExpScoresMeth1;



for iter=1:nIters
    
    disp(num2str(iter))
    clear spaceWhale examp
    close all
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
    
    
    for jj = 1:length(SimThresh1)
        examp.Cluster_id =[];
        examp.cutoff = SimThresh1(jj);
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=exp(TimeThresh(ii));
            examp.updateChains;
            examp.getRand();
            ExpScoresMeth1(iter, ii,jj) = examp.AdjRand;
            
        end
        
    end
    
    figure
    imagesc(squeeze(ExpScoresMeth1(1,:,:)))
    title('Method 1')
    
    
    % Second method, Ideal
    examp.clearCalcValues();
    examp.time_cut = 1500;
    examp.simMatIdeal();
    
    
    for jj = 1:length(SimThresh2)
        examp.Cluster_id =[];
        examp.cutoff = SimThresh2(jj);
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=exp(TimeThresh(ii));
            examp.updateChains;
            examp.getRand();
            ExpScoresMeth2(iter, ii,jj) = examp.AdjRand;
            
        end
        
    end
    
    figure
    imagesc(squeeze(ExpScoresMeth2(1,:,:)))
    title('Method 2')
    
    % Third method, ad hoc
    examp.clearCalcValues();
    examp.time_cut = 1500;
    examp.simMatadHoc();
    
    
    for jj = 1:length(SimThresh3)
        examp.Cluster_id =[];
        examp.cutoff = SimThresh3(jj);
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=exp(TimeThresh(ii));
            examp.updateChains;
            examp.getRand();
            ExpScoresMeth3(iter, ii,jj) = examp.AdjRand;
            
        end
        
    end
    
    figure
    imagesc(squeeze(ExpScoresMeth3(1,:,:)))
    title('Method 3')
    
    
    % Fourth method, baseline
    examp.clearCalcValues();
    examp.toaOnlyCluster();
    examp.getRand();
    
    for jj = 1:length(SimThresh3)
        examp.Cluster_id =[];
        examp.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=exp(TimeThresh(ii));
             examp.toaOnlyCluster();
            examp.getRand();
            ExpScoresBaseline(iter, ii,jj) = examp.AdjRand;
            
        end
        
    end
    
    
    
end



% 
% imagesc(ExpScores)
% colorbar
% 
% xticklabs =(num2str(round(SimThresh(1:4:end),2)));
% 
% xticks((1:4:length(SimThresh)))
% xticklabels({round(SimThresh(1:4:end),2)})
% yticks((1:2:length(TimeThresh)))
% yticklabels({round(TimeThresh(1:2:end),2)})

%%
% First method, TDOA only
    
    stClassifierPerf(examp)
    Meth1_aRandDrift_60(ii)=examp.AdjRand;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Second Method, ad hoc
    examp.clearCalcValues();
    simMatadHoc(examp);
    examp.getRand();
    Meth2_aRandDrift_60(ii)=examp.AdjRand;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Third method, ideal
    examp.clearCalcValues();
    simMatIdeal(examp);
    examp.getRand();
    Meth3_aRandDrift_60(ii)= examp.AdjRand;    
    
    % Fourth model Baseline
    examp.clearCalcValues();
    toaOnlyCluster(examp);
    examp.getRand();
    Meth4_aRandDrift_60(ii)= examp.AdjRand; 


