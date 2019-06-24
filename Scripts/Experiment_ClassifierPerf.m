%% Experiment- Classifier performance
% Setup for the classifier performance experiment. Simulation looks at how
% classifier improvment given current (binary) classifier performance

% Run experiments
close all; clear all
clear classes; clc

whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2})
load(whereAmI{3})


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
clear whereAmI
%%
% Median correct classification rate for species A is 0.8, 0.85, 0.95
%score_mean = [1.39 1.74 2.95];
score_mean = [1 3 6];


score_sd = [4.75 4.75 4.75];
nRuns = 100;

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
    [spaceWhale] =  createRandomSpaceWhale(0.75,7, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth);
    
    % Populate data and parameters
    examp = simulationClass();
    examp.spaceWhale = spaceWhale;
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.spaceWhale= spaceWhale;
    examp.time_cut = 10*60;
    examp.randomMiss =0;

    % Method 1- TDOA only
    examp.clearCalcValues();
    simMatTDOAonly(examp);
    examp.cutoff = .5;
    examp.time_cut = 490;

    examp.SppCorrRate = score_mean(1);
    examp.SppCorrSd=score_sd(1);
    perf = examp.estClassifierPerf;
    perf_meth1(ii).Low = perf;
    examp.truthAndpreds =[];
    
    
    examp.SppCorrRate = score_mean(2);
    examp.SppCorrSd=score_sd(2);
    perf = examp.estClassifierPerf;
    perf_meth1(ii).Med = perf;
    examp.truthAndpreds=[];
    
    examp.SppCorrRate = score_mean(3);
    examp.SppCorrSd=score_sd(3);
    perf = examp.estClassifierPerf;
    perf_meth1(ii).High = perf;
    examp.truthAndpreds=[];
    

    % Second Method, ideal
    examp.clearCalcValues();
    examp.simMatIdeal();
    examp.cutoff = .49;
    examp.time_cut = 1020;
    
    examp.SppCorrRate = score_mean(1)
    examp.SppCorrSd=score_sd(1);
    perf = examp.estClassifierPerf;
    perf_meth2(ii).Low = perf;
    examp.truthAndpreds =[];
    
    examp.SppCorrRate = score_mean(2)
    examp.SppCorrSd=score_sd(2);
    perf = examp.estClassifierPerf;
    perf_meth2(ii).Med = perf;
    examp.truthAndpreds=[];
    
    examp.SppCorrRate = score_mean(3)
    examp.SppCorrSd=score_sd(3);
    perf = examp.estClassifierPerf;
    perf_meth2(ii).High = perf;
    examp.truthAndpreds=[];    
    
    
    % Third Method, ad hoc
    examp.clearCalcValues();
    simMatadHoc(examp);
    examp.cutoff = .15;
    examp.time_cut = 677;

    examp.SppCorrRate = score_mean(1);
    examp.SppCorrSd=score_sd(1);
    updateChains(examp)
    perf = examp.estClassifierPerf;
    perf_meth3(ii).Low = perf;
    examp.truthAndpreds =[];
    
    examp.SppCorrRate = score_mean(2);chil
    examp.SppCorrSd=score_sd(2);
    perf = examp.estClassifierPerf;
    perf_meth3(ii).Med = perf;
    examp.truthAndpreds=[];
    
    examp.SppCorrRate = score_mean(3);
    examp.SppCorrSd=score_sd(3);
    perf = examp.estClassifierPerf;
    perf_meth3(ii).High = perf;
    examp.truthAndpreds=[];

    % Fourth Method, baseline
    examp.clearCalcValues();
    examp.time_cut = 81;
    examp.toaOnlyCluster();
    examp.getRand();

    examp.SppCorrRate = score_mean(1);
    examp.SppCorrSd=score_sd(1);
    perf = examp.estClassifierPerf;
    perf_methbaseline(ii).Low = perf;
    examp.truthAndpreds =[];
    
    examp.SppCorrRate = score_mean(2);
    examp.SppCorrSd=score_sd(2);
    perf = examp.estClassifierPerf;
    perf_methbaseline(ii).Med = perf;
    examp.truthAndpreds=[];
    
    examp.SppCorrRate = score_mean(3);
    examp.SppCorrSd=score_sd(3);
    perf = examp.estClassifierPerf;
    perf_methbaseline(ii).High = perf;
    examp.truthAndpreds=[];
    

end

% 
% perf1 = cell2mat({perf_meth2.Low([1:end]).totCorrectClassMeth})'-cell2mat({perf_meth2.Low([1:end]).totCorrect}');
% perf2 = cell2mat({perf_meth2.Med([1:end]).totCorrectClassMeth})'-cell2mat({perf_meth2.Med([1:end]).totCorrect}');
% perf3 = cell2mat({perf_meth2.High([1:end]).totCorrectClassMeth})'-cell2mat({perf_meth2.High([1:end]).totCorrect}');
% figure
% hist([perf1 perf2 perf3],20)
% legend('Low', 'Medium', 'High','Location','best')
% 


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
