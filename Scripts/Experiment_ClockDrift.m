
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
n_hrs = 0.75;
n_agents = 7;
nRuns = 100;

% Number of calls included in the analysis
n_calls = zeros(1, nRuns);

Meth1_aRand = zeros(1, nRuns);
Meth2_aRand = zeros(1, nRuns);
Meth3_aRand = zeros(1, nRuns);
Meth4_aRand = zeros(1, nRuns);


Meth1_aRandDrift_2 = zeros(1, nRuns);
Meth2_aRandDrift_2 = zeros(1, nRuns);
Meth3_aRandDrift_2 = zeros(1, nRuns);
Meth4_aRandDrift_2 = zeros(1, nRuns);


Meth1_aRandDrift_10 = zeros(1, nRuns);
Meth2_aRandDrift_10 = zeros(1, nRuns);
Meth3_aRandDrift_10 = zeros(1, nRuns);
Meth4_aRandDrift_10 = zeros(1, nRuns);


for ii =1:nRuns
    disp(num2str(ii))

    clear spaceWhale examp
    close all
%     
    spaceWhale=[];
    % Create new agents
    [spaceWhale] =  createRandomSpaceWhale(n_hrs,n_agents, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth, [5, 1,2]);
    
    % Populate data and parameters
    examp = simulationClass();
    examp.spaceWhale = spaceWhale;
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.spaceWhale= spaceWhale;
    examp.time_cut = 7*60;
    examp.randomMiss =0;
    
    jj =1;
    
    disp(['Run ', num2str(ii), ' exp ',num2str(jj), ' of 12'])
    
    % First method, TDOA only
    examp.time_cut = 7*60;
    examp.cutoff =0.95;
    simMatTDOAonly(examp);
    examp.getRand();
    Meth1_aRand(ii)=examp.AdjRand;
    

    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    
    % Third Method, ad hoc
    examp.clearCalcValues();
    examp.cutoff =0.95;
    examp.time_cut = 7*60;
    simMatadHoc(examp);
    examp.getRand();
    Meth3_aRand(ii)=examp.AdjRand;
    
    
   jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    
    % Second method, ideal
    examp.clearCalcValues();
    examp.cutoff =0.6;
    examp.time_cut = 7*60;
    simMatIdeal(examp);
    examp.getRand();
    Meth2_aRand(ii)= examp.AdjRand;
    
    
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
    
    

    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % First method, TDOA only
    examp.clearCalcValues();
    examp.cutoff =0.95;
    examp.time_cut = 7*60;
    simMatTDOAonly(examp);
    examp.getRand();
    Meth1_aRandDrift_2(ii)=examp.AdjRand;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Third Method, ad hoc
    examp.clearCalcValues();
    examp.cutoff =0.95;
    examp.time_cut = 7*60;
    simMatadHoc(examp);
    examp.getRand();
    Meth3_aRandDrift_2(ii)=examp.AdjRand;

    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Second method, ideal
    examp.clearCalcValues();
    examp.cutoff =0.6;
    examp.time_cut = 7*60;
    simMatIdeal(examp);
    examp.getRand();
    Meth2_aRandDrift_2(ii)= examp.AdjRand;  
    
    % Fourth model Baseline
    examp.clearCalcValues();
    toaOnlyCluster(examp);
    examp.getRand();
    Meth4_aRandDrift_2(ii)= examp.AdjRand; 
    
    
    %%%%%%%%%% Severe Drift Present %%%%%%%%%%%%%%%
     
    % Shift by 10 seconds allow for 15 seconds association time
    examp.assSec =15;
    examp.drift=10;
      
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % First method, TDOA only
    simMatTDOAonly(examp);
    examp.cutoff =0.95;
    examp.time_cut = 7*60;
    examp.getRand();
    Meth1_aRandDrift_10(ii)=examp.AdjRand;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Third Method, ad hoc
    examp.clearCalcValues();
    examp.cutoff =0.95;
    examp.time_cut = 7*60;
    simMatadHoc(examp);
    examp.getRand();
    Meth3_aRandDrift_10(ii)=examp.AdjRand;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Second method, ideal
    examp.clearCalcValues();
    examp.cutoff =0.6;
    examp.time_cut = 7*60;
    simMatIdeal(examp);
    examp.getRand();
    Meth2_aRandDrift_10(ii)= examp.AdjRand;    
    
    % Fourth model Baseline
    examp.clearCalcValues();
    toaOnlyCluster(examp);
    examp.getRand();
    Meth4_aRandDrift_10(ii)= examp.AdjRand; 
    

    n_calls(ii) = length(examp.arrivalArray);
end

%%

edges = 0.5:.05:.9;


h1 = histcounts(Meth1_aRand,edges,  'Normalization', 'probability');
h2 = histcounts(Meth2_aRand,edges, 'Normalization', 'probability');
h3 = histcounts(Meth3_aRand,edges, 'Normalization', 'probability');
h10 = histcounts(Meth4_aRand,edges, 'Normalization', 'probability');

h4 = histcounts(Meth1_aRandDrift_2,edges, 'Normalization', 'probability');
h5 = histcounts(Meth2_aRandDrift_2,edges, 'Normalization', 'probability');
h6 = histcounts(Meth3_aRandDrift_2,edges, 'Normalization', 'probability');
h11 = histcounts(Meth4_aRandDrift_2,edges, 'Normalization', 'probability');


h7 = histcounts(Meth1_aRandDrift_10,edges, 'Normalization', 'probability');
h8 = histcounts(Meth2_aRandDrift_10,edges, 'Normalization', 'probability');
h9 = histcounts(Meth3_aRandDrift_10,edges, 'Normalization', 'probability');
h12 = histcounts(Meth4_aRandDrift_10,edges, 'Normalization', 'probability');


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
scatter(n_calls, Meth1_aRandDrift_2, '.') 
scatter(n_calls, Meth2_aRandDrift_2, '.')
scatter(n_calls, Meth3_aRandDrift_2, '.')
scatter(n_calls, Meth4_aRandDrift_2, '.')
title('Random Association 4 Sec')
ylabel('Adjusted Rand Index')
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline')

subplot(3,2,6)
hold on
scatter(n_calls, Meth1_aRandDrift_10, '.') 
scatter(n_calls, Meth2_aRandDrift_10, '.')
scatter(n_calls, Meth3_aRandDrift_10, '.')
scatter(n_calls, Meth4_aRandDrift_10, '.')
ylabel('Adjusted Rand Index')
xlabel('Number of Calls in the Dataset')
title('Random Association 60 Sec')
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline')
