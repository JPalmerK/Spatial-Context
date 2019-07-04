
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



%% Experiment 1 Ideal method assuming no drift

% Number of hours the experiment runs and the number of agents included
n_hrs = 0.75;
n_agents = 6;
nRuns = 20;

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


perf_meth1 = struct;
perf_meth2 = struct;
perf_meth3 = struct;
perf_methbaseline = struct;


%%

for ii =11:20
    disp(num2str(ii))

    clear spaceWhale examp
    close all
%     
    spaceWhale=[];
    % Create new agents
    [spaceWhale] =  createRandomSpaceWhale(n_hrs,n_agents, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth, array_struct.slave([1,2,3]));
    
    % Populate data and parameters
    examp = simulationClass();
    examp.spaceWhale = spaceWhale;
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.spaceWhale= spaceWhale;
    examp.cutoff = .95;
    examp.time_cut = 5*60;
    examp.randomMiss =0;
    examp.child_idx = [1,2,3];

    examp.betaParm1 = 20;
    examp.betaParm2 = 17;
    
    jj =1;
    
    disp(['Run ', num2str(ii), ' exp ',num2str(jj), ' of 12'])
    
    % First method, TDOA only
    simMatTDOAonly(examp);
    examp.getRand();
    examp.maxEltTime =  quantile(diff(examp.arrivalArray(:,1)), .85);
    Meth1_aRand(ii)=examp.AdjRand;
    perf = examp.estClassifierPerf;
    perf_meth1(ii).ideal = perf;
    

    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    
    % Third Method, ad hoc
    examp.clearCalcValues();
    simMatadHoc(examp);
    examp.getRand();
    Meth3_aRand(ii)=examp.AdjRand;
    perf = examp.estClassifierPerf;
    perf_meth3(ii).ideal = perf;
    
    
   jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    
    % Second method, ideal
    examp.clearCalcValues();
    simMatIdeal(examp);
    examp.getRand();
    Meth2_aRand(ii)= examp.AdjRand;
    perf = examp.estClassifierPerf;
    perf_meth2(ii).ideal = perf;
    
    
    %Fourth model Baseline
    examp.clearCalcValues();
    toaOnlyCluster(examp);
    examp.getRand();
    Meth4_aRand(ii)= examp.AdjRand;
    perf = examp.estClassifierPerf;
    perf_methbaseline(ii).ideal = perf;

    %%%%%%%%%%%% Slight Drift Present %%%%%%%%%%%%%%%
    

    % Shift by two seconds allow for 4 seconds of drift (default)
    examp.randomMiss =1;
    examp.assSec = 4;
    examp.drift=2;
    
    

    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % First method, TDOA only
    examp.clearCalcValues();
    simMatTDOAonly(examp);
    examp.getRand();
    Meth1_aRandDrift_2(ii)=examp.AdjRand;
    perf = examp.estClassifierPerf;
    perf_meth1(ii).drift4 = perf;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Third Method, ad hoc
    examp.clearCalcValues();
    simMatadHoc(examp);
    examp.getRand();
    Meth3_aRandDrift_2(ii)=examp.AdjRand;
    perf = examp.estClassifierPerf;
    perf_meth3(ii).drift4 = perf;
    
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Second method, ideal
    examp.clearCalcValues();
    simMatIdeal(examp);
    examp.getRand();
    Meth2_aRandDrift_2(ii)= examp.AdjRand;  
    perf = examp.estClassifierPerf;
    perf_meth2(ii).drift4 = perf;
    
    %Fourth model Baseline
    examp.clearCalcValues();
    toaOnlyCluster(examp);
    examp.getRand();
    Meth4_aRandDrift_2(ii)= examp.AdjRand; 
    perf = examp.estClassifierPerf;
    perf_methbaseline(ii).drift4 = perf;   
    
    %%%%%%%%%% Severe Drift Present %%%%%%%%%%%%%%%
     
    % Shift by 10 seconds allow for 15 seconds association time
    examp.assSec =7;
    examp.drift=5;
      
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % First method, TDOA only
    simMatTDOAonly(examp);
    examp.getRand();
    Meth1_aRandDrift_10(ii)=examp.AdjRand;
    perf = examp.estClassifierPerf;
    perf_meth1(ii).drift10 = perf;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Third Method, ad hoc
    examp.clearCalcValues();
    examp.cutoff =0.5;
    simMatadHoc(examp);
    examp.getRand();
    Meth3_aRandDrift_10(ii)=examp.AdjRand;
    perf = examp.estClassifierPerf;
    perf_meth3(ii).drift10 = perf;
    
    jj =jj+1;
    disp(['Run ', num2str(ii), ' exp ',num2str(jj) , ' of 12'])
    
    % Second method, ideal
    examp.clearCalcValues();
    simMatIdeal(examp);
    examp.getRand();
    Meth2_aRandDrift_10(ii)= examp.AdjRand;    
    perf = examp.estClassifierPerf;
    perf_meth2(ii).drift10 = perf;
    
    % Fourth model Baseline
    examp.clearCalcValues();
    toaOnlyCluster(examp);
    examp.getRand();
    Meth4_aRandDrift_10(ii)= examp.AdjRand; 
    perf = examp.estClassifierPerf;
    perf_methbaseline(ii).drift10 = perf;   

    n_calls(ii) = length(examp.arrivalArray);
end

%%

edges = 0:.1:.9;
h1 = histcounts(Meth1_aRand,edges,  'Normalization', 'probability');
h2 = histcounts(Meth2_aRand,edges, 'Normalization', 'probability');
h3 = histcounts(Meth3_aRand,edges, 'Normalization', 'probability');
h10 = histcounts(Meth4_aRand,edges, 'Normalization', 'probability');

edges1 = 0:.1:.4
h4 = histcounts(Meth1_aRandDrift_2,edges1, 'Normalization', 'probability');
h5 = histcounts(Meth2_aRandDrift_2,edges1, 'Normalization', 'probability');
h6 = histcounts(Meth3_aRandDrift_2,edges1, 'Normalization', 'probability');
h11 = histcounts(Meth4_aRandDrift_2,edges1, 'Normalization', 'probability');


h7 = histcounts(Meth1_aRandDrift_10,edges1, 'Normalization', 'probability');
h8 = histcounts(Meth2_aRandDrift_10,edges1, 'Normalization', 'probability');
h9 = histcounts(Meth3_aRandDrift_10,edges1, 'Normalization', 'probability');
h12 = histcounts(Meth4_aRandDrift_10,edges1, 'Normalization', 'probability');


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
b = bar(edges1(2:end),[h4; h5; h6; h11]', 'FaceColor','flat')
for k = 1:4
    b(k).CData = k;
end
ylim([0 .6])
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline', 'Location','best')
ylabel('Proportion of Runs')
title('Random Association 4 Sec')


subplot(3,1,3)
b = bar(edges1(2:end),[h7; h8; h9; h12]', 'FaceColor','flat')
for k = 1:4
    b(k).CData = k;
end
xlabel('Percent of Calls Correctly Classified')
ylabel('Proportion of Runs')
title('Random Association 60 Sec')
xlabel('Adjusted Rand Index')
legend('TDOA Only','AdHoc', 'Idealized Spatial Model', 'Baseline','Location','best')
ylim([0 .6])

figure;
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

%% Classifier Performance


correct_classMeth1 =table(zeros(nRuns,1), zeros(nRuns,1), zeros(nRuns,1),zeros(nRuns,1), 'VariableNames', {'Low', 'Medium', 'High', 'Method'});
correct_classMeth4 =correct_classMeth1;
correct_classMeth3 = correct_classMeth1;
correct_classMeth2 = correct_classMeth1;

for ii=1:length(perf_meth1)
    mod = struct2table(perf_meth1(ii).ideal);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth1.Low(ii) = prct_improvement;
    
    
    mod = struct2table(perf_meth1(ii).drift4  );
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth1.Medium(ii) = prct_improvement;
    
    % 1s correctly classified as 1
    mod = struct2table(perf_meth1(ii).drift10);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth1.High(ii) = prct_improvement;
    correct_classMeth1.Method(ii) =1;
end



for ii=1:length(perf_meth2)
    
    mod = struct2table(perf_meth2(ii).ideal);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth2.Low(ii) = prct_improvement;
    
    % 1s correctly classified as 1
    mod = struct2table(perf_meth2(ii).drift4);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth2.Medium(ii) = prct_improvement;
    
        % 1s correctly classified as 1
    mod = struct2table(perf_meth2(ii).drift10);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth2.High(ii) = prct_improvement;
    correct_classMeth2.Method(ii) =2;
end





for ii=1:length(perf_meth3)
    mod = struct2table(perf_meth3(ii).ideal);    
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth3.Low(ii) = prct_improvement;
    
    
    % 1s correctly classified as 1
    mod = struct2table(perf_meth3(ii).drift4);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth3.Medium(ii) = prct_improvement;
    
    mod = struct2table(perf_meth3(ii).drift10);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth3.High(ii) = prct_improvement;
    correct_classMeth3.Method(ii) = 3;
end



for ii=1:length(perf_methbaseline)
    mod = struct2table(perf_methbaseline(ii).ideal);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth4.Low(ii) = prct_improvement;
    
    
    % 1s correctly classified as 1
    mod = struct2table(perf_methbaseline(ii).drift4);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth4.Medium(ii) = prct_improvement;
    
    
    % 1s correctly classified as 1
    mod = struct2table(perf_methbaseline(ii).drift10);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth4.High(ii) = prct_improvement;
    correct_classMeth4.Method(ii) =4;
end





aa =[correct_classMeth1; correct_classMeth2; correct_classMeth3; correct_classMeth4];
aa.Method(aa.Method==0)=2

xaxtitle = {['TDOA'], ['Idealize Spatial'], ['AdHoc Spatial'], ['Baseline']};


figure
subplot(3,1,1)
Y = (aa.Low*100);
Y = reshape(Y, [], 4)
violin(Y,'xlabel',xaxtitle, 'edgecolor', [], 'mc', [], 'quan', '--k', 'plotlegend',[])
fix_xticklabels(gca,0.1,{'FontSize',9});
textvals = strcat(repmat(['Median = '],4,1),(num2str(round(median(Y),2)')));
text([1:4]-.17', ones(1,4)*70, textvals)
title('No Clock Drift')
ylabel('Percent Improvement')
ylim([-60 80])


subplot(3,1,2)
Y = (aa.Medium*100);
Y = reshape(Y, [], 4)
violin(Y,'xlabel',xaxtitle, 'edgecolor', [], 'mc', [], 'quan', '--k', 'plotlegend',[])
fix_xticklabels(gca,0.1,{'FontSize',9});
textvals = strcat(repmat(['Median = '],4,1),(num2str(round(median(Y),2)')));
text([1:4]-.17', ones(1,4)*70, textvals)
title('2s  Drift, 4s Random Association')
ylabel('Percent Improvement')
ylim([-60 80])


subplot(3,1,3)
Y = (aa.High*100);
Y = reshape(Y, [], 4)
violin(Y,'xlabel',xaxtitle, 'edgecolor', [], 'mc', [], 'quan', '--k', 'plotlegend',[])
fix_xticklabels(gca,0.1,{'FontSize',9});
textvals = strcat(repmat(['Median = '],4,1),(num2str(round(median(Y),2)')));
text([1:4]-.17', ones(1,4)*70, textvals)
title('10s  Drift, 20s Random Association')
ylabel('Percent Improvement')
ylim([-60 80])
























