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
score_mean = [1.39 1.74 2.95];
betaParm1= [8, 10, 10]
% betaParm2=[1.5, 1,.6]    
% figure;
% for ii=1:3
% 
% 
% subplot(3,1,ii)
% hist(betarnd(betaParm1(ii),betaParm2(ii),[1 1000]), 20)
% 
% end


nRuns = 20;



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
    [spaceWhale] =  createRandomSpaceWhale(0.75,8, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth, [array_struct.master, array_struct.slave([2,3,4])]);
    
    % Populate data and parameters
    examp = simulationClass();
    examp.spaceWhale = spaceWhale;
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.spaceWhale= spaceWhale;
    examp.randomMiss =0;
    examp.child_idx = [2,3,4];

    % Method 1- TDOA only
    examp.clearCalcValues();
    examp.cutoff = .75;
    examp.time_cut = 5*60;
    simMatTDOAonly(examp);
    examp.updateClusterID;



    examp.betaParm1 = betaParm1(1);
    examp.betaParm2=betaParm2(1);
    perf = examp.estClassifierPerf;
    perf_meth1(ii).Low = perf;

    
    examp.betaParm1 = betaParm1(2);
    examp.betaParm2=betaParm2(2);
    perf = examp.estClassifierPerf;
    perf_meth1(ii).Med = perf;

    
    examp.betaParm1 = betaParm1(3);
    examp.betaParm2=betaParm2(3);
    perf = examp.estClassifierPerf;
    perf_meth1(ii).High = perf;

    
    % Second Method, ideal
    examp.clearCalcValues();
    examp.simMatIdeal();
    examp.cutoff = .5;
    examp.time_cut = 5*60;
    
    examp.betaParm1 = betaParm1(1);
    examp.betaParm2=betaParm2(1);
    perf = examp.estClassifierPerf;
    perf_meth2(ii).Low = perf;
  
    
    examp.betaParm1 = betaParm1(2);
    examp.betaParm2=betaParm2(2);
    perf = examp.estClassifierPerf;
    perf_meth2(ii).Med = perf;

    examp.betaParm1 = betaParm1(3);
    examp.betaParm2=betaParm2(3);
    perf = examp.estClassifierPerf;
    perf_meth2(ii).High = perf;

    
    
    
    % Third Method, ad hoc
    examp.clearCalcValues();
    simMatadHoc(examp);
    examp.cutoff = .5;
    examp.time_cut = 5*60;

    examp.betaParm1 = betaParm1(1);
    examp.betaParm2=betaParm2(1);
    updateClusterID(examp)
    perf = examp.estClassifierPerf;
    perf_meth3(ii).Low = perf;

    examp.betaParm1 = betaParm1(2);
    examp.betaParm2=betaParm2(2);
    perf = examp.estClassifierPerf;
    perf_meth3(ii).Med = perf;

    
    examp.betaParm1 = betaParm1(3);
    examp.betaParm2=betaParm2(3);
    perf = examp.estClassifierPerf;
    perf_meth3(ii).High = perf;

    
    
    
    % Fourth Method, baseline
    examp.clearCalcValues();
    examp.time_cut = 60;
    examp.toaOnlyCluster();
    examp.getRand();

    examp.betaParm1 = betaParm1(1);
    examp.betaParm2=betaParm2(1);
    perf = examp.estClassifierPerf;
    perf_methbaseline(ii).Low = perf;

    examp.betaParm1 = betaParm1(2);
    examp.betaParm2=betaParm2(2);
    perf = examp.estClassifierPerf;
    perf_methbaseline(ii).Med = perf;

    examp.betaParm1 = betaParm1(3);
    examp.betaParm2=betaParm2(3);
    perf = examp.estClassifierPerf;
    perf_methbaseline(ii).High = perf;

    
    
    
end



%% Compare Results
correct_classMeth1 =table(zeros(nRuns,1), zeros(nRuns,1), zeros(nRuns,1),zeros(nRuns,1), 'VariableNames', {'Low', 'Medium', 'High', 'Method'});
correct_classMeth4 =correct_classMeth1;
correct_classMeth3 = correct_classMeth1;
correct_classMeth2 = correct_classMeth1;

for ii=1:length(perf_meth1)
    mod = struct2table(perf_meth1(ii).Low);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth1.Low(ii) = prct_improvement;
    
    
    mod = struct2table(perf_meth1(ii).Med);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth1.Medium(ii) = prct_improvement;
    
    % 1s correctly classified as 1
    mod = struct2table(perf_meth1(ii).High);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth1.High(ii) = prct_improvement;
    correct_classMeth1.Method(ii) =1;
end



for ii=1:length(perf_meth2)
    
    mod = struct2table(perf_meth2(ii).Low);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth2.Low(ii) = prct_improvement;
    
    % 1s correctly classified as 1
    mod = struct2table(perf_meth2(ii).Med);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth2.Medium(ii) = prct_improvement;
    
        % 1s correctly classified as 1
    mod = struct2table(perf_meth2(ii).High);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth2.High(ii) = prct_improvement;
    correct_classMeth2.Method(ii) =2;
end





for ii=1:length(perf_meth3)
    mod = struct2table(perf_meth3(ii).Low);    
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth3.Low(ii) = prct_improvement;
    
    
    % 1s correctly classified as 1
    mod = struct2table(perf_meth3(ii).Med);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth3.Medium(ii) = prct_improvement;
    
    mod = struct2table(perf_meth3(ii).High);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth3.High(ii) = prct_improvement;
    correct_classMeth3.Method(ii) = 3;
end



for ii=1:length(perf_methbaseline)
    mod = struct2table(perf_methbaseline(ii).Low);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth4.Low(ii) = prct_improvement;
    
    
    % 1s correctly classified as 1
    mod = struct2table(perf_methbaseline(ii).Med);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) + length(find(mod.ClusterScore<0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1)) + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
    correct_classMeth4.Medium(ii) = prct_improvement;
    
    
    % 1s correctly classified as 1
    mod = struct2table(perf_methbaseline(ii).High);
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
title('Low Performance Classifier')
ylabel('Percent Improvement')
ylim([-60 80])


subplot(3,1,2)
Y = (aa.Medium*100);
Y = reshape(Y, [], 4)
violin(Y,'xlabel',xaxtitle, 'edgecolor', [], 'mc', [], 'quan', '--k', 'plotlegend',[])
fix_xticklabels(gca,0.1,{'FontSize',9});
textvals = strcat(repmat(['Median = '],4,1),(num2str(round(median(Y),2)')));
text([1:4]-.17', ones(1,4)*70, textvals)
title('Mid Performance Classifier')
ylabel('Percent Improvement')
ylim([-60 80])


subplot(3,1,3)
Y = (aa.High*100);
Y = reshape(Y, [], 4)
violin(Y,'xlabel',xaxtitle, 'edgecolor', [], 'mc', [], 'quan', '--k', 'plotlegend',[])
fix_xticklabels(gca,0.1,{'FontSize',9});
textvals = strcat(repmat(['Median = '],4,1),(num2str(round(median(Y),2)')));
text([1:4]-.17', ones(1,4)*70, textvals)
title('High Performance Classifier')
ylabel('Percent Improvement')
ylim([-60 80])


