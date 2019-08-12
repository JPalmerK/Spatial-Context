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

child_idx = [1 2 3];
%% Create the agents



[spaceWhale] =   createRandomSpaceWhale(0.75, 6, hyd_arr,...
    array_struct,hydrophone_struct, ssp, grid_depth,...
    [array_struct.master, array_struct.slave(child_idx)]);

simStruct=struct();
simStruct.spaceWhale=spaceWhale;
simStruct.array_struct=array_struct;
simStruct.truncateKm=10;
simStruct.s = 8
simStruct.PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
% Position/soundspeed uncertainty, see Eva's paper
% receiver position = 5 m -> 5/1500 = 0.0004 sec
% height LSQ peak (check with Tyler) = .1 sec
% sound speed profile = 50 m/s -> 50/1500 =0.3 sec
simStruct.drift=0;
simStruct.maxEltTime =60;
simStruct.cutoff = .25;
simStruct.child_idx= child_idx;
simStruct.randomMiss=0;
simStruct.c =1500;
simStruct.arrivalTable = UpdateArrTable(simStruct);
simStruct.arrivalArray= UpdateArrArray(simStruct);
simStruct.TDOA_vals = UpdateTDOA(simStruct);




%% Sensitivity Experiment
TimeThresh=linspace(1, 30,20);
SimThresh = linspace(0,1,30);
SimThresh1 =linspace(0,1,30);

ExpScoresMeth_out2D = zeros(length(TimeThresh), length(SimThresh), 20)/0;
ExpScoresMeth_outTDOA = ExpScoresMeth_out2D;
ExpScoresMeth_outMaxMean = ExpScoresMeth_out2D;

allSimStructs =[];


tic

parfor ii=1:20
    % Replace the space whale component
    [spaceWhale] =   createRandomSpaceWhale(1, 6, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth,...
        [array_struct.master, array_struct.slave(child_idx)]);
    
    % Copy for each of the methods
    simStructNew = simStruct;
    
    
    % Update the arrival array and simulation matrix
    simStructNew.spaceWhale=spaceWhale;
    simStructNew.arrivalTable = UpdateArrTable(simStructNew);
    simStructNew.arrivalArray= UpdateArrArray(simStructNew);
    simStructNew.TDOA_vals = UpdateTDOA(simStructNew);
    
    
%     % Create copy for TDOA method
%     simStructTDOA = simStructNew;
%     
%     % TDOA only method
%     simStructTDOA.Sim_mat = simMatTDOAonly(simStructTDOA);
%     
%     %Run the sensitivity loop
%     [ExpScoresMethTDOA, ~] = runSensitivtyLp(simStructTDOA,TimeThresh,SimThresh);
%     ExpScoresMeth_outTDOA(:,:,ii) = ExpScoresMethTDOA;
%     
%     
%     % 2D XCORR METHOD
%     simStructNew.Sim_mat =simMatIdealXcorrDist(simStructNew);
%     
%     % Run the sensitivity loop
%     [ExpScoresMeth2D nAgents] = runSensitivtyLp(simStructNew,TimeThresh,SimThresh1);
%     ExpScoresMeth_out2D(:,:,ii) = ExpScoresMeth2D;
%     
      % Max of mean 
    simStructMaxMean = simStructNew;
    
    % Mean method
    %simStructMaxMean.Sim_mat = simMatMaxofMean(simStructMaxMean);
    
    % Product Method
    simStructMaxMean.Sim_mat = simMatMaxofProd(simStructMaxMean)
    
    
    %Run the sensitivity loop
    [ExpScoresMethMaxMean, ~] = runSensitivtyLp(simStructMaxMean,TimeThresh,SimThresh);
    ExpScoresMeth_outMaxMean(:,:,ii) = ExpScoresMethMaxMean;
    
    ii
   
end
toc


figure(1)
% subplot(3,1,1)
% imagesc(SimThresh1nanmedian(ExpScoresMeth_out2D,3)) ,  axis xy, colorbar
% title('Median ARI Max of Mean')
% 
% 
% subplot(3,1,2)
% imagesc(TimeThresh,SimThresh, nanmedian(ExpScoresMeth_outTDOA,3)) ,  axis xy, colorbar
% title('TDOA only')
% 


subplot(3,1,3)
imagesc(TimeThresh,SimThresh1, nanmean(ExpScoresMeth_outMaxMean,3)) ,  axis xy, colorbar
title('TDOA only')
xlabel('Time Threshold (s)')
ylabel('Similarity Trhreshold')

%% Classifier Performance Experiment
aa =14
betaParm1= aa;
betaParm2=[aa-1 aa-2 aa-3 aa-4 aa-5 aa-6];

perf_out = [];

parfor ii=1:length(betaParm2)
    perf_row =struct();
    for jj =1:100
    % Replace the space whale component
    [spaceWhale] =   createRandomSpaceWhale(1, 4, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth,...
        [array_struct.master, array_struct.slave(child_idx)]);
    
    % Copy for each of the methods
    simStructNew = simStruct;
    simStructNew.spaceWhale = spaceWhale;
    
    simStructNew.arrivalTable = UpdateArrTable(simStructNew);
    simStructNew.arrivalArray= UpdateArrArray(simStructNew);
    simStructNew.TDOA_vals = UpdateTDOA(simStructNew);
    simStructNew.betaParm1 = betaParm1;
    simStructNew.betaParm2 = betaParm2(ii);
    simStructNew.cutoff =.95;
    simStructNew.maxEltTime = 15;
    
    
    %     % Update the arrival array and simulation matrix
    %     simStructNew.spaceWhale=spaceWhale;
    %     simStructNew.arrivalTable = UpdateArrTable(simStructNew);
    %     simStructNew.arrivalArray= UpdateArrArray(simStructNew);
    %     simStructNew.TDOA_vals = UpdateTDOA(simStructNew);
    %
    % Max of mean
    simStructNew.Sim_mat = simMatMaxofProd(simStructNew);
    
    simStructNew.chains =updateChainsEncounterFirst(simStructNew);
    simStructNew.Cluster_id= updateClusterID(simStructNew);
    perf = estClassifierPerf(simStructNew);
    
    perf_row(jj).perf = perf;
    
    
    end
    perf_out=[perf_out;perf_row];
    
end
figure



for jj=1:length(betaParm2)
    perf =[];
    for ii =1:50
        
        [prct_improvement,~, ~]= extractClassiferMetrics(perf_out(jj,ii));
        perf(ii)=prct_improvement;
        
    end
    
    subplot(length(betaParm2),1,jj)
    hist(perf);
    title(['Prop Runs Improvement ', num2str(num2str(sum(perf>0)/length(perf))),...
        ' Prop Worse ',  num2str(num2str(sum(perf<0)/length(perf)))])
    ylabel('Simulation Runs')
    xlabel(['Change in Error Rate: Initial Error ', num2str(round(...
        betacdf(.5, betaParm1, betaParm2(jj)),2))])
end



%% Clock Drift Experiment
simStruct.c = 1500;
simStruct.drift=1.5;
simStruct.assSec=2;
simStruct.hydrophone_struct = hydrophone_struct;
simStruct.arrivalArray= UpdateArrArray(simStruct);
simStruct.betaParm1 = betaParm1(1);
simStruct.betaParm2=betaParm2(1);
perf = estClassifierPerf(simStruct);

