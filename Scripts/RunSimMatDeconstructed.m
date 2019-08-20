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
simStruct.s = 8;
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

% Time threshold needs to be ordered maximum to minimum in order for the
% agent sensnsitivity to work. This is because the simulation matrix where
% each row represents an acoustic encounter is only calculated once. In the
% sensisitivy loop, acoustic encounters are extending beyond the maximum
% elapsed time are trimmed such that it's ok to calculate the sensitivy
% threshold when the encounter is too large, but will break when it's too
% small
nIters =200;
TimeThresh=fliplr(linspace(1, 120, 50));
SimThresh = linspace(0,1,30);
SimThresh1 =linspace(0,1,30);

ExpScoresMeth_out2D = zeros(length(TimeThresh), length(SimThresh), nIters)/0;
ExpScoresMeth_outTDOA = ExpScoresMeth_out2D;
ExpScoresMeth_outMaxProd = ExpScoresMeth_out2D;

allSimStructs =[];

% Number of agents (experiemnts)
nAgents = [3,6,9];

% Structure for output
AgentExp = struct();


tic

for jj =1:length(nAgents)
    nAgent = nAgents(jj);
    parfor ii=1:nIters
        % Replace the space whale component
        [spaceWhale] =   createRandomSpaceWhale(1, nAgent, hyd_arr,...
            array_struct,hydrophone_struct, ssp, grid_depth,...
            [array_struct.master, array_struct.slave(child_idx)]);
        
        % Copy for each of the methods
        simStructNew = simStruct;
        
        
        % Update the arrival array and simulation matrix
        simStructNew.spaceWhale=spaceWhale;
        simStructNew.arrivalTable = UpdateArrTable(simStructNew);
        simStructNew.arrivalArray= UpdateArrArray(simStructNew);
        simStructNew.TDOA_vals = UpdateTDOA(simStructNew);
        
        
        % Create copy for TDOA method
        simStructTDOA = simStructNew;
        
        % TDOA only method
        simStructTDOA.Sim_mat = simMatTDOAonly(simStructTDOA);
        
        %Run the sensitivity loop
        [ExpScoresMethTDOA, ~] = runSensitivtyLp(simStructTDOA,TimeThresh,SimThresh);
        ExpScoresMeth_outTDOA(:,:,ii) = ExpScoresMethTDOA;
        %
        % Max of mean
        simStructMaxMean = simStructNew;
        
        
        % Product Method
        simStructMaxMean.Sim_mat = simMatMaxofProd(simStructMaxMean);
        
        
        %Run the sensitivity loop
        [ExpScoresMethMaxMean, ~] = runSensitivtyLp(simStructMaxMean,TimeThresh,SimThresh);
        ExpScoresMeth_outMaxProd(:,:,ii) = ExpScoresMethMaxMean;
        
        ii
        
    end
    AgentExp(jj).TDOA =ExpScoresMeth_outTDOA;
    AgentExp(jj).MaxProd = ExpScoresMeth_outMaxProd;
    
    
    
end
toc


figure(1)
% subplot(3,1,1)
% imagesc(SimThresh1nanmedian(ExpScoresMeth_out2D,3)) ,  axis xy, colorbar
% title('Median ARI Max of Mean')
%
%
subplot(2,1,1)
imagesc(TimeThresh,SimThresh, nanmedian(ExpScoresMeth_outTDOA,3)) ,  axis xy  colorbar
title('TDOA only')



subplot(2,1,2)
imagesc(TimeThresh,SimThresh1, nanmean(ExpScoresMeth_outMaxProd,3)) ,  axis xy  colorbar
title('Spatial Method')
xlabel('Time Threshold (s)')
ylabel('Similarity Trhreshold')

%% Classifier Performance Experiment
aa =14;
betaParm1= aa;
betaParm2=[aa-2 aa-5 aa-7];

perf_out = [];
nruns =500;

timethresh = 10:10:60;

parfor ii=1:length(betaParm2)
    perf_row =struct();
    beta2 = betaParm2(ii);
    
    for jj =1:nruns
        % Replace the space whale component
        [spaceWhale] =   createRandomSpaceWhale(1, 6, hyd_arr,...
            array_struct,hydrophone_struct, ssp, grid_depth,...
            [array_struct.master, array_struct.slave(child_idx)]);
        
        % Copy for each of the methods
        simStructNew = simStruct;
        simStructNew.spaceWhale = spaceWhale;
        simStructNew.arrivalTable = UpdateArrTable(simStructNew);
        simStructNew.arrivalArray= UpdateArrArray(simStructNew);
        simStructNew.TDOA_vals = UpdateTDOA(simStructNew);
        simStructNew.betaParm1 = betaParm1;
        simStructNew.betaParm2 = beta2;
        simStructNew.cutoff =.85;
        simStructNew.maxEltTime = 20;
        
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

bins = linspace(-.8, .6, 20);

for jj=1:length(betaParm2)
    perf =[];
    for ii =1:nruns
        
        [prct_improvement,~, ~]= extractClassiferMetrics(perf_out(jj,ii));
        perf(ii)=prct_improvement;
        
    end
    
    subplot(length(betaParm2),1,jj)
    perf = perf(perf~=0);
    hist(perf, bins);
    title(['Prop Runs Improvement ', num2str(num2str(sum(perf>0)/length(perf))),...
        ' Prop Worse ',  num2str(num2str(sum(perf<0)/length(perf)))])
    ylabel('Simulation Runs')
    xlabel(['Change in Error Rate: Initial Error ', num2str(round(...
        betacdf(.5, betaParm1, betaParm2(jj)),2))])
    xlim([-.8 .8])
    
    (sum(perf>0)/length(perf))-(sum(perf<0)/length(perf))
    
end



%% Clock Drift Experiment

simStruct.c = 1500;
simStruct.betaParm1 = betaParm1(1);
simStruct.betaParm2=betaParm2(1);
simStruct.cutoff =.85;
simStruct.maxEltTime = 20;
simStruct.hydrophone_struct = hydrophone_struct;
clock_drift = [1 2 3 4];
association_sec = [0 4 6];


perf_assocExp =struct();
nRuns = 10;
for ii=1:length(association_sec)
    
    assoc = association_sec(ii);
    perf_driftExp =[];
    for jj = 1:length(clock_drift)
        cdrift = clock_drift(jj);
        perf_row =struct();
        for kk = 1:nRuns
            % Replace the space whale component
            [spaceWhale] =   createRandomSpaceWhale(1, 6, hyd_arr,...
                array_struct,hydrophone_struct, ssp, grid_depth,...
                [array_struct.master, array_struct.slave(child_idx)]);
            
            % Copy for each of the methods
            simStructNew = simStruct;
            simStructNew.spaceWhale = spaceWhale;
            simStructNew.drift=cdrift;
            simStructNew.assSec=assoc;
            
            simStructNew.arrivalTable = UpdateArrTable(simStructNew);
            simStructNew.arrivalArray= UpdateArrArray(simStructNew);
            simStructNew.TDOA_vals = UpdateTDOA(simStructNew);
            
            
            % Max of mean
            simStructNew.Sim_mat = simMatMaxofProd(simStructNew);
            simStructNew.chains =updateChainsEncounterFirst(simStructNew);
            simStructNew.Cluster_id= updateClusterID(simStructNew);
            perf = estClassifierPerf(simStructNew);
            
            perf_row(kk).perf = perf;
            
        end
        perf_driftExp=[perf_driftExp; perf_row];
        
    end
    
    perf_assocExp(ii).AssocExp = perf_driftExp;
end

















