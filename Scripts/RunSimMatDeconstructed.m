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
ssp = 1500;
grid_depth = 5;

child_idx = [1:8];
%% Create the agents



[spaceWhale] =   createRandomSpaceWhale(0.5, 6, hyd_arr,...
    array_struct,hydrophone_struct, ssp, grid_depth,...
    [parent, array_struct.slave(child_idx)]);

simStruct=struct();
simStruct.spaceWhale=gather(spaceWhale);
simStruct.array_struct=gather(array_struct);
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
simStruct.arrivalTable = gather(UpdateArrTable(simStruct));
simStruct.arrivalArray= gather(UpdateArrArray(simStruct));
simStruct.TDOA_vals = gather(UpdateTDOA(simStruct));




%% Sensitivity Experiment

%load('ExperimentCallDensityElipseFit.mat')

% Time threshold needs to be ordered maximum to minimum in order for the
% agent sensnsitivity to work. This is because the simulation matrix where
% each row represents an acoustic encounter is only calculated once. In the
% sensisitivy loop, acoustic encounters are extending beyond the maximum
% elapsed time are trimmed such that it's ok to calculate the sensitivy
% threshold when the encounter is too large, but will break when it's too
% small
nIters = 200;
TimeThresh=fliplr(linspace(5, 120, 30));
SimThresh = linspace(.01,.99,20);


ExpScoresMeth_out2D = zeros(length(TimeThresh),...
    length(SimThresh), nIters, 'gpuArray')/0;
ExpScoresMeth_outTDOA = ExpScoresMeth_out2D;
ExpScoresMeth_outMaxProd = ExpScoresMeth_out2D;
ExpScoresMeth_outBaseline = zeros(length(TimeThresh), 1, nIters, 'gpuArray')/0;
allSimStructs =[];

%%
% Number of agents (experiemnts)
nAgents = [3,6,9];

% Structure for output
%AgentExp = struct();
tic

for jj =1:length(nAgents)

    %ExpScoresMeth_out2D = zeros(length(TimeThresh), length(SimThresh), nIters, 'gpuArray')/0;
    ExpScoresMeth_outTDOA = ExpScoresMeth_out2D;
   % ExpScoresMeth_outMaxProd = ExpScoresMeth_out2D;
    nAgent = nAgents(jj);
    parfor ii=1:nIters
        % Replace the space whale component
        [spaceWhale] =   createRandomSpaceWhale(0.5, nAgent, hyd_arr,...
            array_struct,hydrophone_struct, ssp, grid_depth,...
            [array_struct.master, array_struct.slave(child_idx)]);
        
        % Copy for each of the methods
        simStructNew = simStruct;
        simStructNew.maxEltTime = max(TimeThresh);
        
        % Update the arrival array and simulation matrix
        simStructNew.spaceWhale=spaceWhale;
        simStructNew.arrivalTable = UpdateArrTable(simStructNew);
        simStructNew.arrivalArray= UpdateArrArray(simStructNew);
        simStructNew.TDOA_vals = UpdateTDOA(simStructNew);
        
        % Baseline model
        simStructBaseline = simStructNew;
        [ExpScoresMethBaseline, ~] = runSensitivtyLp(simStructBaseline,TimeThresh);
         ExpScoresMeth_outBaseline(:,:,ii) = ExpScoresMethBaseline;
        
        
        % Create copy for TDOA method
        simStructTDOA = simStructNew;
        
        % TDOA only method
        simStructTDOA.Sim_mat = simMatTDOAonly(simStructTDOA);
        
        %Run the sensitivity loop
        [ExpScoresMethTDOA, ~] = runSensitivtyLp(simStructTDOA,TimeThresh,SimThresh);
        ExpScoresMeth_outTDOA(:,:,ii) = ExpScoresMethTDOA;

%         % Baseline model
%         simStructBaseline = simStructNew;
%         [ExpScoresMethBaseline, ~] = runSensitivtyLp(simStructBaseline,TimeThresh);
%         ExpScoresMeth_outBaseline(:,:,ii) = ExpScoresMethBaseline;
        
        
        % Create copy for TDOA method
        simStructTDOA = simStructNew;
        
        % TDOA only method
        simStructTDOA.Sim_mat = simMatTDOAonly(simStructTDOA);
        
        
        %Run the sensitivity loop
        [ExpScoresMethTDOA, ~] = runSensitivtyLp(simStructTDOA,TimeThresh,SimThresh);
        ExpScoresMeth_outTDOA(:,:,ii) = ExpScoresMethTDOA;
        
%         % Max of prod
%         simStructMaxProd = simStructNew;
%         
%   
%         % Product Method
%         simStructMaxProd.Sim_mat = simMatMaxofProd(simStructMaxProd);
%    
%         
%         %Run the sensitivity loop
%         [ExpScoresMethMaxMean, ~] = runSensitivtyLp(simStructMaxProd,TimeThresh,SimThresh);
%         ExpScoresMeth_outMaxProd(:,:,ii) = ExpScoresMethMaxMean;
        
        ii
        
    end
      AgentExp(jj).TDOA =ExpScoresMeth_outTDOA;
%      AgentExp(jj).MaxProd = ExpScoresMeth_outMaxProd;
%      AgentExp(jj).Baseline = ExpScoresMeth_outBaseline;
    
    
    
end
toc
%%
save('ExperimentCallDensityElipseFit.mat', AgentExp)

%% Plot experiment 1
% 3 Agents

%output table, experiment
base_agents =gather(nanmedian(gather(squeeze(AgentExp(1).Baseline),2)));
tdoa_agents = gather(nanmedian(gather(squeeze(AgentExp(1).TDOA),2)));
spatial_agents = gather(nanmedian(gather(squeeze(AgentExp(1).MaxProd),2)));

quantile(base_agents, .5)
quantile(tdoa_agents(:), .5)
quantile(spatial_agents(:), .5)

iqr(base_agents)
iqr(tdoa_agents(:))
iqr(spatial_agents(:))


vals = cat(3, nanmedian(squeeze(AgentExp(1).Baseline),2),...
    nanmedian(squeeze(AgentExp(1).TDOA),2),...
    nanmedian(squeeze(AgentExp(1).MaxProd),2));
figure(1)
subplot(2,2,1)
plot(TimeThresh, nanmedian(squeeze(AgentExp(1).Baseline),2))
title('Baseline')
xlabel('Time Threshold (s)')
ylabel('Adjusted Rand Index')


subplot(2,2,2)
imagesc(TimeThresh,SimThresh, nanmedian(AgentExp(1).TDOA,3)),  axis xy , colorbar
caxis([-.1 .9])
caxis(gather([min(vals(:)) max(vals(:))]))
title('TDOA only')
xlabel('Time Threshold (s)')
ylabel('Similarity Trhreshold')
c = colorbar;
c.Label.String = 'Adjusted Rand Index';

subplot(2,2,3)
imagesc(TimeThresh, SimThresh, nanmedian(AgentExp(1).MaxProd,3)),  axis xy,  colorbar
caxis([-.1 .9])
caxis(gather([min(vals(:)) max(vals(:))]))
title('Spatial Method')
xlabel('Time Threshold (s)')
ylabel('Similarity Trhreshold')
c = colorbar;
c.Label.String = 'Adjusted Rand Index';



% 6 Agents
base_agents =gather(nanmedian(gather(squeeze(AgentExp(2).Baseline),2)));
tdoa_agents = gather(nanmedian(gather(squeeze(AgentExp(2).TDOA),2)));
spatial_agents = gather(nanmedian(gather(squeeze(AgentExp(2).MaxProd),2)));

quantile(base_agents, .5)
quantile(tdoa_agents(:), .5)
quantile(spatial_agents(:), .5)

iqr(base_agents)
iqr(tdoa_agents(:))
iqr(spatial_agents(:))
vals = cat(3, nanmedian(squeeze(AgentExp(2).Baseline),2),...
    nanmedian(squeeze(AgentExp(2).TDOA),2),...
    nanmedian(squeeze(AgentExp(2).MaxProd),2));

figure(2)
subplot(2,2,1)
plot(TimeThresh, nanmedian(squeeze(AgentExp(2).Baseline),2))
title('Baseline')
xlabel('Time Threshold (s)')
ylabel('Adjusted Rand Index')

subplot(2,2,2)
imagesc(TimeThresh,SimThresh, nanmedian(AgentExp(2).TDOA,3)),  axis xy , colorbar
caxis([-.1 .6])
caxis(gather([min(vals(:)) max(vals(:))]))
title('TDOA only')
xlabel('Time Threshold (s)')
ylabel('Similarity Trhreshold')
c = colorbar;
c.Label.String = 'Adjusted Rand Index';


subplot(2,2,3)
imagesc(TimeThresh,SimThresh,  nanmedian(AgentExp(2).MaxProd,3)),  axis xy,  colorbar
caxis(gather([min(vals(:)) max(vals(:))]))
title('Spatial Method')
xlabel('Time Threshold (s)')
ylabel('Similarity Trhreshold')
c = colorbar;
c.Label.String = 'Adjusted Rand Index';


% 9 Agents

base_agents =gather(nanmedian(gather(squeeze(AgentExp(3).Baseline),2)));
tdoa_agents = gather(nanmedian(gather(squeeze(AgentExp(3).TDOA),2)));
spatial_agents = gather(nanmedian(gather(squeeze(AgentExp(3).MaxProd),2)));

quantile(base_agents, .5)
quantile(tdoa_agents(:), .5)
quantile(spatial_agents(:), .5)

iqr(base_agents)
iqr(tdoa_agents(:))
iqr(spatial_agents(:))
vals = cat(3, nanmedian(squeeze(AgentExp(3).Baseline),2),...
    nanmedian(squeeze(AgentExp(3).TDOA),2),...
    nanmedian(squeeze(AgentExp(3).MaxProd),2));

figure(3)
subplot(2,2,1)
plot(TimeThresh, nanmedian(squeeze(AgentExp(3).Baseline),2))
title('Baseline')
xlabel('Time Threshold (s)')
c.Label.String = 'Adjusted Rand Index';


subplot(2,2,2)
imagesc(TimeThresh,SimThresh, nanmedian(AgentExp(3).TDOA,3)),  axis xy , colorbar
caxis(gather([min(vals(:)) max(vals(:))]))
title('TDOA only')
xlabel('Time Threshold (s)')
ylabel('Similarity Trhreshold')
c = colorbar;
c.Label.String = 'Adjusted Rand Index';


subplot(2,2,3)
imagesc(TimeThresh,SimThresh,  nanmedian(AgentExp(3).MaxProd,3)),  axis xy,  colorbar
caxis(gather([min(vals(:)) max(vals(:))]))
title('Spatial Method')
xlabel('Time Threshold (s)')
ylabel('Similarity Trhreshold')
c = colorbar;
c.Label.String = 'Adjusted Rand Index';


%% Classifier Performance Experiment


%load('ExperimentClassifyPerfElipseFit.mat')
nIters =50;
TimeThresh=fliplr(linspace(5, 120, 30));
SimThresh = linspace(.01,.99,20);

aa =14;
betaParm1= 7;
betaParm2=[2.5, 2, 1.5];

perf_out_baseline1 = zeros(length(TimeThresh), length(SimThresh), nIters);
perf_out_TDOA1 =perf_out_baseline1;
perf_out_MaxProd1 = perf_out_baseline1;

perf_out_baseline2 = perf_out_baseline1;
perf_out_TDOA2 =perf_out_baseline1;
perf_out_MaxProd2 = perf_out_baseline1;


perf_out_baseline3 = perf_out_baseline1;
perf_out_TDOA3 =perf_out_baseline1;
perf_out_MaxProd3 = perf_out_baseline1;


ClassifierPerfExp =struct();
% Create filter grids to knock out ambiguity surfaces 
simStruct.filtGrid = createDetRangeFiltGrid(simStruct, hydrophone_struct);

 for jj =1:nIters
        % Replace the space whale component
        [spaceWhale] =   createRandomSpaceWhale(0.5, 6, hyd_arr,...
            array_struct,hydrophone_struct, ssp, grid_depth,...
            [array_struct.master, array_struct.slave(child_idx)]);
        
        % Copy for each of the methods
        simStructNew = simStruct;
        simStructNew.spaceWhale = spaceWhale;
        simStructNew.arrivalTable = gather(UpdateArrTable(simStructNew));
        simStructNew.arrivalArray= gather(UpdateArrArray(simStructNew));
        simStructNew.TDOA_vals = gather(UpdateTDOA(simStructNew));
        simStructNew.maxEltTime = gather(max(TimeThresh));
        
        % Baseline
        simStructBaseline =simStructNew;
        
        % TDOA
        simStructTDOA =simStructNew;
        simStructTDOA.Sim_mat = simMatTDOAonly(simStructTDOA);
        
        % Max of mean
        simStructMaxProd =simStructNew;
        simStructMaxProd.Sim_mat = simMatMaxofProd(simStructMaxProd);
        
        b2 =betaParm2(1);
        [aa, bb, cc] = classifierLpFx(simStructMaxProd,simStructTDOA,simStructBaseline,...
            TimeThresh, SimThresh, betaParm1, b2);
        perf_out_baseline1(:,jj) = aa;
        perf_out_TDOA1(:,:,jj) =bb;
        perf_out_MaxProd1(:,:,jj) = cc;
        
        
        
        b2 =betaParm2(2);
        [aa, bb, cc] = classifierLpFx(simStructNew,simStructTDOA,simStructBaseline,...
            TimeThresh, SimThresh, betaParm1, b2);
        perf_out_baseline2(:,jj) = aa;
        perf_out_TDOA2(:,:,jj) =bb;
        perf_out_MaxProd2(:,:,jj) = cc;
        
        
        simStructNew = simStructTDOA;
        b2 =betaParm2(3);
        [aa, bb, cc] = classifierLpFx(simStructNew,simStructTDOA,simStructBaseline,...
            TimeThresh, SimThresh, betaParm1, b2);
        perf_out_baseline3(:,jj) = aa;
        perf_out_TDOA3(:,:,jj) =bb;
        perf_out_MaxProd3(:,:,jj) = cc;
        
        
        
        

 end
    
 %%
   ClassifierPerfExp(1).MaxProd =gather(perf_out_MaxProd1);
   ClassifierPerfExp(2).MaxProd =gather(perf_out_MaxProd2);
   ClassifierPerfExp(3).MaxProd =gather(perf_out_MaxProd3);
   
   ClassifierPerfExp(1).TDOA =gather(perf_out_TDOA1);
   ClassifierPerfExp(2).TDOA =gather(perf_out_TDOA2);
   ClassifierPerfExp(3).TDOA =gather(perf_out_TDOA3);
   
   
   ClassifierPerfExp(1).Baseline =gather(perf_out_baseline1);
   ClassifierPerfExp(2).Baseline =gather(perf_out_baseline2);
   ClassifierPerfExp(3).Baseline =gather(perf_out_baseline3);
    



%save('ExperimentClassifyPerfAssocElipseFit_prctileET_new.mat', 'ClassifierPerfExp')


%load('ExperimentClassifyPerfAssocElipseFit_prctileET_new.mat')
aa = quantile(ClassifierPerfExp(1).MaxProd,.5,3);
bb = quantile(ClassifierPerfExp(2).MaxProd,.5,3);
cc = quantile(ClassifierPerfExp(3).MaxProd,.5,3);


figure;
subplot(2,2,1)
contourf(SimThresh,TimeThresh, aa,5),axis xy , colorbar
xlabel('Time Threshold (s)'); ylabel('Similarity Threshold');
caxis([-.6  .3])
title('Spatial Method- Low Q Classifier')
subplot(2,2,2)
contourf(SimThresh,TimeThresh,bb, 5),axis xy , colorbar
xlabel('Time Threshold (s)'); ylabel('Similarity Threshold');
caxis([-.6  .3]);
title('Spatial Method- Mod Q Classifier')
subplot(2,2,3)
contourf(SimThresh,TimeThresh,cc,5),axis xy , colorbar
xlabel('Time Threshold (s)'); ylabel('Similarity Threshold');
caxis([-.6  .3]);
title('Spatial Method- High Q Classifier')


aa = quantile(ClassifierPerfExp(1).TDOA,.5,3);
bb = quantile(ClassifierPerfExp(2).TDOA,.5,3);
cc = quantile(ClassifierPerfExp(3).TDOA,.5,3);

figure;
subplot(2,2,1)
contourf(SimThresh,TimeThresh, aa,5),axis xy , colorbar
xlabel('Time Threshold (s)'); ylabel('Similarity Threshold');
caxis([-.6  .3])
title('TDOA Method- Low Q Classifier')
subplot(2,2,2)
contourf(SimThresh,TimeThresh,bb,5),axis xy , colorbar
xlabel('Time Threshold (s)'); ylabel('Similarity Threshold');
caxis([-.6  .3])
title('TDOA Method- Mod Q Classifier')
subplot(2,2,3)
contourf(SimThresh,TimeThresh,cc,5),axis xy , colorbar
caxis([-.6  .3])
xlabel('Time Threshold (s)'); ylabel('Similarity Threshold');
title('TDOA Method- High Q Classifier')

aa = quantile(ClassifierPerfExp(1).Baseline(:,1:nIters),[ .5 .25 .75],2);
bb = quantile(ClassifierPerfExp(2).Baseline(:,1:nIters),[.5 .25 .75],2);
cc = quantile(ClassifierPerfExp(3).Baseline(:,1:nIters),[.5 .25 .75],2);

figure;
subplot(2,2,1)
plot_ci(TimeThresh,aa, 'PatchColor', 'r', 'PatchAlpha', 0.2, ...
    'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', 'b', ...
    'LineWidth', 1.5, 'LineStyle','--', 'LineColor', 'k')
xlabel('Time Threshold (s)'); ylabel('Change in Accuracy')
title('Baseline Method- Low Q Classifier')
subplot(2,2,2)
plot_ci(TimeThresh,bb, 'PatchColor', 'r', 'PatchAlpha', 0.2, ...
    'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', 'b', ...
    'LineWidth', 1.5, 'LineStyle','--', 'LineColor', 'k')
title('Baseline Method- Mod Q Classifier')
xlabel('Time Threshold (s)'); ylabel('Change in Accuracy')
subplot(2,2,3)

plot_ci(TimeThresh,cc, 'PatchColor', 'r', 'PatchAlpha', 0.2, ...
    'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', 'b', ...
    'LineWidth', 1.5, 'LineStyle','--', 'LineColor', 'k')
xlabel('Time Threshold (s)'); ylabel('Change in Accuracy')
title('Baseline Method- High Q Classifier')




%%)
%%
bins = linspace(-1, 1, 40);
for jj=1:length(betaParm2)
    perfMaxProd =[];
    perfTDOA =[];
    perfBaseline=[];
    for ii =1:nIters
        
        [prct_improvement,~, ~]= extractClassiferMetrics(ClassifierPerfExp(jj).Baseline(ii));
        perfBaseline(ii)=prct_improvement;
        
        [prct_improvement,~, ~]= extractClassiferMetrics(ClassifierPerfExp(jj).TDOA(ii));
        perfTDOA(ii)=prct_improvement;

        [prct_improvement,~, ~]= extractClassiferMetrics(ClassifierPerfExp(jj).MaxProd(ii));
        perfMaxProd(ii)=prct_improvement;
        
        
        
    end
    
    median_error =betacdf(.5, betaParm1, betaParm2(jj));
    
    violinData = flipud([perfBaseline'; perfTDOA'; perfMaxProd']);
    violinCell = flipud([repmat({'Baseline'}, size(perfBaseline'));...
        repmat({'TDOA'}, size(perfTDOA'));...
        repmat({'Spatial'}, size(perfMaxProd'))]);
    
    violinData(isinf(violinData)) =0;
    
    % Proportion of Runs Improved
     IMP(1)= sum(perfBaseline>0)/length(perfBaseline)*100;
     IMP(2)= sum(perfTDOA>0)/length(perfTDOA)*100;
     IMP(3)= sum(perfMaxProd>0)/length(perfMaxProd)*100;
     
     IMP = round(IMP);
     txt{1} = {[num2str(IMP(1)) '% Runs']; 'Improved Accuracy'};
     txt{2} = {[num2str(IMP(2)) '% Runs']; 'Improved Accuracy'}
     txt{3} = {[num2str(IMP(3)) '% Runs']; 'Improved Accuracy'}
    
    texty = max(violinData)+.5;
    figure(9)
    subplot(3,1,jj)
    violinplot(violinData, violinCell, 'ViolinColor', [.5 .5 .5]);
    title(['Median Classification Error: ' num2str(median_error)])
    text(1:3,[texty,texty,texty], txt,'HorizontalAlignment','center')
    ylim([-2 texty+.7])
    
    
    
%     
%     figure(4)
%     subplot(length(betaParm2),1,jj)
%     %perfBaseline = perfBaseline(perfBaseline~=0);
%     hist(perfBaseline, bins);
%     title(['Prop Runs Improvement Baseline', num2str(num2str(sum(perfBaseline>0)/length(perfBaseline))),...
%         ' Prop Worse ',  num2str(num2str(sum(perfBaseline<0)/length(perfBaseline)))])
%     ylabel('Simulation Runs')
%     xlabel(['Change in Error Rate: Initial Error ', num2str(round(...
%         betacdf(.5, betaParm1, betaParm2(jj)),2))])
%     xlim([-.8 .8])
%     
%     
%     figure(5)
%     subplot(length(betaParm2),1,jj)
%     %perfTDOA = perfTDOA(perfTDOA~=0);
%     hist(perfTDOA, bins);
%     title(['Prop Runs Improvement TDOA', num2str(num2str(sum(perfTDOA>0)/length(perfTDOA))),...
%         ' Prop Worse ',  num2str(num2str(sum(perfTDOA<0)/length(perfTDOA)))])
%     ylabel('Simulation Runs')
%     xlabel(['Change in Error Rate: Initial Error ', num2str(round(...
%         betacdf(.5, betaParm1, betaParm2(jj)),2))])
%     xlim([-.8 .8])
%     
%     
%     
%     figure(6)
%     subplot(length(betaParm2),1,jj)
%     %perfMaxProd = perfMaxProd(perfMaxProd~=0);
%     hist(perfMaxProd, bins);
%     title(['Prop Runs Improvement MaxProd', num2str(num2str(sum(perfMaxProd>0)/length(perfMaxProd))),...
%         ' Prop Worse ',  num2str(num2str(sum(perfMaxProd<0)/length(perfMaxProd)))])
%     ylabel('Simulation Runs')
%     xlabel(['Change in Error Rate: Initial Error ', num2str(round(...
%         betacdf(.5, betaParm1, betaParm2(jj)),2))])
%     xlim([-.8 .8])
%     
%     (sum(perfMaxProd>0)/length(perfMaxProd))-(sum(perfMaxProd<0)/length(perfMaxProd))
%         
    
end



%% Clock Drift Experiment
load('ExperimentClassifyPerfAssocElipseFit.mat')

simStruct.c = 1500;
simStruct.betaParm1 = betaParm1(1);
simStruct.betaParm2=betaParm2(1);
simStruct.cutoff =.85;
simStruct.maxEltTime = 30;
simStruct.hydrophone_struct = hydrophone_struct;
clock_drift = [1 2 3 4];
association_sec = [0 0 4 6];


perf_assocExp =struct();
nIters = 200;

perf_out_baseline = struct();
perf_out_TDOA = struct();
perf_out_MaxProd = struct();


for ii=1:length(association_sec)
    
    assoc = association_sec(ii);
    cdrift = clock_drift(ii);
    
    parfor kk = 1:nIters
        % Replace the space whale component
        [spaceWhale] =   createRandomSpaceWhale(0.5, 6, hyd_arr,...
            array_struct,hydrophone_struct, ssp, grid_depth,...
            [array_struct.master, array_struct.slave(child_idx)]);
        
        simStructNew = simStruct;
        simStructNew.spaceWhale = spaceWhale;
        simStructNew.drift=cdrift;
        simStructNew.assSec=assoc;
        simStructNew.arrivalTable = UpdateArrTable(simStructNew);
        simStructNew.arrivalArray= UpdateArrArray(simStructNew);
        simStructNew.TDOA_vals = UpdateTDOA(simStructNew);
        
        
        
        % Baseline
        simStructBaseline =simStructNew;
        simStructBaseline.Cluster_id = acEnc(simStructBaseline);
        perf_out_baseline(kk).Perf =estClassifierPerf(simStructBaseline);    
        
        

        % TDOA
        simStructTDOA =simStructNew;
        simStructTDOA.Sim_mat = simMatTDOAonly(simStructTDOA);
        simStructTDOA.chains =updateChainsEncounterFirst(simStructTDOA);
        simStructTDOA.Cluster_id= updateClusterID(simStructTDOA);
        perf_out_TDOA(kk).Perf= estClassifierPerf(simStructTDOA);
        
        % Max of mean
        simStructNew.Sim_mat = simMatMaxofProd(simStructNew);
        simStructNew.chains =updateChainsEncounterFirst(simStructNew);
        simStructNew.Cluster_id= updateClusterID(simStructNew);
        perf_out_MaxProd(kk).Perf= estClassifierPerf(simStructNew);
        
        
        
    end
    perf_assocExp(ii).Baseline = perf_out_baseline;
    perf_assocExp(ii).TDOA = perf_out_TDOA;
    perf_assocExp(ii).MaxProd = perf_out_MaxProd;
end


save('ExperimentClassifyPerfAssocElipseFit.mat', perf_assocExp)

%%

bins = linspace(-.8, .9, 40);
for jj=1:length(clock_drift)
    perfMaxProd =[];
    perfTDOA =[];
    perfBaseline=[];
    
    for ii =1:nIters
        
        [prct_improvement,~, ~]= extractClassiferMetrics(perf_assocExp(jj).Baseline(ii));
        perfBaseline(ii)=prct_improvement;
        
        [prct_improvement,~, ~]= extractClassiferMetrics(perf_assocExp(jj).TDOA(ii));
        perfTDOA(ii)=prct_improvement;

        [prct_improvement,~, ~]= extractClassiferMetrics(perf_assocExp(jj).MaxProd(ii));
        perfMaxProd(ii)=prct_improvement;
        
        
        
    end
    
    figure(7)
    subplot(length(clock_drift),1,jj)
    %perfBaseline = perfBaseline(perfBaseline~=0);
    hist(perfBaseline, bins);
    title(['Prop Runs Improvement Baseline', num2str(num2str(sum(perfBaseline>0)/length(perfBaseline))),...
        ' Prop Worse ',  num2str(num2str(sum(perfBaseline<0)/length(perfBaseline)))])
    ylabel('Simulation Runs')
    xlabel(['Change in Error Rate: Initial Error '])
    xlim([-.8 .8])
    
    
    figure(8)
    subplot(length(clock_drift),1,jj)
    %perfTDOA = perfTDOA(perfTDOA~=0);
    hist(perfTDOA, bins);
    title(['Prop Runs Improvement TDOA', num2str(num2str(sum(perfTDOA>0)/length(perfTDOA))),...
        ' Prop Worse ',  num2str(num2str(sum(perfTDOA<0)/length(perfTDOA)))])
    ylabel('Simulation Runs')
    xlabel(['Change in Error Rate: Initial Error '])
    xlim([-.8 .8])
    
    
    
    figure(9)
    subplot(length(clock_drift),1,jj)
    %perfMaxProd = perfMaxProd(perfMaxProd~=0);
    hist(perfMaxProd, bins);
    title(['Prop Runs Improvement MaxProd', num2str(num2str(sum(perfMaxProd>0)/length(perfMaxProd))),...
        ' Prop Worse ',  num2str(num2str(sum(perfMaxProd<0)/length(perfMaxProd)))])
    ylabel('Simulation Runs')
    xlabel(['Change in Error Rate: Initial Error '])
    xlim([-.8 .8])
    
        
    
end













