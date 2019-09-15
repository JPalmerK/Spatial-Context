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


% Create filter grids to knock out ambiguity surfaces
simStruct.filtGrid = createDetRangeFiltGrid(simStruct, hydrophone_struct);


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
TimeThresh=fliplr(linspace(2, 100, 30));
SimThresh = linspace(.01,.99,20);


ExpScoresMeth_out2D = zeros(length(TimeThresh),...
    length(SimThresh), nIters)/0;
ExpScoresMeth_outTDOA = ExpScoresMeth_out2D;
ExpScoresMeth_outMaxProd = ExpScoresMeth_out2D;
ExpScoresMeth_outBaseline = zeros(length(TimeThresh), 1, nIters)/0;
allSimStructs =[];

%%
% Number of agents (experiemnts)
nAgents = [3,6,9];

nIters = 4;
% Structure for output
AgentExp = struct();
tic

for jj =1:length(nAgents)

    ExpScoresMeth_out2D = zeros(length(TimeThresh), length(SimThresh), nIters, 'gpuArray')/0;
    ExpScoresMeth_outTDOA = ExpScoresMeth_out2D;
    ExpScoresMeth_outMaxProd = ExpScoresMeth_out2D;
    nAgent = nAgents(jj);
    parfor ii=1:nIters
        % Replace the space whale component
        [spaceWhale] =   createRandomSpaceWhale(0.5, nAgent, hyd_arr,...
            array_struct,hydrophone_struct, ssp, grid_depth,...
            [array_struct.master, array_struct.slave(child_idx)]);
        
        % Copy for each of the methods
        simStructNew = simStruct;
        simStructNew.spaceWhale = spaceWhale;
        simStructNew.arrivalTable = gather(UpdateArrTable(simStructNew));
        simStructNew.arrivalArray= gather(UpdateArrArray(simStructNew));
        simStructNew.TDOA_vals = gather(UpdateTDOA(simStructNew));
        simStructNew.maxEltTime = gather(max(TimeThresh));
        simStructNew.truthTable = createTruthTable(simStructNew);
        
       
        
        % Baseline model
        simStructBaseline = simStructNew;
        [ExpScoresMethBaseline, ~] = runSensitivtyLp(simStructBaseline,TimeThresh);
         ExpScoresMeth_outBaseline(:,:,ii) = ExpScoresMethBaseline;

        % TDOA only method
        simStructTDOA = simStructNew;
        simStructTDOA.Sim_mat = simMatTDOAonly(simStructTDOA);
        [ExpScoresMethTDOA, ~] = runSensitivtyLp(simStructTDOA,TimeThresh,SimThresh);
        ExpScoresMeth_outTDOA(:,:,ii) = ExpScoresMethTDOA;
        
        % Max of prod
        simStructMaxProd = simStructNew;
        simStructMaxProd.Sim_mat = simMatMaxofProd(simStructMaxProd);
        [ExpScoresMethMaxProd, ~] = runSensitivtyLp(simStructMaxProd,TimeThresh,SimThresh);
        ExpScoresMeth_outMaxProd(:,:,ii) = ExpScoresMethMaxProd;
        
        ii
        
    end
    AgentExp(jj).TDOA =ExpScoresMeth_outTDOA;
    AgentExp(jj).MaxProd = ExpScoresMeth_outMaxProd;
    AgentExp(jj).Baseline = ExpScoresMeth_outBaseline;
    
    
    
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
nIters =20;
TimeThresh=fliplr(linspace(2, 100, 30));
SimThresh = linspace(.01,.99,20);

aa =14;
betaParm1= 5;
betaParm2_1=2.5;
betaParm2_2=2;
betaParm2_3= 1.5;

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

UnAidedPerf1= zeros(nIters, 1);
UnAidedPerf2= zeros(nIters, 1);
UnAidedPerf3= zeros(nIters, 1);

parfor jj =1:nIters
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
    simStructNew.truthTable = createTruthTable(simStructNew);
    
    % Baseline
    simStructBaseline =simStructNew;
    
    % TDOA
    simStructTDOA =simStructNew;
    simStructTDOA.Sim_mat = simMatTDOAonly(simStructTDOA);
    
    % Max of mean
    simStructMaxProd =simStructNew;
    simStructMaxProd.Sim_mat = simMatMaxofProd(simStructMaxProd);
    
    [aa, bb, cc, UnAidedPerf] = classifierLpFx(simStructMaxProd,...
        simStructTDOA,simStructBaseline,...
        TimeThresh, SimThresh, betaParm1, betaParm2_1);
    perf_out_baseline1(:,jj) = aa;
    perf_out_TDOA1(:,:,jj) =bb;
    perf_out_MaxProd1(:,:,jj) = cc;
    UnAidedPerf1(jj)= unique(UnAidedPerf);
    
    
    [aa, bb, cc, UnAidedPerf] = classifierLpFx(simStructMaxProd,simStructTDOA,simStructBaseline,...
        TimeThresh, SimThresh, betaParm1, betaParm2_2);
    perf_out_baseline2(:,jj) = aa;
    perf_out_TDOA2(:,:,jj) =bb;
    perf_out_MaxProd2(:,:,jj) = cc;
    UnAidedPerf2(jj)=unique(UnAidedPerf);
    
    
    [aa, bb, cc, UnAidedPerf] = classifierLpFx(simStructMaxProd,simStructTDOA,simStructBaseline,...
        TimeThresh, SimThresh, betaParm1, betaParm2_3);
    perf_out_baseline3(:,jj) = aa;
    perf_out_TDOA3(:,:,jj) =bb;
    perf_out_MaxProd3(:,:,jj) = cc;
    UnAidedPerf3(jj)=unique(UnAidedPerf);
    

    
end

save('ExperimentClassifyPerfElipseFit.mat','perf_out_baseline1',...
    'perf_out_baseline2','perf_out_baseline3','perf_out_TDOA1','perf_out_TDOA2',...
    'perf_out_TDOA3', 'perf_out_MaxProd1','perf_out_MaxProd2','perf_out_MaxProd3',...
    'UnAidedPerf1','UnAidedPerf2','UnAidedPerf3')


 %% Plot the performance experiments
 
 %load('ExperimentClassifyPerfElipseFit.mat')
 
 close all;
   ClassifierPerfExp(1).MaxProd =gather(perf_out_MaxProd1);
   ClassifierPerfExp(2).MaxProd =gather(perf_out_MaxProd2);
   ClassifierPerfExp(3).MaxProd =gather(perf_out_MaxProd3);
   
   ClassifierPerfExp(1).TDOA =gather(perf_out_TDOA1);
   ClassifierPerfExp(2).TDOA =gather(perf_out_TDOA2);
   ClassifierPerfExp(3).TDOA =gather(perf_out_TDOA3);
   
   
   ClassifierPerfExp(1).Baseline =gather(perf_out_baseline1);
   ClassifierPerfExp(2).Baseline =gather(perf_out_baseline2);
   ClassifierPerfExp(3).Baseline =gather(perf_out_baseline3);
    

aa_error = 1 -ClassifierPerfExp(1).MaxProd;
bb_error = 1 -ClassifierPerfExp(2).MaxProd;
cc_error = 1 -ClassifierPerfExp(3).MaxProd; 
UnaidError1 = 1- UnAidedPerf1;
UnaidError2 = 1- UnAidedPerf2;
UnaidError3 = 1- UnAidedPerf3;

% Method one change inerror
DeltErr1= mean(bsxfun(@minus,...
    reshape(UnaidError1,[1,1, nIters]), aa_error)./aa_error, 3);

DeltErr2= mean(bsxfun(@minus,...
    reshape(UnaidError2,[1,1, nIters]), bb_error)./bb_error, 3);

DeltErr3= mean(bsxfun(@minus,...
    reshape(UnaidError3,[1,1, nIters]), cc_error)./cc_error, 3);


figure;
axisVals =[min([DeltErr1(:); DeltErr2(:); DeltErr3(:)]),...
    max([DeltErr1(:); DeltErr2(:); DeltErr3(:)])];
subplot(2,2,1)
imagesc(SimThresh,TimeThresh, DeltErr1),axis xy , colorbar
ylabel('Time Threshold (s)'); xlabel('Similarity Threshold');
caxis(axisVals)


title('Spatial Method- Low Q Classifier')
subplot(2,2,2)
imagesc(SimThresh,TimeThresh,DeltErr2),axis xy , colorbar
ylabel('Time Threshold (s)'); xlabel('Similarity Threshold');
caxis(axisVals)


title('Spatial Method- Mod Q Classifier')
subplot(2,2,3)
imagesc(SimThresh,TimeThresh,DeltErr3),axis xy , colorbar
ylabel('Time Threshold (s)'); xlabel('Similarity Threshold');
title('Spatial Method- High Q Classifier')
caxis(axisVals)



aa_error = 1 -ClassifierPerfExp(1).TDOA;
bb_error = 1 -ClassifierPerfExp(2).TDOA;
cc_error = 1 -ClassifierPerfExp(3).TDOA; 
UnaidError1 = 1- UnAidedPerf1;
UnaidError2 = 1- UnAidedPerf2;
UnaidError3 = 1- UnAidedPerf3;

% Method one change inerror
DeltErr1= mean(bsxfun(@minus,...
    reshape(UnaidError1,[1,1, nIters]), aa_error)./aa_error, 3);

DeltErr2= mean(bsxfun(@minus,...
    reshape(UnaidError2,[1,1, nIters]), bb_error)./bb_error, 3);

DeltErr3= mean(bsxfun(@minus,...
    reshape(UnaidError3,[1,1, nIters]), cc_error)./cc_error, 3);


axisVals =[min([DeltErr1(:); DeltErr2(:); DeltErr3(:)]),...
    max([DeltErr1(:); DeltErr2(:); DeltErr3(:)])];

figure;
subplot(2,2,1)
imagesc(SimThresh,TimeThresh,DeltErr1),axis xy , colorbar
ylabel('Time Threshold (s)'); xlabel('Similarity Threshold');
caxis(axisVals)
title('TDOA Method- Low Q Classifier')

subplot(2,2,2)
imagesc(SimThresh,TimeThresh,DeltErr2),axis xy , colorbar
ylabel('Time Threshold (s)'); xlabel('Similarity Threshold');
title('TDOA Method- Mod Q Classifier')
subplot(2,2,3)
imagesc(SimThresh,TimeThresh,DeltErr3),axis xy , colorbar
caxis(axisVals)
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




%% Clock Drift Experiment

simStruct.c = 1500;
simStruct.betaParm1 = betaParm1(1);
simStruct.betaParm2=betaParm2(1);
simStruct.hydrophone_struct = hydrophone_struct;
clock_drift = [1 2 3 4];
association_sec = [0 0 4 6];

perf_assocExp =struct();
nIters = 4;

perf_out_baseline = struct();
perf_out_TDOA = struct();
perf_out_MaxProd = struct();


for ii=1:length(association_sec)
    
    assoc = association_sec(ii);
    cdrift = clock_drift(ii);
    perf_TDOA = zeros(length(TimeThresh), length(SimThresh), nIters);
   
    perf_Max = perf_TDOA;
    Baseline_perf = zeros(length(TimeThresh), nIters);
    UnaidedPerfRes = zeros([1 nIters]);
    
    parfor kk = 1:nIters
        % Replace the space whale component
        [spaceWhale] =   createRandomSpaceWhale(0.5, 6, hyd_arr,...
            array_struct,hydrophone_struct, ssp, grid_depth,...
            [array_struct.master, array_struct.slave(child_idx)]);
        

        % Copy for each of the methods
        simStructNew = simStruct;
        simStructNew.spaceWhale = spaceWhale;
        simStructNew.betaParm1=5;
        simStructNew.betaParm2=2;
        simStructNew.drift=cdrift;
        simStructNew.assSec=assoc;
        simStructNew.arrivalTable = gather(UpdateArrTable(simStructNew));
        simStructNew.arrivalArray= gather(UpdateArrArray(simStructNew));
        simStructNew.TDOA_vals = gather(UpdateTDOA(simStructNew));
        simStructNew.maxEltTime = gather(max(TimeThresh));
        simStructNew.truthTable = createTruthTable(simStructNew);
        simStructNew.truthTable=createSpeciesPreds(simStructNew);
        
        
        
        % Baseline
        simStructBaseline =simStructNew;
        

        % TDOA
        simStructTDOA =simStructNew;
        simStructTDOA.Sim_mat = simMatTDOAonly(simStructTDOA);
        
        
        % Max of mean
        simStructMaxProd = simStructNew;
        simStructMaxProd.Sim_mat = simMatMaxofProd(simStructNew);
        
        [aa bb cc, UnAidedPerf] = clockDriftLpFx(simStructMaxProd, simStructTDOA,simStructBaseline,...
            TimeThresh, SimThresh);
        
        
        Baseline_perf(:,kk) = aa;
        perf_TDOA(:,:,kk) =bb;
        perf_Max(:,:,kk) = cc;
        UnaidedPerfRes(kk)=unique(UnAidedPerf);
        
        
    end
    perf_assocExp(ii).Baseline = Baseline_perf;
    perf_assocExp(ii).TDOA = perf_TDOA;
    perf_assocExp(ii).MaxProd = perf_Max;
    perf_assocExp(ii).UnaidedPerf =UnaidedPerfRes;
end


%% Plot Clock Drift Experiment Results
clock_drift = [2 2 3 4];
association_sec = [0 3 4 6];

for jj=1:length(clock_drift)
    perfMaxProd =[];
    perfTDOA =[];
    perfBaseline=[];
    

   
    UnaidError1 =1- perf_assocExp(jj).UnaidedPerf;
    
    aa_error = 1 -perf_assocExp(jj).Baseline;
    DeltErr1= mean(bsxfun(@minus,...
        reshape(UnaidError1,[1, nIters]), aa_error)./aa_error, 3);
    
    DeltErr1 = 1-quantile(DeltErr1,[ .5 .25 .75],2);
    

    
    
    figure
    subplot(2,2,1)
    plot_ci(TimeThresh, DeltErr1, 'PatchColor', 'r', 'PatchAlpha', 0.2, ...
        'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', 'b', ...
        'LineWidth', 1.5, 'LineStyle','--', 'LineColor', 'k')
    xlabel('Time Threshold (s)'); ylabel('Change in Accuracy')
    title('Baseline Method- Low Q Classifier')
    
    
    
    bb_error = 1 -perf_assocExp(jj).TDOA;
    cc_error = 1 -perf_assocExp(jj).MaxProd;
    
    UnaidError1 = 1- UnAidedPerf1;
    
    % Method one change inerror
    
    DeltErr2= mean(bsxfun(@minus,...
        reshape(UnaidError1,[1,1, nIters]), bb_error)./bb_error, 3);
    
    DeltErr3= mean(bsxfun(@minus,...
        reshape(UnaidError1,[1,1, nIters]), cc_error)./cc_error, 3);
    
    subplot(2,2,2)
    imagesc(SimThresh,TimeThresh, DeltErr2), axis xy
    ylabel('Time Threshold (s)'); xlabel('Similarity Threshold');
    title([num2str(clock_drift(jj)), ' s Clock Drift', num2str(association_sec(jj)),...
        ' s Association- TDOA'])
    
    subplot(2,2,3)
    imagesc(SimThresh,TimeThresh, DeltErr3), axis xy
    ylabel('Time Threshold (s)'); xlabel('Similarity Threshold');
    title([num2str(clock_drift(jj)), ' s Clock Drift', num2str(association_sec(jj)),...
        ' s Association- MaxProd'])

    
    
    
        
    
end




