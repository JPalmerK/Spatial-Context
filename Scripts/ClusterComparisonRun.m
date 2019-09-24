
% Run experiments
close all; clear all
clear classes; clc

whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2})
load(whereAmI{3})
load(whereAmI{4})
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

% Merge the array structure with the localize structure

for i = 1:length(localize_struct.hyd)
    localize_struct.hyd(i).array = array_struct_data(i).array;
end


fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;

% number of true positive RW dets in the untrimmed dataset (206 hyd 5 and 348 hyd 8)
TotalTp = 348;
parent = 8;
examp = struct();
examp.hydrophone_struct = hydrophone_struct;
examp.randomMiss =0;
examp.s = 8;
examp.child_idx = [1:8];
examp.fs =2000;
examp.PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
examp.drift = 2;
examp.c =1500;
examp.truncateKm=15;
examp.array_struct = localize_struct.hyd(parent).array;

% Create filter grids to knock out ambiguity surfaces
examp.filtGrid = createDetRangeFiltGrid(examp, hydrophone_struct);

% Remove all the delays where the TDOA is greater than the array geometry
localize_struct1 = trimlocalize_structTDOA(localize_struct, hyd)

% Remove calls where cross correlation is 1
localize_struct_temp = trimlocalize_structCrossCorr(localize_struct1, hyd)


%% Load truth dataset
% Find the indicies that coinced with the truth table
% Load the calls
truth = readtable('C:\Users\Kaitlin\Desktop\NOPP6_EST_20090328.selections.txt');
%truth = readtable('/home/kpalmer/Desktop/NOPP6_EST_20090328.selections.txt');
truth.mstart = datenum('20090328', 'yyyymmdd')+truth.BeginTime_s_/86400;
truth.mend = datenum('20090328', 'yyyymmdd')+truth.EndTime_s_/86400;
truth.mid = (truth.mend+truth.mstart)/2;

% Get the total number of detections and TP for precision/recall


%% Create the structure
simThreshs = linspace(0.01,.99, 10);
TimeThreshs = fliplr(linspace(2,20,5));
corrThresh = fliplr(linspace(0.1,.8,15));

ExperimentTDOA = struct();
ExperimentMaxProd = struct();
ExperimentUnaided = struct();

% Pre allocate the output matrix
PrecisionMatSpatial = zeros(length(simThreshs), length(TimeThreshs), length(corrThresh))/0;
RecallMatSpatial = PrecisionMatSpatial;
PrecisionMatTDOA = PrecisionMatSpatial;
RecallMatTDOA = PrecisionMatSpatial;
PrecisionMatBaseline =PrecisionMatSpatial;
RecallMatBaseline =PrecisionMatSpatial;
PrecisionMatAcEnc= PrecisionMatSpatial;
RecallMatAcEnc =PrecisionMatSpatial;
% Create the basic structure, this needs to be done for every correlation
% threshold


%%
exampNew=examp;
localize_struct_trimmed = trimlocalize_struct(localize_struct_temp, hyd, 0.1);


examp= exampNew;
% Clip the detections by the correlation threshold
examp.localize_struct =localize_struct_trimmed;
examp.maxEltTime = max(TimeThreshs);
examp.callParm =  hyd(parent).detection.parm;
examp.arrivalTable = updateArrTableGPL(examp); % Run this first!
examp.arrivalArray = UpdateArrArrayRealData(examp);
examp.TDOA_vals = UpdateTDOA(examp);
examp.maxEltTime = max(TimeThreshs);

exampSpatial = examp;
exampTDOA = examp;
exampBaseline = examp;

exampSpatial.Sim_mat= simMatMaxofProd(exampSpatial);
exampTDOA.Sim_mat = simMatTDOAonly(examp);
exampBaseline.Cluster_id = acEnc(examp);

Results=struct();

for ii = 1:length(simThreshs)
    simThresh = simThreshs(ii);
    exampSpatial.cutoff = simThresh;
    exampTDOA.cutoff = simThresh;
    
    for jj = 1:length(TimeThreshs)
        % Create precision recall curves for each similarity and time threshold
        
            exampSpatial.maxEltTime = TimeThreshs(jj);
            exampTDOA.maxEltTime = TimeThreshs(jj);
            exampBaseline.maxEltTime = TimeThreshs(jj);
            
            
            % Run the clustering spatial
            exampSpatial.chains = updateChainsEncounterFirst(exampSpatial);
            exampSpatial.Cluster_id= updateClusterID(exampSpatial);
            
            %Run the clustering TDOA
            exampTDOA.chains = updateChainsEncounterFirst(exampTDOA);
            exampTDOA.Cluster_id= updateClusterID(exampTDOA);
            
            % Baseline clustering, time only
            exampBaseline.Cluster_id = acEnc(examp);
            
            % Spatial and baseline
            [ClusteredSpatial, Detector] = PrecisionRecall(exampSpatial, parent, hyd, truth);
            Results.Spatial(ii,jj) = ClusteredSpatial;
            
            if ii == jj ==1
            Results.Baseline(1,1) = Detector;
            end
            
            % TDOA precision/recall
            [ClusteredTDOA, ~] = PrecisionRecall(exampTDOA, parent, hyd, truth);
            Results.TDOA(ii,jj) = ClusteredTDOA;
            
            
            % Baseline clustering precision recall
            [ClusteredBase, ~] = PrecisionRecall(exampBaseline, parent, hyd, truth);
            Results.AcousticEncounters(ii,jj) = ClusteredBase;
        
    end
    disp('Finish Sim Thres Iter')
end

close all

timeidx =4;
figure;
for ii =1:10
    
    titlestr =[num2str(TimeThreshs(timeidx)), ' sec TimeThresh ',...
        num2str(num2str(simThreshs(ii),'%0.2f')), ' SimThresh'];
    
    subplot(3,4,ii)
    
    hold on
    plot(Results.Baseline(1,1).Recall, Results.Baseline(1,1).Precision, '.-')
    plot(Results.TDOA(ii,timeidx).Recall, Results.TDOA(ii,timeidx).Precision, '.-')
    plot(Results.Spatial(ii,timeidx).Recall, Results.Spatial(ii,timeidx).Precision, '.-')
    plot(Results.AcousticEncounters(ii,timeidx).Recall,...
        Results.AcousticEncounters(ii,timeidx).Precision, '.-')
    xlabel('Recall'); ylabel('Precision');
    %legend('Detector Only', 'TDOA', 'Spatial', 'Acoustic Encounters Only')
    title(titlestr)


end






disp('Finish Corr Iter')

figure; hold on;
plot(squeeze(RecallMatBaseline(5,3,:)), squeeze(PrecisionMatBaseline(5,3,:)), '.-')
plot(squeeze(RecallMatTDOA(5,3,:)), squeeze(PrecisionMatTDOA(5,3,:)),'.-')
plot(squeeze(RecallMatAcEnc(5,3,:)), squeeze(PrecisionMatAcEnc(5,3,:)), '.-')
plot(squeeze(RecallMatSpatial(5,3,:)), squeeze(PrecisionMatSpatial(5,3,:)),'.-')
legend('Baseline', 'TDOA', 'Time Only Clustering', 'Spatial')
ylabel('Precision'); xlabel('Recall')

%% Error Rates - channel 5 as parent
% Populate data and parameters
corrThresh =0.5;

% Make the summary table
GPLSummary = GPLRavenSummaryComp(hyd, truth, corrThresh, 5)

GPLSummary = GPLRavenSummaryComp(localize_struct, hyd, truth, corrThresh, 5)


% Trim the localize structure by correlation score
hydTrimmed = trimCallsCorrThresh(hyd, corrThresh)
localize_struct_trimmed = trimlocalize_struct(localize_struct,...
    hydTrimmed, corrThresh);

%calls_arrivals =struct2table(hyd(parent).detection.calls);

%%
examp = struct();
examp.array_struct = localize_struct_trimmed.hyd(parent).array;
examp.hydrophone_struct = hydrophone_struct;
examp.randomMiss =0;
examp.s = 8;
examp.child_idx = [1:8];
examp.localize_struct =localize_struct_trimmed;
examp.limitTime =24*60*60;
examp.maxEltTime = max(TimeThresh);
%examp.calls =calls_arrivals;
examp.callParm =  hyd(parent).detection.parm;
examp.fs =2000;
examp.PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
examp.drift = 2;
examp.c =1500;

examp.truncateKm=15;

examp.arrivalTable = updateArrTableGPL(examp); % Run this first!
examp.arrivalArray = UpdateArrArrayRealData(examp);
examp.TDOA_vals = UpdateTDOA(examp);


% Create filter grids to knock out ambiguity surfaces
examp.filtGrid = createDetRangeFiltGrid(examp, hydrophone_struct);

% Spatial Method
exampSpatial = examp;
exampSpatial.Sim_mat= simMatMaxofProd(exampSpatial);


% Do the same for the TDOA only method
exampTDOA = examp;
exampTDOA.Sim_mat= simMatTDOAonly(exampTDOA);

% TOA/acoustic encounters only
exampBaseline = examp;

%% Plot Results


Propout8 =[];
NClustout8 =[];


PropoutTDOA8 =[];
NClustoutTDOA8 =[];


PropoutBaseline8 =[];
NClustoutBaseline8 =[];



for jj = 1:length(simThreshs)
    for ii=1:length(TimeThresh)
        simThresh = simThreshs(jj);
        
        exampSpatial.cutoff = simThresh;
        exampSpatial.maxEltTime = TimeThresh(ii);
        exampSpatial.chains = updateChainsEncounterFirst(exampSpatial);
        exampSpatial.Cluster_id= updateClusterID(exampSpatial);
        [propCorrect, NClusters] = RunComp(exampSpatial, parent, hyd, truth,simThresh);
        Propout8(ii,jj)=propCorrect;
        NClustout8(ii,jj)=length(exampSpatial.Cluster_id)/length(exampSpatial.chains);
        
        
        exampTDOA.maxEltTime = TimeThresh(ii);
        exampTDOA.cutoff = simThresh;
        exampTDOA.chains = updateChainsEncounterFirst(exampTDOA);
        exampTDOA.Cluster_id= updateClusterID(exampTDOA);
        [propCorrect, NClusters] = RunComp(exampTDOA, parent, hyd, truth, simThresh);
        PropoutTDOA8(ii,jj)=propCorrect;
        NClustoutTDOA8(ii,jj)= length(exampTDOA.Cluster_id)/length(exampTDOA.chains);
        
        
        exampBaseline.maxEltTime= TimeThresh(ii);
        exampBaseline.cutoff = simThresh;
        exampBaseline.Cluster_id = acEnc(exampBaseline);
        [propCorrectBase, NClusters] = RunComp(exampBaseline, parent, hyd, truth, simThresh);
        PropoutBaseline8(ii,jj) = propCorrectBase;
        meanClusterSize = length(exampBaseline.Cluster_id)/length(unique(exampBaseline.Cluster_id));
        NClustoutBaseline8(ii,jj)= meanClusterSize;
        
    end
    jj
end

figure (3)
subplot(3,2,1)
imagesc(simThreshs, TimeThresh, Propout8), axis xy, colorbar
title('Spatial Method Hyd 8')
xlabel('Similarity Threshold')
ylabel('Time Threshold')


subplot(3,2,3)
imagesc(simThreshs, TimeThresh, PropoutTDOA8), axis xy, colorbar
xlabel('Similarity Threshold')
ylabel('Time Threshold')
title('TDOA Method Hyd 8')

subplot(3,2,5)
imagesc(simThreshs, TimeThresh, PropoutBaseline8), axis xy, colorbar
xlabel('Similarity Threshold')
ylabel('Time Threshold')
title('Baseline Method Hyd 8')


%% Error Rates - channel 5 as parent
% Populate data and parameters
parent = 5;
calls_arrivals =struct2table(hyd(parent).detection.calls);
examp = struct();
examp.array_struct = localize_struct_trimmed.hyd(parent).array;
examp.hydrophone_struct = hydrophone_struct;
examp.randomMiss =0;
examp.s = 8;
examp.child_idx = [1:8];
examp.localize_struct =localize_struct_trimmed;
examp.limitTime =3*60*60;
examp.maxEltTime = max(TimeThresh);
examp.calls =calls_arrivals;
examp.callParm =  hyd(parent).detection.parm;
examp.fs =2000;
examp.PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
examp.drift = 2;
examp.c =1500;

examp.arrivalTable = updateArrTableGPL(examp); % Run this first!
examp.arrivalArray = UpdateArrArrayRealData(examp);
examp.TDOA_vals = UpdateTDOA(examp);
examp.truncateKm = 15;
% Create filter grids to knock out ambiguity surfaces
examp.filtGrid = createDetRangeFiltGrid(examp, hydrophone_struct);

% Spatial Method
exampSpatial = examp;
exampSpatial.Sim_mat= simMatMaxofProd(exampSpatial);


% Do the same for the TDOA only method
exampTDOA = examp;
exampTDOA.Sim_mat= simMatTDOAonly(exampTDOA);

% TOA/acoustic encounters only
exampBaseline = examp;

%% Plot Results


Propout5 =[];
NClustout5 =[];

PropoutTDOA5 =[];
NClustoutTDOA5 =[];

PropoutBaseline5 =[];
NClustoutBaseline5 =[];



for jj = 1:length(simThreshs)
    
    simThresh =simThreshs(jj);
    
    Propout5_row =zeros(size(TimeThresh,1));
    NClustout5_row = Propout5_row;
    PropoutTDOA5_row = Propout5_row;
    NClustoutTDOA5_row =Propout5_row;
    PropoutBaseline5_row =Propout5_row;
    NClustoutBaseline5_row =Propout5_row;
    
    
    parfor ii=1:length(TimeThresh)
        exampSpatial1 =exampSpatial;
        exampSpatial1.cutoff = simThresh;
        exampSpatial1.maxEltTime = TimeThresh(ii);
        exampSpatial1.chains = updateChainsEncounterFirst(exampSpatial1);
        exampSpatial1.Cluster_id= updateClusterID(exampSpatial1);
        [propCorrect, ~] = RunComp(exampSpatial1, parent, hyd, truth,simThresh);
        Propout5_row(ii)=propCorrect;
        NClustout5_row(ii)=length(exampSpatial1.Cluster_id)/length(exampSpatial1.chains);
        
        exampTDOA1=exampTDOA;
        exampTDOA1.maxEltTime = TimeThresh(ii);
        exampTDOA1.cutoff = simThresh;
        exampTDOA1.chains = updateChainsEncounterFirst(exampTDOA1);
        exampTDOA1.Cluster_id= updateClusterID(exampTDOA1);
        [propCorrect, ~] = RunComp(exampTDOA1, parent, hyd, truth,simThresh);
        PropoutTDOA5_row(ii)=propCorrect;
        NClustoutTDOA5_row(ii)= length(exampTDOA1.Cluster_id)/length(exampTDOA1.chains);
        
        exampBaseline1=exampBaseline;
        exampBaseline1.maxEltTime= TimeThresh(ii);
        exampBaseline1.cutoff = simThresh;
        exampBaseline1.Cluster_id = acEnc(exampBaseline1);
        [propCorrectBase, ~] = RunComp(exampBaseline1, parent, hyd, truth, simThresh);
        PropoutBaseline5_row(ii) = propCorrectBase;
        meanClusterSize = length(exampBaseline1.Cluster_id)/length(unique(exampBaseline1.Cluster_id));
        NClustoutBaseline5_row(ii)= meanClusterSize;
        
    end
    Propout5(jj,:)= Propout5_row;
    NClustout5(jj,:)=NClustout5_row;
    PropoutTDOA5(jj,:)=PropoutTDOA5_row;
    NClustoutTDOA5(jj,:)=NClustoutTDOA5_row;
    PropoutBaseline5(jj,:)=PropoutBaseline5_row;
    NClustoutBaseline5(jj,:)=NClustoutBaseline5_row;
    
end




%%

figure (3)
subplot(3,2,2)
imagesc(simThreshs, TimeThresh, Propout5), axis xy, colorbar
xlabel('Similarity Threshold')
ylabel('Time Threshold')
title('Spatial Method Hyd 8')


subplot(3,2,4)
imagesc(simThreshs, TimeThresh, PropoutTDOA5), axis xy, colorbar
xlabel('Similarity Threshold')
ylabel('Time Threshold')
title('TDOA Method Hyd 8')

subplot(3,2,6)
imagesc(simThreshs, TimeThresh, PropoutBaseline5), axis xy, colorbar
xlabel('Similarity Threshold')
ylabel('Time Threshold')
title('Baseline Method Hyd 8')


%
% figure(3)
% subplot(2,1,1)
% hold on
% plot((TimeThresh), (PropoutTDOA8), '-+')
% plot((TimeThresh), (Propout8), '-*')
% plot((TimeThresh), (PropoutBaseline8), '-k')
% xlabel('Elapsed Time Between Acoustic Encounters')
% ylabel('Proportion of Calls Correctly Classified')
% title('Hydrophone 5')
% legend('TDOA Method', 'Spatial Method', 'Baseline')
% ylim([.75 1])
%
%
%
%
%
% figure(3)
% subplot(2,1,2)
% hold on
% legend('TDOA Method', 'Spatial Method', 'Baseline')
% ylim([.75 1])














