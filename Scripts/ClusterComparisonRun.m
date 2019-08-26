
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


% Trim the LSQ space
localize_struct = trimlocalize_struct(localize_struct, .5, 5)

localize_struct = trimlocalize_struct(localize_struct, .5, 8)


%% Pick a parent and go

fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;






%% Load truth dataset
% Find the indicies that coinced with the truth table
% Load the calls
%truth = readtable('C:\Users\Kaitlin\Desktop\NOPP6_EST_20090328.selections.txt');
truth = readtable('/home/kpalmer/Desktop/NOPP6_EST_20090328.selections.txt');
truth.mstart = datenum('20090328', 'yyyymmdd')+truth.BeginTime_s_/86400;
truth.mend = datenum('20090328', 'yyyymmdd')+truth.EndTime_s_/86400;
truth.mid = (truth.mend+truth.mstart)/2;


%% Error Rates - channel 5 as parent
% Populate data and parameters
parent =8

calls_arrivals =struct2table(hyd(parent).detection.calls);
examp = struct();
examp.array_struct = localize_struct.hyd(parent).array;
examp.hydrophone_struct = hydrophone_struct;
examp.randomMiss =0;
examp.s = 12;
examp.child_idx = [1:8];
examp.localize_struct =localize_struct;
examp.limitTime =3*60*60;
examp.maxEltTime = max(TimeThresh);
examp.calls =calls_arrivals;
examp.callParm =  hyd(parent).detection.parm;
examp.fs =2000;
examp.PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
examp.drift = 3;
examp.c =1500;
examp.cutoff = .2;

examp.arrivalTable = updateArrTableGPL(examp); % Run this first!
examp.arrivalArray = UpdateArrArray(examp);
examp.TDOA_vals = UpdateTDOA(examp);

% Spatial Method
exampSpatial = examp;
exampSpatial.Sim_mat= simMatMaxofProd(exampSpatial);


% Do the same for the TDOA only method
exampTDOA = examp;
exampTDOA.Sim_mat= simMatTDOAonly(exampTDOA);

% TOA/acoustic encounters only
exampBaseline = examp;

%% Plot Results

simThreshs = linspace(0.01,.99, 30);
TimeThresh = fliplr(linspace(2,60*2,20));

Propout8 =[];
NClustout8 =[];


PropoutTDOA8 =[];
NClustoutTDOA8 =[];


PropoutBaseline8 =[];
NClustoutBaseline8 =[];



for jj = 1:length(simThreshs)
    for ii=1:length(TimeThresh)
        
        exampSpatial.cutoff = simThreshs(jj);
        exampSpatial.maxEltTime = TimeThresh(ii);
        chains =updateChainsEncounterFirst(exampSpatial);
        exampSpatial.chains = chains;
        exampSpatial.Cluster_id= updateClusterID(exampSpatial);
        [propCorrect, NClusters] = RunComp(exampSpatial, parent, hyd, truth,simThresh);
        Propout8(ii,jj)=propCorrect;
        NClustout8(ii,jj)=length(exampSpatial.Cluster_id)/length(exampSpatial.chains);
        
        
        exampTDOA.maxEltTime = TimeThresh(ii);
        exampTDOA.cutoff = simThreshs(jj);
        chains =updateChainsEncounterFirst(exampTDOA);
        exampTDOA.chains = chains;
        exampTDOA.Cluster_id= updateClusterID(exampTDOA);
        [propCorrect, NClusters] = RunComp(exampTDOA, parent, hyd, truth,simThresh);
        PropoutTDOA8(ii,jj)=propCorrect;
        NClustoutTDOA8(ii,jj)= length(exampTDOA.Cluster_id)/length(exampTDOA.chains);
        
        
        exampBaseline.maxEltTime= TimeThresh(ii);
        exampBaseline.cutoff = simThreshs(jj);
        exampBaseline.Cluster_id = acEnc(exampBaseline);
        [propCorrectBase, NClusters] = RunComp(exampBaseline, parent, hyd, truth, simThresh);
        PropoutBaseline8 = [PropoutBaseline8,propCorrectBase];
        meanClusterSize = length(exampBaseline.Cluster_id)/length(unique(exampBaseline.Cluster_id));
        NClustoutBaseline8(ii,jj)= meanClusterSize;
        
    end
end




figure(3)
subplot(2,1,2)
hold on
plot((TimeThresh), (PropoutTDOA8), '-+')
plot((TimeThresh), (Propout8), '-*')
plot((TimeThresh), (PropoutBaseline8), '-k')
xlabel('Elapsed Time Between Acoustic Encounters')
ylabel('Proportion of Calls Correctly Classified')
title('Hydrophone 8')
legend('TDOA Method', 'Spatial Method', 'Baseline')
ylim([.75 1])





%% Error Rates - channel 5 as parent
% Populate data and parameters
TimeThresh = fliplr(linspace(2,60*2,10));
parent =5

calls_arrivals =struct2table(hyd(parent).detection.calls);
examp = struct();
examp.array_struct = localize_struct.hyd(parent).array;
examp.hydrophone_struct = hydrophone_struct;
examp.randomMiss =0;
examp.s = 12;
examp.child_idx = [1:8];
examp.localize_struct =localize_struct;
examp.limitTime =3*60*60;
examp.maxEltTime = max(TimeThresh);
examp.calls =calls_arrivals;
examp.callParm =  hyd(parent).detection.parm;
examp.fs =2000;
examp.PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
examp.drift = 3;
examp.c =1500;
examp.arrivalTable = updateArrTableGPL(examp); % Run this first!
examp.arrivalArray = UpdateArrArray(examp);
examp.TDOA_vals = UpdateTDOA(examp);

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


for ii=1:length(TimeThresh)
    
    exampSpatial.cutoff = simThresh;
    exampSpatial.maxEltTime = TimeThresh(ii);
    chains =updateChainsEncounterFirst(exampSpatial);
    exampSpatial.chains = chains;
    exampSpatial.Cluster_id= updateClusterID(exampSpatial);
    [propCorrect, NClusters] = RunComp(exampSpatial, parent, hyd, truth,simThresh);
    Propout5(ii)=propCorrect;
    NClustout5(ii)=length(exampSpatial.Cluster_id)/length(exampSpatial.chains);
    

    exampTDOA.maxEltTime = TimeThresh(ii);
    exampTDOA.cutoff = simThresh;
    chains =updateChainsEncounterFirst(exampTDOA)
    exampTDOA.chains = chains;
    exampTDOA.Cluster_id= updateClusterID(exampTDOA);
    [propCorrect, NClusters] = RunComp(exampTDOA, parent, hyd, truth,simThresh);
    PropoutTDOA5(ii)=propCorrect;
    NClustoutTDOA5(ii)= length(exampTDOA.Cluster_id)/length(exampTDOA.chains);


    exampBaseline.maxEltTime= TimeThresh(ii);
    exampBaseline.cutoff = simThresh;
    exampBaseline.Cluster_id = acEnc(exampBaseline);
    [propCorrectBase, NClusters] = RunComp(exampBaseline, parent, hyd, truth, simThresh);
    PropoutBaseline5 = [PropoutBaseline8,propCorrectBase];
    meanClusterSize = length(exampBaseline.Cluster_id)/length(unique(exampBaseline.Cluster_id));
    NClustoutBaseline5= [NClustoutBaseline8, meanClusterSize];

end




figure(3)
subplot(2,1,1)
hold on
plot((TimeThresh), (PropoutTDOA8), '-+')
plot((TimeThresh), (Propout8), '-*')
plot((TimeThresh), (PropoutBaseline8), '-k')
xlabel('Elapsed Time Between Acoustic Encounters')
ylabel('Proportion of Calls Correctly Classified')
title('Hydrophone 5')
legend('TDOA Method', 'Spatial Method', 'Baseline')
ylim([.75 1])





figure(3)
subplot(2,1,2)
hold on
legend('TDOA Method', 'Spatial Method', 'Baseline')
ylim([.75 1])














