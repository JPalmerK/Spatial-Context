% Run experiments
close all; clear all
clear classes; clc

cd('/home/kpalmer/AnacondaProjects/Localisation/Scripts')
% Load metadata to create Hydrophone Structur
dclde_2013_meta = xlsread(strcat('/cache/kpalmer/quick_ssd/data/',...
    'DCLDE_2013_10Channel/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx'));
% Load the DCLDE array structure (just for the array structure, depth, and SSP)
load('/home/kpalmer/AnacondaProjects/Localisation/Scripts/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_array_struct1_671.mat')
load('/home/kpalmer/AnacondaProjects/Localisation/Scripts/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_localize_struct_671.mat')


% % Load metadata to create Hydrophone Structur
% dclde_2013_meta = xlsread('C:\Users\Kaitlin\Desktop\DCL2013_NEFSC_SBNMS_200903_metadata.xlsx');
% load('D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_array_struct1_671.mat')
% load('D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_localize_struct_671.mat')


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


%% Example

% Create the agents

% Number of hours the experiment runs and the number of agents included
n_hrs = 0.75;
n_agents = 7;

[spaceWhale] =  createRandomSpaceWhale(n_hrs, n_agents, hyd_arr,...
    array_struct,hydrophone_struct, ssp, grid_depth);
% Create the prob loc spaces for all calls


% Create the example class
examp = simulationClass();
examp.spaceWhale = spaceWhale;
examp.array_struct = array_struct;
examp.hydrophone_struct = hydrophone_struct;
examp.spaceWhale= spaceWhale;
examp.time_cut = 15*60;
examp.randomMiss = 0;


% Baseline
toaOnlyCluster(examp)
examp.getRand()
examp.AdjRand

examp.clearCalcValues()
simMatIdeal(examp)
examp.getRand()
examp.AdjRand

examp.clearCalcValues()
simMatadHoc(examp)
examp.getRand()
examp.AdjRand

examp.clearCalcValues()
simMatTDOAonly(examp)
examp.getRand()
examp.AdjRand



%%

UpdateprojSpace(examp)

% Figure of LSQ spaces
figure
s1 = subplot(2,3,1);
imagesc(array_struct.longrid, array_struct.latgrid,...
    (squeeze(examp.projSpace(1).projection(:,:,1))))
caxis([0 .4])
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')
scatter(hyd_arr([5, 1, 2],2), hyd_arr([5, 1, 2],1), 80, 'r', 'filled', 'd')
title_str= ['Call 1'];
set(gca,'XTickLabel',[]);
ylabel('Latitude')
title(title_str)

s2 = subplot(2,3,5);
imagesc(array_struct.longrid, array_struct.latgrid,...
    (squeeze(examp.projSpace(1).projection(:,:,2))))
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')
scatter(hyd_arr([5, 1, 2],2), hyd_arr([5, 1, 2],1), 80, 'r', 'filled', 'd')
caxis([0 .4])
title_str= ['Projected ',...
    num2str(round(examp.arrivalArray(2,1)-examp.arrivalArray(1,1))), ' sec'];
title(title_str)
set(gca,'YTickLabel',[]);
xlabel('Longitude')

s3 = subplot(2,3,6);
imagesc(array_struct.longrid, array_struct.latgrid,...
    (squeeze(examp.projSpace(1).projection(:,:,3))))
caxis([0 .4])
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')
scatter(hyd_arr([5, 1, 2],2), hyd_arr([5, 1, 2],1), 80, 'r', 'filled', 'd')
title_str= ['Projected ',...
    num2str(round(examp.arrivalArray(3,1)-examp.arrivalArray(1,1))), ' sec'];
hb = colorbar('location','eastoutside');
title(title_str)
xlabel('Longitude')
set(gca,'YTickLabel',[]);


s4 = subplot(2,3,4);
imagesc(array_struct.longrid, array_struct.latgrid,...
    (squeeze(examp.projSpace(1).projection(:,:,1))))
caxis([0 .4])
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')
scatter(hyd_arr([5, 1, 2],2), hyd_arr([5, 1, 2],1), 80, 'r', 'filled', 'd')

title_str= ['Call 1'];
title(title_str)
ylabel('Latitude')
xlabel('Longitude')

s5 = subplot(2,3,2);
imagesc(array_struct.longrid, array_struct.latgrid,...
    (squeeze(examp.projSpace(6).projection(:,:,1))))
caxis([0 .4])
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')
scatter(hyd_arr([5, 1, 2],2), hyd_arr([5, 1, 2],1), 80, 'r', 'filled', 'd')
title_str= ['Call 6'];
title(title_str)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


s6 = subplot(2,3,3);
imagesc(array_struct.longrid, array_struct.latgrid,...
    (squeeze(examp.projSpace(11).projection(:,:,1))))
caxis([0 .4])
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')
scatter(hyd_arr([5, 1, 2],2), hyd_arr([5, 1, 2],1), 80, 'r', 'filled', 'd')
title_str= ['Call 11'];
title(title_str)
hb = colorbar('location','eastoutside');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

% # Solution:
s1Pos = get(s2,'position');
s2Pos = get(s6,'position');
s2Pos(3:4) = [s1Pos(3:4)];
set(s6,'position',s2Pos);

s1Pos = get(s2,'position');
s2Pos = get(s3,'position');
s2Pos(3:4) = [s1Pos(3:4)];
set(s3,'position',s2Pos);


% Update similarity matrix and create a figure
simMatLowMemory(examp)
val = examp.estClassifierPerf()


% run the example with different similarity threshold
examp.clearCalcValues()
val1 = examp.estClassifierPerf()




%% Look at with and without accurate correlation

nRuns = 1;

% Run it 50 times
exp1Rand=zeros(1,nRuns);
exp0Rand =exp1Rand;
perf_struct_ideal=struct();
perf_struct_randAss=struct();

%Create the agents
n_hrs = .75;
n_agents = 12;


for ii =1:nRuns
    disp(num2str(ii))
    
    close all
    clear spaceWhale
    
    
    % Create new agents
    [spaceWhale] =  createRandomSpaceWhale(n_hrs,n_agents, hyd_arr,...
        array_struct,hydrophone_struct, ssp, grid_depth);
    
    % Populate data and parameters
    exp1 = simulationClass();
    exp1.spaceWhale = spaceWhale;
    exp1.array_struct = array_struct;
    exp1.hydrophone_struct = hydrophone_struct;
    exp1.spaceWhale= spaceWhale;
    exp1.time_cut = 15*60;
    exp1.randomMiss =1;
    
    
    
    % Run experiment until get chains
    simMatLowMemory(exp1)
    updateClusterID(exp1)
    getRand(exp1)
    exp0Rand(ii)= exp1.AdjRand;
    estClassifierPerf(exp1)
    perf_struct_ideal(ii).perfStructure=estClassifierPerf(exp1);
    %exp1.estClassifierPerf()
    
    % Rerun using random association in calls (simulate bad cross
    % correlations)
    runRandom(exp1, 1)
    updateClusterID(exp1)
    getRand(exp1)
    exp1Rand(ii)= exp1.AdjRand;
    perf_struct_randAss(ii).perfStructure=estClassifierPerf(exp1);
%     
end


%% Pull out accuracy values and do some plotting
close all
acc_ideal =zeros(1,nRuns);
acc_worst = zeros(1,nRuns);
acc_unaided = zeros(1,nRuns);

for ii =1:length(perf_struct_ideal)
    
    acc_unaided(ii) = perf_struct_ideal(ii).perfStructure.totCorrect;
    acc_ideal(ii) = perf_struct_ideal(ii).perfStructure.totCorrectClassMeth;
    acc_worst(ii) = perf_struct_randAss(ii).perfStructure.totCorrectClassMeth;
    
end



edges = 0.5:.05:1;

h1 = histcounts(acc_unaided,edges);
h2 = histcounts(acc_ideal,edges);
h3 = histcounts(acc_worst,edges);





figure
b = bar(edges(1:end-1),[h1; h2; h3]', 'FaceColor','flat')

for k = 1:3
    b(k).CData = k;
end


legend('Unaided Classification','Ideal Association', 'Random Association')
xlabel('Percent of Calls Correctly Classified')
ylabel('Number of Model Iterations')


% CDF plots
figure
hold on;
cdfplot(1-acc_unaided)
cdfplot(1-acc_ideal)
cdfplot(1-acc_worst)
legend('Unaided Classification','Ideal Association', 'Random Association')
xlabel('Percentage of Calss Misclassified')
ylabel('Proportion of Model Runs')



figure
hold on
ecdf(1-acc_unaided,'bounds','on');
ecdf(1-acc_ideal,'bounds','on');
ecdf(1-acc_worst,'bounds','on');
legend('Unaided Classification','Ideal Association', 'Random Association')
xlabel('Percent of Calss Misclassified')
ylabel('Proportion of Model Runs')
