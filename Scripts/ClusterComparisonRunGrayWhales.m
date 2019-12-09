
% Code to pre-process Regina's detection files.

% 1) Load the the required files

% ) re-run GPL v3 to get the localization structure
% ) Compare the resulting localization structure to the validated
% ('pruned') detection files to apply validation information

% Re-run localization v3 on Regina's gray whale data
close all; clear all; clc

% Load the array grids, hydrophone structure, and localization paramter file
load('D:\ReginaGrayWhaleDetections\ReginaGPL_Erdata\parameter structures\regina_localize_struct_090116.mat')
load('D:\ReginaGrayWhaleDetections\ReginaGPL_Erdata\parameter structures\regina_array_struct_090116.mat')
load('D:\ReginaGrayWhaleDetections\ReginaGPL_Erdata\parameter structures\regina_hydrophone_struct.mat')


% Modify the localize struct to work on calls only detected by two
% instruments

localize_struct.parm.min_pairs =1;

% Load the detection file
fileLoc='D:\ReginaGrayWhaleDetections\GraniteCanyonNonValid\';
fdir = dir(strcat(fileLoc, '*.mat'));
hyd_out = struct();

for ii=1:length(fdir)

    load(strcat(fileLoc, fdir(ii).name));
    
    test=num2cell(extractfield(hyd(1).detection.calls,'true_start_time'));
    [hyd.detection.calls(1:length(test)).start_time] = test{:};

    test=num2cell(extractfield(hyd(1).detection.calls,'true_end_time'));
    [hyd.detection.calls(1:length(test)).end_time] = test{:};
    hyd_out(ii).detection= hyd.detection;
    clear hyd
end

hyd = hyd_out;
clear hyd_out



parm = hyd(1).detection.parm
% Modify the array structure to make wider appature
latgrid = linspace(36.38, 36.5, 110);
longrid = linspace(-121.99, -121.92, 90);

array_struct.latgrid = latgrid;
array_struct.longrid = longrid;

% Recompute the tdoa grids
[array_struct] = setup_TDOA_grid(hydrophone_struct, array_struct,...
    localize_struct,parm);


[localize_struct] = localize_cross_corr_index_check(array_struct,...
    hyd,localize_struct,1);
[localize_struct] = localize_LSQ_2D(array_struct, ...
    hydrophone_struct,hyd,localize_struct,1);

filename = 'D:\Anaconda Projects\Spatial-Context\Scripts\ReginaLocalizeStructExpanded.mat'
save(filename,'localize_struct','-v7.3')


filename = 'D:\Anaconda Projects\Spatial-Context\Scripts\ReginaArrayStructExpanded.mat'
save(filename,'array_struct','-v7.3')


% %% Get summary statistics on level of interleaveness in reginas data
% 
% 
% % Pruned file loc
% flocValid = 'D:\ReginaGrayWhaleDetections\ReginaGPL_Erdata\Pruned_detection_files\'
% 
% % Load the first detection file with validation (prune) notes
% fdir = dir(strcat(flocValid, '*.mat'));
% 
% for ii=1:length(fdir)
% load(strcat(flocValid, fdir(ii).name));
% detections_valid = struct2table(hyd.detection.calls);
% detections_valid.start_time = detections_valid.julian_start_time;
% detections_valid.end_time = detections_valid.julian_end_time;
% detections_valid = detections_valid(detections_valid.prune ~=0,:);
% 
% 
% detections_valid = sortrows(detections_valid,{'start_time'},{'ascend'});
% 
% 
% end
% 
% % Get rid of false positives


%% Run the system

% Maximum verified calls 105917, number of calls in the detector ~103000

close all
clear all
clc
%load('D:\ReginaGrayWhaleDetections\ReginaGPL_Erdata\parameter structures\regina_array_struct_090116.mat')
%load('D:\Anaconda Projects\Spatial-Context\Scripts\ReginaLocalizeStruct.mat')
%load('D:\Anaconda Projects\Spatial-Context\Scripts\GWArrayStruct.mat')

load('D:\Anaconda Projects\Spatial-Context\Scripts\ReginaLocalizeStructExpanded.mat')
load('D:\Anaconda Projects\Spatial-Context\Scripts\ReginaArrayStructExpanded.mat')
load('D:\ReginaGrayWhaleDetections\ReginaGPL_Erdata\parameter structures\regina_hydrophone_struct.mat')

%%
% Pruned file loc
flocValid = 'D:\ReginaGrayWhaleDetections\ReginaGPL_Erdata\Pruned_detection_files\'

% time (in days) allowable for assoication
wiggleroom = .5/60/60/24;

% Load the first detection file with validation (prune) notes
fdir = dir(strcat(flocValid, '*.mat'));
load(strcat(flocValid, fdir(1).name));
detections_valid = struct2table(hyd.detection.calls);
detections_valid.start_time = detections_valid.julian_start_time;
detections_valid.end_time = detections_valid.julian_end_time;
currentPartitian = 1;
pruned = zeros(length(localize_struct.hyd.dex),1)/0;
mDates = pruned;
det_score = pruned;

% 'Pruned' is what GPL/Localization code uses to denote annotations. A
% detection that has been 'pruned' as type 1 means the user ascribed the
% detection to call type 1. For this work call type 1 is the gray whale M3
% call (what we are after). pruned value 0 means noise and pruned value 2
% are detections annotated as fin whales


% Location of the unverified detection files
flocUnverf = 'D:\ReginaGrayWhaleDetections\GraniteCanyonNonValid\';
det_unverf = load([flocUnverf,...
    'detections_GPL_v2_20_100_GCNE_01_141123_1225_95peroverlap.mat']);

validScoresCount =0;

for ii =1: length(localize_struct.hyd.dex)
    
    % Index of the current localized call from the non-validated
    % hydrophones
    currentDex = localize_struct.hyd.dex(ii);
    
    % Get the arrival time for the non-validated call
    currentTime = det_unverf.hyd.detection.calls(currentDex).julian_start_time;
    
    % check if dex is into the next hydrophone, if so load the data
    if currentDex>max(detections_valid.Entry_Number);
        currentPartitian = currentPartitian + 1;
       
        disp('loading next hydrophone')
        load(strcat(flocValid, fdir(currentPartitian).name));
        
        detections_valid = struct2table(hyd.detection.calls);
        detections_valid.start_time = detections_valid.julian_start_time;
        detections_valid.end_time = detections_valid.julian_end_time;
       
        detections_valid = detections_valid(detections_valid.prune ~=0,:);
        validScoresCount
        validScoresCount=0;
    end
    
    % Start time of the localized call
    vallvals = find(...
        detections_valid.julian_start_time< (currentTime+wiggleroom)  &...
        detections_valid.julian_start_time> (currentTime-wiggleroom));
    
    
    if ~isempty(vallvals)
        
%         if length(vallvals)>1
%             disp('blarg')
%         end
        scores = detections_valid.calltype_1_score(vallvals);
        [scoremax, idx] = max(scores);
        sppid =vallvals(idx);
        pruned(ii)= detections_valid.prune(sppid);
        mDates(ii) = detections_valid.julian_start_time(sppid);
%         % Create fake scores
%         if pruned(ii)==1
%             scoremax = betarnd(8,4);
%         else
%             scoremax = betarnd(4,8);
%         end
        
        det_score(ii) = scoremax;
        validScoresCount=validScoresCount+1;
       
    end

end



localize_struct.hyd.pruned = pruned';
localize_struct.hyd.detectorScore = det_score';
localize_struct.hyd.mdate = mDates;


% Clean out all calls where there are not detector scores
noninterleavedIdx = find(isnan(pruned));
localize_struct.hyd.delays(noninterleavedIdx,:)=[];
localize_struct.hyd.cross_score(noninterleavedIdx,:)=[];
localize_struct.hyd.coord_time(noninterleavedIdx,:)=[];
localize_struct.hyd.rtimes(noninterleavedIdx)=[];
localize_struct.hyd.dex(noninterleavedIdx)=[];
localize_struct.hyd.coordinates(:,:,noninterleavedIdx)=[];
localize_struct.hyd.score(:,noninterleavedIdx)=[];
localize_struct.hyd.detectorScore(noninterleavedIdx)=[];
localize_struct.hyd.pruned(noninterleavedIdx)=[];
localize_struct.hyd.mdate(noninterleavedIdx)=[];


localize_struct = trimlocalize_structCrossCorr(localize_struct, hyd, .8)


labelSpp=[{'Noise'}, {'Gray Whale'}, {'Fin Whale'}]
idxvals =[0 1 2]
for ii=1:length(idxvals)
    
    idx = find(localize_struct.hyd.pruned == idxvals(ii));
    scores = localize_struct.hyd.detectorScore(idx);
    
    subplot(1,3,ii); hist((scores),[0:.05:1]);
    title([labelSpp(ii)]);
    xlabel('Template Detection Scores');
end



%% Precompute the ambiguity surfaces


% Just trim for prelim data
arr_times_sec = localize_struct.hyd(1).rtimes/2000;
noninterleavedIdx = [1:2000,3000:length(localize_struct.hyd.dex)];

localize_struct_trimmed = localize_struct;

localize_struct_trimmed.hyd.delays(noninterleavedIdx,:)=[];
localize_struct_trimmed.hyd.cross_score(noninterleavedIdx,:)=[];
localize_struct_trimmed.hyd.coord_time(noninterleavedIdx,:)=[];
localize_struct_trimmed.hyd.rtimes(noninterleavedIdx)=[];
localize_struct_trimmed.hyd.dex(noninterleavedIdx)=[];
localize_struct_trimmed.hyd.coordinates(:,:,noninterleavedIdx)=[];
localize_struct_trimmed.hyd.score(:,noninterleavedIdx)=[];
localize_struct_trimmed.hyd.detectorScore(noninterleavedIdx)=[];

localize_struct_trimmed.hyd.pruned(noninterleavedIdx)=[];
localize_struct_trimmed.hyd.mdate(noninterleavedIdx)=[];

figure
hold on

loc = cat(1,hydrophone_struct.location)

localize_struct_trimm = trimlocalize_struct(localize_struct, hyd, .2)
x = squeeze(localize_struct_trimmed.hyd.coordinates(1,1,:));
y = squeeze(localize_struct_trimmed.hyd.coordinates(1,2,:));
spp = localize_struct_trimmed.hyd.pruned;
detScore = localize_struct_trimmed.hyd.detectorScore;
figure
categories = {'Noise', 'Gray', 'Fin'};
catnum = [0, 1, 2];
symbol = ['s', 'd', 'v'];  % noise/species symbol
times = cat(1, localize_struct_trimm.hyd.rtimes);
range = minmax(times);
colorsN = 10000;
timestamp = min(colorsN, floor((times - range(1)) / (range(2) - range(1)) * colorsN)+1);
color = hsv(colorsN);


for ii = 1:length(catnum)
    hold on
    catcalls = find(localize_struct_trimmed.hyd.pruned == catnum(ii));
    
    h = scatter(y(catcalls), x(catcalls),40,...
         'filled')
    set(h, 'Marker', symbol(ii))
    

end 
scatter(loc(:,2), loc(:,1), 60,'k', 'filled')
xlabel('Longitude')
ylabel('Latitude')
legend('Noise', 'Gray', 'Fin', 'hydrophones')



%%

fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;

parent = 1;
examp = struct();
examp.hydrophone_struct = hydrophone_struct;
examp.randomMiss =0;
examp.s = 2;
examp.child_idx = [1:3]; % local array for 5, [1,2,3,4,6], for 8 master 4:8
examp.fs =2000;
examp.PosUncertsigma = 0.0004^2 +.01^2 + .03^2; % seconds see EM Nosal 20; % seconds see EM Nosal 20
examp.drift = .01;
examp.c =1500; 
examp.truncateKm=10;
examp.array_struct = array_struct;
examp.propAmbLoc = 'D:\ReginaGrayWhaleDetections\PreComputedAmbSurfs\'; 
%location of the ambiguity surfaces



examp.localize_struct =localize_struct_trimmed;
examp.maxEltTime =45*60;
examp.callParm =  hyd(parent).detection.parm;
examp.arrivalTable = updateArrTableGPL(examp); % Run this first!
examp.arrivalArray = UpdateArrArrayRealData(examp);
examp.TDOA_vals = UpdateTDOA(examp);
% Create filter grids to knock out ambiguity surfaces where calls can't
% have come from
examp.filtGrid = createDetRangeFiltGrid(examp, hydrophone_struct);


exampSpatial = examp;
exampSpatial.Sim_mat= simMatMaxofProdPreComputed(examp);
% Normalize to use geometric mean
exampBaseline = examp;

exampTDOA= examp;
exampTDOA.Sim_mat = simMatTDOAonly(examp);
exampBaseline.Sim_mat = ones(size(exampTDOA.Sim_mat));
%%

simThreshs = linspace(.1,.98, 9);
TimeThreshs = fliplr(linspace(30,300,10));

ExperimentTDOA = struct();
ExperimentMaxProd = struct();
ExperimentUnaided = struct();

% Pre allocate the output matrix
PrecisionMatSpatial = inf(length(simThreshs), length(TimeThreshs),...
    length(corrThresh));
RecallMatSpatial = PrecisionMatSpatial;
PrecisionMatTDOA = PrecisionMatSpatial;
RecallMatTDOA = PrecisionMatSpatial;
PrecisionMatBaseline =PrecisionMatSpatial;
RecallMatBaseline =PrecisionMatSpatial;
PrecisionMatAcEnc= PrecisionMatSpatial;
RecallMatAcEnc =PrecisionMatSpatial;
% Create the basic structure, this needs to be done for every correlation
% threshold

nmi_baseline = zeros(length(simThreshs), length(TimeThreshs));
nmi_tdoa = nmi_baseline;
nmi_spatial = nmi_baseline;

for ii = 1:length(simThreshs)
    simThresh = simThreshs(ii);
    exampSpatial.cutoff = simThresh;
    exampTDOA.cutoff = simThresh;
    exampBaseline.cutoff=.5;
    disp('next simThresh')
    
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
            exampBaseline.chains = updateChainsEncounterFirst(exampBaseline);
            exampBaseline.Cluster_id= updateClusterID(exampBaseline);
            
            % Spatial and baseline
            [ClusteredSpatial, Detector,...
                gpldatout_chan, nmi ] = ...
                PrecisionRecallGPL(exampSpatial, parent, hyd);
            Results.Spatial(ii,jj) = ClusteredSpatial;
            Results.SpatialClusters(ii,jj) = {exampSpatial.Cluster_id};
            Results.GPLTable(ii,jj) = {gpldatout_chan};
            nmi_spatial(ii,jj) =nmi;
            
            if ii == jj ==1
                Results.Baseline(1,1) = Detector;
            end
            
            % TDOA precision/recall
            [ClusteredTDOA, ~,~,nmi] = PrecisionRecallGPL(exampTDOA, parent,...
                hyd);
            Results.TDOA(ii,jj) = ClusteredTDOA;
            Results.TDOAClusters(ii,jj) = {exampTDOA.Cluster_id};
            nmi_tdoa(ii,jj) =nmi;
            
            
            % Baseline clustering precision recall
            [ClusteredBase, ~,~,nmi] = PrecisionRecallGPL(exampBaseline, parent,...
                hyd);
            Results.AcousticEncounters(ii,jj) = ClusteredBase;
            Results.AcousticEncountersClusters(ii,jj) = ...
                {exampBaseline.Cluster_id};
            nmi_baseline(ii,jj) =nmi;
            disp('next correlatinthresh')
        
    end

    disp('Finish Sim Thres Iter')
end

%% Create a nice plot

simIdx =2;
figure
hold on

for jj=1:length(TimeThreshs)
    timeidx =jj;
    
        
        
        subplot(1,3,1)
        hold on
        plot(Results.TDOA(simIdx,timeidx).Recall,...
            Results.TDOA(simIdx,timeidx).Precision, '.-',...
            'color',[230/255 159/255 0/255])
        plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
         xlabel('Recall'); ylabel('Precision');
        title('TDOA')
        
        subplot(1,3,2)
        hold on
        plot(Results.Spatial(simIdx,timeidx).Recall,...
            Results.Spatial(simIdx,timeidx).Precision, '.-',...
            'color',[86/255 180/255 233/255])
        plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
         xlabel('Recall'); ylabel('Precision');
        title('Ambiguity Surface Method')
        
        subplot(1,3,3)
        hold on
        plot(Results.AcousticEncounters(simIdx,timeidx).Recall,...
            Results.AcousticEncounters(simIdx,timeidx).Precision, '.-',...
            'color',[0/255 158/255 115/255])
        plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        xlabel('Recall'); ylabel('Precision');
        title('Acoustic Encounters')
        
        
    
end
%%
close all
jitterAmount = 0.000;

for jj=1:length(TimeThreshs)
    timeidx =jj;
    figure;
    for ii =1:length(simThreshs)
        
        titlestr =[num2str(TimeThreshs(timeidx),'%0.2f'), ' sec TimeThresh ',...
            num2str(num2str(simThreshs(ii),'%0.2f')), ' SimThresh'];
        
        subplot(3,3,ii)
        
        hold on
        
        % Jitter the baseline detector so it shows up
        jitterValuesX = 2*(rand(size(Results.Baseline(1,1)))-0.5)...
            *jitterAmount;   % +/-jitterAmount max
        jitterValuesY = 2*(rand(size(Results.Baseline(1,1).Precision))-0.5)...
            *jitterAmount;   % +/-jitterAmount max
        
        
        plot(Results.TDOA(ii,timeidx).Recall,...
            Results.TDOA(ii,timeidx).Precision, '.-')
        plot(Results.Spatial(ii,timeidx).Recall,...
            Results.Spatial(ii,timeidx).Precision, '.-')
        plot(Results.AcousticEncounters(ii,timeidx).Recall,...
            Results.AcousticEncounters(ii,timeidx).Precision, '.-')
        
        plot(Results.Baseline(1,1).Recall+jitterValuesX,...
            Results.Baseline(1,1).Precision+jitterValuesY, '.-')
        xlabel('Recall'); ylabel('Precision');
        
        title(titlestr)
        
        
    end
end

%% Make plots showing cluster accuracy

 %legend('TDOA', 'Spatial', 'Acoustic Encounters Only','Detector Only')
 
 
figure
subplot(2,2,1); imagesc(TimeThreshs, simThreshs,nmi_baseline); colorbar
title('Acoustic Encounters')
xlabel('Elapsed Time Threshold (s)'); ylabel('Similarity Threshold');
%caxis([.2 .27])

subplot(2,2,2); imagesc(TimeThreshs, simThreshs,nmi_tdoa); colorbar
title('TDOA')
xlabel('Elapsed Time Threshold (s)'); ylabel('Similarity Threshold');
%caxis([.2 .27])

subplot(2,2,3); imagesc(TimeThreshs, simThreshs,nmi_spatial); colorbar
title('Spatial')
xlabel('Elapsed Time Threshold (s)'); ylabel('Similarity Threshold');
%caxis([.2 .27]) 
 
 
%%

simIndex =1;
figure; 
jitterAmount = 0.000;
for ii =1:3
    
    titlestr =[num2str(TimeThreshs(timeidx)), ' sec TimeThresh ',...
        num2str(num2str(simThreshs(ii),'%0.2f')), ' SimThresh'];
    
    subplot(3,4,ii)
    
    hold on
    
    % Jitter the baseline detector so it shows up
    jitterValuesX = 2*(rand(size(Results.Baseline(1,1)))-0.5)...
        *jitterAmount;   % +/-jitterAmount max
    jitterValuesY = 2*(rand(size(Results.Baseline(1,1).Precision))-0.5)...
    *jitterAmount;   % +/-jitterAmount max

    
    plot(Results.TDOA(ii,timeidx).Recall,...
        Results.TDOA(ii,timeidx).Precision, '.-')
    plot(Results.Spatial(ii,timeidx).Recall,...
        Results.Spatial(ii,timeidx).Precision, '.-')
    plot(Results.AcousticEncounters(ii,timeidx).Recall,...
        Results.AcousticEncounters(ii,timeidx).Precision, '.-')
    
    plot(Results.Baseline(1,1).Recall+jitterValuesX,...
        Results.Baseline(1,1).Precision+jitterValuesY, '.-')
    xlabel('Recall'); ylabel('Precision');
  
    title(titlestr)


end

 legend('TDOA', 'Spatial', 'Acoustic Encounters Only','Detector Only')

%% Area under the precision recall curve
maxA1 =0; maxA2 =0; maxA3=0;
figure
simIdx = 1
for ii=1:length(TimeThreshs)

    subplot(1,3,1)
    hold on
    area =Results.TDOA(simIdx,ii).Recall.*...
        Results.TDOA(simIdx,ii ).Precision;
     maxA1 = max([area, maxA1]);
     plot(area, '.-')
    title('TDOA')
    ylabel('Area under the Precision/Recall Curve')
    plot(...
        Results.Baseline.Recall.*Results.Baseline.Precision,...
        'k-')
    subplot(1,3,2)
    hold on
    
    area =Results.Spatial(simIdx,ii).Recall.*...
        Results.Spatial(simIdx,ii ).Precision;
    maxA2 = max([area, maxA2]);
    plot(area, '.-')
        plot(...
        Results.Baseline.Recall.*Results.Baseline.Precision,...
        'k-')
    %ylim([.4 .8])
    
    title('Spatial')
    
    subplot(1,3,3)
    hold on
    
    area =Results.AcousticEncounters(simIdx,ii).Recall.*...
        Results.AcousticEncounters(simIdx,ii ).Precision;
     maxA3 = max([area, maxA3]);
     plot(area, '.-')
        plot(...
        Results.Baseline.Recall.*Results.Baseline.Precision,...
        'k-')
    %ylim([.4 .8])
    
    title('Acoustic Encounters')
    
    maxOut =[maxA1 maxA2 maxA3];

end


