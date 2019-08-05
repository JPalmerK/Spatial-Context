% Run experiments
close all; clear all; clc

% Load the workspace (not KJP)
load('NewCompMethod.mat')

% Load the workspace(KJP)
whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2});

% Expected tdoa surfaces from hydrophone pairs 5,1 and 5,2
sigma = 0.0004^2 +.1^2 + .3^2; % 10 m deployment uncertainty, etc. etc. etc.
grid = cell2mat( array_struct.toa_diff(2));
grid2 =  cell2mat( array_struct.toa_diff(3));


% Get distance in meters between the lower and upper right
grid_v = vdist(min(array_struct.latgrid),...
    min(array_struct.longrid),...
    max(array_struct.latgrid),...
    min(array_struct.longrid));


% Get distance in meters between the lower left and lower right
grid_h = vdist(min(array_struct.latgrid),...
    min(array_struct.longrid),...
    min(array_struct.latgrid),...
    max(array_struct.longrid));

% Grid spacing in meters
deltalat_space = grid_v/ (length(array_struct.latgrid)-1);
deltalon_space = grid_h/ (length(array_struct.longrid)-1);

% How many grid squares per second can the whale move
s =8; % maximum swim speed of the animal
lat_persec = s / deltalat_space; % Distance (m) travelled in 1 sec
lon_persec = s / deltalon_space;


%% Example 1 two well localized calls
tdoaObs = [0 3]; % call 1 tdoa observations
gridAdj1 = grid-tdoaObs(1);
gridAdj2 = grid2-tdoaObs(2);

gridCall1 = cat(3,gridAdj1, gridAdj2);
gridCall1 =  normpdf(gridCall1, 0, 2*sigma)./normpdf(0,0,2*sigma);

AmbCall1 = nanmean(gridCall1,3);
normpdf(0,0,2*sigma)

% Delta TDOA, simulate second agent moving
deltatdoa1 = -1:.1:1;
score =[]; % similarity scores
dists=[]; % distance travelled

for ii =1:length(deltatdoa1)
    % second call, shifted a bit
    tdoaObs =[0+deltatdoa1(ii), 3+deltatdoa1(ii)];
    gridAdj1 = grid-tdoaObs(1);
    gridAdj2 = grid2-tdoaObs(2);
    gridCall1 = cat(3,gridAdj1, gridAdj2);
    gridCall1 =  normpdf(gridCall1, 0, 2*sigma)./normpdf(0,0,2*sigma);
    
    % Ambiguity surface call 2
    AmbCall2 = nanmean(gridCall1,3);
    
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    
    figure; subplot(1,2,1);imagesc(AmbCall1);subplot(1,2,2);imagesc(AmbCall2);
    score =[score, corrScore];
    dists=[dists, dist];
end
figure;
plot(deltatdoa1, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, Two Hyd- Moving Animal')

%% Example 1- Time projection only
close all
timegaps = linspace(0,100,10);
tdoaObs = [0 3];
gridAdj1 = grid-tdoaObs(1);
gridAdj2 = grid2-tdoaObs(2);
gridCall1 = cat(3,gridAdj1, gridAdj2);

% Ambiguity Surface for call 1 (projected in time)
gridCall1 =  normpdf(gridCall1, 0, 2*sigma)./normpdf(0,0,2*sigma);
AmbCall1 = nanmean(gridCall1,3);

% Ambiguity surface for call two that will not be projected in time
AmbCall2 = AmbCall1;

% Size of the filters 
filt_size_lat = ceil(lat_persec * timegaps);
filt_size_lon = ceil(lon_persec * timegaps);
filt_size = ([filt_size_lat; filt_size_lon])';
score=[];

for jj= 1:length(timegaps)
    
    
    % Grow the likelihood space based using image
    % processing max filter. Set the filter size based
    % on the maximum swim speed
    
    AmbCallProj = imdilate(AmbCall1, true(filt_size(jj,:)));
    % Get the comparison value of the projected prob
    % loc space and the next call in the squence
    
    
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    
        score =[score, corrScore];
        dists=[dists, dist];
    
%     Second MR method
%     [corrScore]= crossCorrSimScoresFilt(AmbCall2, AmbCallProj, ...
%         filt_size(jj,1), filt_size(jj,2));
%     
%     score =[score, corrScore];
%     dists=[dists, dist];
    
end


figure;
plot(timegaps, [score])
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, Two Hyd- Time Proj')



%% Example 2- one well localized call one not

tdoaObs = [0 3]; % call 1 tdoa observations
gridAdj1 = grid-tdoaObs(1);
gridAdj2 = grid2-tdoaObs(2);

gridCall1 = cat(3,gridAdj1, gridAdj2);
gridCall1 =  normpdf(gridCall1, 0, 2*sigma)./normpdf(0,0,2*sigma);

AmbCall1 = nanmean(gridCall1,3);
normpdf(0,0,2*sigma)

% Delta TDOA, simulate second agent moving
deltatdoa1 = -1:.1:1;
score =[]; % similarity scores
dists=[]; % distance travelled

for ii =1:length(deltatdoa1)
    
    % second call, shifted a bit
    tdoaObs =[0+deltatdoa1(ii)];
    
    % Call 2 grid 1
    AmbCall2 =  normpdf(grid-tdoaObs, 0, 2*sigma)./normpdf(0,0,2*sigma);
    
    
    figure; subplot(1,2,1);imagesc(AmbCall1);subplot(1,2,2);imagesc(AmbCall2);
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    
    score = [score, corrScore];
    dists=[dists, dist];
end

figure; plot(deltatdoa1, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, One Hyd')

%% Example 2- Time projection only
close all
timegaps = linspace(0,100,10);
tdoaObs = [0 3];
gridAdj1 = grid-tdoaObs(1);
gridAdj2 = grid2-tdoaObs(2);
gridCall1 = cat(3,gridAdj1, gridAdj2);

% Ambiguity Surface for call 1 (projected in time)
gridCall1 =  normpdf(gridCall1, 0, 2*sigma)./normpdf(0,0,2*sigma);
AmbCall1 = nanmean(gridCall1,3);

% Ambiguity surface for call two that will not be projected in time
AmbCall2 = squeeze(gridCall1(:,:,1));

% Size of the filters 
filt_size_lat = ceil(lat_persec * timegaps);
filt_size_lon = ceil(lon_persec * timegaps);
filt_size = ([filt_size_lat; filt_size_lon])';
score=[];

for jj= 1:length(timegaps)
    
   
    % Grow the likelihood space based using image
    % processing max filter. Set the filter size based
    % on the maximum swim speed
    
    AmbCallProj = imdilate(AmbCall1, true(filt_size(jj,:)));
    % Get the comparison value of the projected prob
    % loc space and the next call in the squence
    
  
    figure; subplot(1,2,1);imagesc(AmbCallProj);subplot(1,2,2);imagesc(AmbCall2);  
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    
        score =[score, corrScore];
        dists=[dists, dist];
    
%     Second MR method
%     [corrScore]= crossCorrSimScoresFilt(AmbCall2, AmbCallProj, ...
%         filt_size(jj,1), filt_size(jj,2));
%     
%     score =[score, corrScore];
%     dists=[dists, dist];
    
end


figure;
plot(timegaps, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, Two Hyd')


%% Example 3- Calls detected on only two instruments
close all

tdoaObs = [0 3]; % call 1 tdoa observations
gridAdj1 = grid-tdoaObs(1);
AmbCall1 =  normpdf(gridAdj1, 0, 2*sigma)./normpdf(0,0,2*sigma);

% Delta TDOA, simulate second agent moving
deltatdoa1 = linspace(-4,4,10);
score =[]; % similarity scores
dists=[]; % distance travelled

for ii =1:length(deltatdoa1)
    % second call, shifted a bit (moving animal)
    tdoaObs =[3+deltatdoa1(ii)];
    AmbCall2 = normpdf(grid2-tdoaObs, 0, 2*sigma)./normpdf(0,0,2*sigma);

    figure; subplot(1,2,1);imagesc(AmbCall1);subplot(1,2,2);imagesc(AmbCall2);
    [dist, corrScore]= crossCorrSimScores(AmbCall2,AmbCall1,...
        deltalat_space, deltalon_space);
    score = [score, corrScore];
    dists=[dists, dist];
end

figure
plot(deltatdoa1, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('One Hyd, One Hyd')

%% Example 3- Time projection only
close all
timegaps = linspace(0,100,10);
tdoaObs = [0 3];
gridAdj1 = grid-tdoaObs(1);

% Ambiguity Surface for call 1 (projected in time)
AmbCall1 =  normpdf(gridAdj1, 0, 2*sigma)./normpdf(0,0,2*sigma);


% Ambiguity surface for call two that will not be projected in time
AmbCall2 =  normpdf(grid2-tdoaObs(2), 0, 2*sigma)./normpdf(0,0,2*sigma);

% Size of the filters 
filt_size_lat = ceil(lat_persec * timegaps);
filt_size_lon = ceil(lon_persec * timegaps);
filt_size = ([filt_size_lat; filt_size_lon])';
score=[];

for jj= 1:length(timegaps)
    
   
    % Grow the likelihood space based using image
    % processing max filter. Set the filter size based
    % on the maximum swim speed
    
    AmbCallProj = imdilate(AmbCall1, true(filt_size(jj,:)));
    % Get the comparison value of the projected prob
    % loc space and the next call in the squence
    
  
    figure; subplot(1,2,1);imagesc(AmbCallProj);subplot(1,2,2);imagesc(AmbCall2);  
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    
        score =[score, corrScore];
        dists=[dists, dist];
    
%     Second MR method
%     [corrScore]= crossCorrSimScoresFilt(AmbCall2, AmbCallProj, ...
%         filt_size(jj,1), filt_size(jj,2));
%     
%     score =[score, corrScore];
%     dists=[dists, dist];
    
end


figure;
plot(timegaps, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, Two Hyd')



