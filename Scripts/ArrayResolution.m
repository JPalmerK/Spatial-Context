% Run experiments
close all; clear all
clear classes; clc

whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2})
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
ssp = 1500
grid_depth = 5
child_idx = 1:8;
PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
s = 8;
truncateKm=12;

array_struct = array_struct_data(5).array;
clear array_struct_data;

%% For each pair of hydrophones, determine where calls can be detected by
%  two or more MARUs

grid_v = vdist(min(array_struct.latgrid),...
    min(array_struct.longrid),...
    max(array_struct.latgrid),...
    min(array_struct.longrid));


% Get distance in meters between the lower left and lower right
grid_h = vdist(min(array_struct.latgrid),...
    min(array_struct.longrid),...
    min(array_struct.latgrid),...
    max(array_struct.longrid));

% Grid spacing in m
grid_dy = grid_v/ (length(array_struct.latgrid)-1);
grid_dx = grid_h/ (length(array_struct.longrid)-1);



% Grid spacing in m
grid_dy = grid_v/ (length(array_struct.latgrid)-1);
grid_dx = grid_h/ (length(array_struct.longrid)-1);

filt_grid = zeros([size(array_struct.toa_diff{2}), 10]);
hyd_idx = zeros(10,2);

hyd_ids = [array_struct.master, array_struct.slave]

% Detection Probability bubbles
for ii =1:10
[~, lat_dy] = min(abs(array_struct.latgrid-hydrophone_struct(ii).location(1)));
[~, lon_dx] = min(abs(array_struct.longrid-hydrophone_struct(ii).location(2)));

swim_filter_x = (grid_dx*lon_dx) - (0:grid_dx:(length(array_struct.longrid)-1)*grid_dx);
swim_filter_y = (grid_dy*lat_dy) - (0:grid_dy:(length(array_struct.latgrid)-1)*grid_dy);

[Fx,Fy] = meshgrid(swim_filter_x, swim_filter_y);
SD = sqrt(Fx.^2 + Fy.^2);

F = 0*SD;
F((SD <= truncateKm*1000)) = 1;

filt_grid(:,:,ii) = F;
hyd_idx(ii,1) = lat_dy;
hyd_idx(ii,2) = lon_dx;

end

clear dist_mat


% Pick a master
master_id = 5;



% Loop through every slave if there is more than one are where the
% probability surfaces overlap

tdoa_pairs = struct();
tdoa_mat = zeros((length(array_struct.longrid)*length(array_struct.latgrid)), 10)/0;


for ii =1:10

    
    % Index of where there are overlapping detection probability spaces
    idx = find(sum(filt_grid(:,:,unique([master_id,ii])),3)>1);
   
    [r c] = ind2sub(size(filt_grid(:,:,1)), idx);
   
    % Matrix for the distance between each call
    dist_mat = zeros(length(idx), length(idx));
   
    % get the distance between each location and the remaining locations
    at1=[];
    at2=[];
    tdoa = zeros(size(idx));
    
    % Depth Difference between the hydrophones
    depth_range = (hydrophone_struct(5).depth...
        - hydrophone_struct(ii).depth)^2 ;
    
    for jj =1:length(idx)
   
        % distance between the hydrophone and each point in the hydrophone
        % grid that falls within the overlapping area
        at1(1) =  grid_dy *(r(jj) - hyd_idx(5, 1));
        at2(1) =  grid_dy *(r(jj) - hyd_idx(ii, 1));
       
        at1(2) =   grid_dx *(c(jj) - hyd_idx(5, 2));
        at2(2) =   grid_dx *(c(jj) - hyd_idx(ii, 2));
        
        
        tot_dist(1) =  sqrt(sum([at1.^2 depth_range]));
        tot_dist(2) =  sqrt(sum([at2.^2 depth_range]));
        
        
        
        
        %TDOA for each point 
        tdoa(jj) = diff(tot_dist)/ssp;
         
       
    end
    tdoa_pairs(ii).Val = tdoa;
    tdoa_pairs(ii).IDX = idx;
    tdoa_pairs(ii).rc = [r, c];
    
    tdoa_mat(idx,ii) =tdoa;
   

end




aa = zeros(size(filt_grid(:,:,1)))/0;
aa(tdoa_pairs(1).IDX) = tdoa_pairs(1).Val;
bb = tdoa_pairs(1).IDX;
cntr_idx = [floor(quantile(tdoa_pairs(1).rc(:,1), .75)), floor(quantile(tdoa_pairs(1).rc(:,2), 0.25))];
centr_val = aa(cntr_idx(1), cntr_idx(2));
imagesc(aa)

dtdoa = tdoa_pairs(1).Val - centr_val;
aa(tdoa_pairs(1).IDX) = dtdoa;
imagesc(aa)


% Step through the TDOAs and calculate probability grids and tdoa vals for
% each location

dtdoa = tdoa_mat - centr_val;
mu = zeros(size(dtdoa));
sigma = (mu+1)*sqrt(PosUncertsigma+2);
x =dtdoa; %values
likelihood = normpdf(x,mu,sigma);

% Normalize

LikelihoodNormFac= normpdf(0,0,sigma);


P_tdoaOnly =likelihood./LikelihoodNormFac;
P_tdoaOnly(:,11) = nanmax(P_tdoaOnly,[],2);

aa = zeros(size(filt_grid(:,:,1)))/0
aa(unique(cat(1,tdoa_pairs(1:end).IDX))) =1;
imagesc(aa)



subplot(3,2,1) % Hydrophones 5 and 1
slave_idx = [1];
aa = zeros(size(filt_grid(:,:,1)))/0
val = reshape(aa,[],1);
val((sum(~isnan(P_tdoaOnly(:,slave_idx)),2)>0)) = 1; % set this for the area where TDOAs are available
val = sum([val, -min(P_tdoaOnly(:,slave_idx),[],2)],2);
val = reshape(val, size(filt_grid(:,:,1)));
imagesc(array_struct.longrid, array_struct.latgrid, 1-val), axis xy, colorbar
caxis(gca,[-.001 1 ]);
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), ones(1,10)*80, 'd', 'filled','k');
scatter(hyd_arr(5,2), hyd_arr(5,1), ones(1,1)*80, 'd', 'filled','r');
scatter(hyd_arr(slave_idx,2), hyd_arr(slave_idx,1), ones(1,1)*80, 'd', 'filled','b');
scatter(array_struct.longrid(cntr_idx(2)), array_struct.latgrid(cntr_idx(1)),ones(1,1)*80,[.5 .5 .5], 'v', 'filled')
title('Maximum Uncertainty 2 Channels')
colormap( [1 1 1;  flipud(parula(1024))])
xlabel('Longitude Deg')
ylabel('Latitude Deg')





subplot(3,2,2) % Hydrophones 5 and 1
slave_idx = [1 2];
aa = zeros(size(filt_grid(:,:,1)))/0
val = reshape(aa,[],1);
val((sum(~isnan(P_tdoaOnly(:,slave_idx)),2)>0)) = 1; % set this for the area where TDOAs are available
val = sum([val, -min(P_tdoaOnly(:,slave_idx),[],2)],2);
val = reshape(val, size(filt_grid(:,:,1)));
imagesc(array_struct.longrid, array_struct.latgrid, 1-val), axis xy, colorbar
caxis(gca,[-.001 1 ]);
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), ones(1,10)*80, 'd', 'filled','k');
scatter(hyd_arr(5,2), hyd_arr(5,1), ones(1,1)*80, 'd', 'filled','r');
scatter(hyd_arr(slave_idx,2), hyd_arr(slave_idx,1), ones(1,1)*80, 'd', 'filled','b');
scatter(array_struct.longrid(cntr_idx(2)), array_struct.latgrid(cntr_idx(1)),ones(1,1)*80,[.5 .5 .5], 'v', 'filled')
title('Maximum Uncertainty 3 Channels')
colormap( [1 1 1;  flipud(parula(1024))])
xlabel('Longitude Deg')
ylabel('Latitude Deg')




subplot(3,2,3) % Hydrophones 5 and 1
slave_idx = [1 2 3];
aa = zeros(size(filt_grid(:,:,1)))/0
val = reshape(aa,[],1);
val((sum(~isnan(P_tdoaOnly(:,slave_idx)),2)>0)) = 1; % set this for the area where TDOAs are available
val = sum([val, -min(P_tdoaOnly(:,slave_idx),[],2)],2);
val = reshape(val, size(filt_grid(:,:,1)));
imagesc(array_struct.longrid, array_struct.latgrid, 1-val), axis xy, colorbar
caxis(gca,[-.001 1 ]);
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), ones(1,10)*80, 'd', 'filled','k');
scatter(hyd_arr(5,2), hyd_arr(5,1), ones(1,1)*80, 'd', 'filled','r');
scatter(hyd_arr(slave_idx,2), hyd_arr(slave_idx,1), ones(1,1)*80, 'd', 'filled','b');
scatter(array_struct.longrid(cntr_idx(2)), array_struct.latgrid(cntr_idx(1)),ones(1,1)*80,[.5 .5 .5], 'v', 'filled')
title('Maximum Uncertainty 4 Channels')
colormap( [1 1 1;  flipud(parula(1024))])
xlabel('Longitude Deg')
ylabel('Latitude Deg')


subplot(3,2,4) % Hydrophones 5 and 1
slave_idx = [1 2 3 4];
aa = zeros(size(filt_grid(:,:,1)))/0
val = reshape(aa,[],1);
val((sum(~isnan(P_tdoaOnly(:,slave_idx)),2)>0)) = 1; % set this for the area where TDOAs are available
val = sum([val, -min(P_tdoaOnly(:,slave_idx),[],2)],2);
val = reshape(val, size(filt_grid(:,:,1)));
imagesc(array_struct.longrid, array_struct.latgrid, 1-val), axis xy, colorbar
caxis(gca,[-.001 1 ]);
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), ones(1,10)*80, 'd', 'filled','k');
scatter(hyd_arr(5,2), hyd_arr(5,1), ones(1,1)*80, 'd', 'filled','r');
scatter(hyd_arr(slave_idx,2), hyd_arr(slave_idx,1), ones(1,1)*80, 'd', 'filled','b');
scatter(array_struct.longrid(cntr_idx(2)), array_struct.latgrid(cntr_idx(1)),ones(1,1)*80,[.5 .5 .5], 'v', 'filled')
title('Maximum Uncertainty 5 Channels')
colormap( [1 1 1;  flipud(parula(1024))])
xlabel('Longitude Deg')
ylabel('Latitude Deg')


subplot(3,2,5) % Hydrophones 5 and 1
slave_idx = [1 2 3 4 6 ];
aa = zeros(size(filt_grid(:,:,1)))/0
val = reshape(aa,[],1);
val((sum(~isnan(P_tdoaOnly(:,slave_idx)),2)>0)) = 1; % set this for the area where TDOAs are available
val = sum([val, -min(P_tdoaOnly(:,slave_idx),[],2)],2);
val = reshape(val, size(filt_grid(:,:,1)));
imagesc(array_struct.longrid, array_struct.latgrid, 1-val), axis xy, colorbar
caxis(gca,[-.001 1 ]);
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), ones(1,10)*80, 'd', 'filled','k');
scatter(hyd_arr(5,2), hyd_arr(5,1), ones(1,1)*80, 'd', 'filled','r');
scatter(hyd_arr(slave_idx,2), hyd_arr(slave_idx,1), ones(1,1)*80, 'd', 'filled','b');
scatter(array_struct.longrid(cntr_idx(2)), array_struct.latgrid(cntr_idx(1)),ones(1,1)*80,[.5 .5 .5], 'v', 'filled')
title('Maximum Uncertainty 6 Channels')
colormap( [1 1 1;  flipud(parula(1024))])
xlabel('Longitude Deg')
ylabel('Latitude Deg')



subplot(3,2,6) % Hydrophones 5 and 1
slave_idx = [1 2 3 4 6 8];
aa = zeros(size(filt_grid(:,:,1)))/0
val = reshape(aa,[],1);
val((sum(~isnan(P_tdoaOnly(:,slave_idx)),2)>0)) = 1; % set this for the area where TDOAs are available
val = sum([val, -min(P_tdoaOnly(:,slave_idx),[],2)],2);
val = reshape(val, size(filt_grid(:,:,1)));
imagesc(array_struct.longrid, array_struct.latgrid, 1-val), axis xy, colorbar
caxis(gca,[-.001 1]);
hold on
scatter(hyd_arr(:,2), hyd_arr(:,1), ones(1,10)*80, 'd', 'filled','k');
scatter(hyd_arr(5,2), hyd_arr(5,1), ones(1,1)*80, 'd', 'filled','r');
scatter(hyd_arr(slave_idx,2), hyd_arr(slave_idx,1), ones(1,1)*80, 'd', 'filled','b');
scatter(array_struct.longrid(cntr_idx(2)), array_struct.latgrid(cntr_idx(1)),ones(1,1)*80,[.5 .5 .5], 'v', 'filled')
title('Maximum Uncertainty 7 Channels')
colormap( [1 1 1;  flipud(parula(4096))])
xlabel('Longitude Deg')
ylabel('Latitude Deg')


%% Do the same thing for the MAX prod simulation




cntr_idx = [floor(quantile(tdoa_pairs(1).rc(:,1), .75)), floor(quantile(tdoa_pairs(1).rc(:,2), 0.25))];
centr_val = aa(cntr_idx(1), cntr_idx(2));
tdoa_matMp = tdoa_mat;
tdoa_matMp(:,[4,5]) =[];

% Create max prod matrix
maxProdmat = zeros(size(tdoa_matMp))/0;

hyds = array_struct.slave

% Step through all the rows in the area matrix and calculate the maxprod
for ii=1:size(tdoa_matMp,1)
    
    row_val = tdoa_matMp(ii,:);
    
    if sum(~isnan(row_val))>0
        
        hyd_ridx = find(~isnan(row_val));
        
        for jj=1:length(hyd_ridx)
        
        
        grid = array_struct.toa_diff{hyd_ridx(jj)+1};
        gridAdj = grid - 
        
        
        end
        
        
        
        
        
        
    end
    
    
    
end













