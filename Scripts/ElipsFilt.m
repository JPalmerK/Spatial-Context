function AS_propagated =  ElipsFilt(simStruct,averageLklhd_space, time_gaps,...
    grid_v,grid_h, idx)
% Function for creating the time projected similarity matricies using
% eliptical filter (EM Nosal)

array_struct = simStruct.array_struct;

% Create the filter for the master and the reference channel for the first
delays = simStruct.TDOA_vals(idx, :);
hyd_det = [simStruct.array_struct.master...
    simStruct.array_struct.slave(simStruct.child_idx(~isnan(delays)))];
filt = prod(simStruct.filtGrid(:,:,hyd_det),3);

% Grid spacing in m
grid_dy = grid_v/ (length(array_struct.latgrid)-1);
grid_dx = grid_h/ (length(array_struct.longrid)-1);

% Maximum swim speed
s = simStruct.s;

% Maximum swim distance
msd = gather(s*time_gaps);
th = 0:0.01:2*pi;

% Prealocate propagated filter size

AS_propagated = ones([...
    length(simStruct.array_struct.latgrid),...
    length(simStruct.array_struct.longrid),...
    length(time_gaps)]);


averageLklhd_space = gather(averageLklhd_space);

% Trim the average likelhiood space
AS_propagated(:,:,1) = (averageLklhd_space);

tic
for ii=1:length(msd)
    % Create the eliptical swim filter
    swim_filter_x = 0:grid_dx:msd(ii)+grid_dx;
    swim_filter_x = [-fliplr(swim_filter_x(2:end)) swim_filter_x];
    swim_filter_y = 0:grid_dy:msd(ii)+grid_dy;
    swim_filter_y = [-fliplr(swim_filter_y(2:end)) swim_filter_y];
    
    % Swim distance from cernter of filter to each grid point
    [Fx,Fy] = meshgrid(swim_filter_x, swim_filter_y);
    SD = sqrt(Fx.^2 + Fy.^2);
    
    
    % Find point in the filter with swim distance less than or equal to what is
    % possible and set them to 1
    % Note: This doesn't account for the probability of an animal swimming the
    % various distances. We might want to include this (i.e. make this weighted
    % toward the middle?)
    F = 0*SD;
    F((SD <= msd(ii))) = 1;

    
    AS_propagated(:,:,ii) = filt.*imdilate(averageLklhd_space, F);  
    
end
toc 
%Elapsed time is 3.750770 seconds.

% hyd_loc = cat(1,simStruct.hydrophone_struct.location);
% hyd_loc(4,:) =[];
% % Figure for paper
% figure
% subplot(1,2,1)
% imagesc(array_struct.longrid, array_struct.latgrid, squeeze(AS_propagated(:,:,1))), axis xy
% c = colorbar;
% c.Label.String = 'Likelihood';
% hold on
% scatter(hyd_loc(:,2), hyd_loc(:,1), 'filled', 'k', 'd')
% scatter(hyd_loc([1:3],2), hyd_loc([1:3],1), 'filled', 'w', 'd')
% scatter(hyd_loc(4,2), hyd_loc(4,1), 'filled', 'r', 'd')
% xlabel('Longitude')
% ylabel('Latitude')
% title('Projected Ambiguity Surface')
%
%
% subplot(1,2,2)
% imagesc(array_struct.longrid, array_struct.latgrid, squeeze(AS_propagated(:,:,end))), axis xy
% c = colorbar;
% c.Label.String = 'Likelihood';
% hold on
% scatter(hyd_loc(:,2), hyd_loc(:,1), 'filled', 'k', 'd')
% scatter(hyd_loc([1:3],2), hyd_loc([1:3],1), 'filled', 'w', 'd')
% scatter(hyd_loc(4,2), hyd_loc(4,1), 'filled', 'r', 'd')
% xlabel('Longitude')
% ylabel('Latitude')
% title('Projected Ambiguity Surface')



end

