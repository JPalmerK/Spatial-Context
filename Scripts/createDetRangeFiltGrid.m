function filt_grid = createDetRangeFiltGrid(obj,hydrophone_struct)

% Create grids indicating the probability of detecting a call at each of
% the TDOA grid locations (currently just half normal, could be replaced by
% more complicated propagation modelling or simplified using the maximum
% detection range (examp.truncateKm)

truncateKm = obj.truncateKm; % not currently used
array_struct = obj.array_struct;

%%  two or more MARUs

% Vertical and horizontal  grid granularity in meters
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

% Pre allocate the filter grid
filt_grid = zeros([size(array_struct.toa_diff{2}), length(hydrophone_struct)]);
hyd_idx = zeros(length(hydrophone_struct),2);

% Define Half normal distance function
pd4 = @(beta,sig,y) 1-exp(-(y/sig).^-beta);

% Detection Probability bubbles
for ii =1:length(hydrophone_struct)
    [~, lat_dy] = min(abs(array_struct.latgrid-hydrophone_struct(ii).location(1)));
    [~, lon_dx] = min(abs(array_struct.longrid-hydrophone_struct(ii).location(2)));
    
    swim_filter_x = (grid_dx*lon_dx) - (0:grid_dx:(length(array_struct.longrid)-1)*grid_dx);
    swim_filter_y = (grid_dy*lat_dy) - (0:grid_dy:(length(array_struct.latgrid)-1)*grid_dy);
    
    [Fx,Fy] = meshgrid(swim_filter_x, swim_filter_y);
    SD = sqrt(Fx.^2 + Fy.^2);
%     
%     
%     F = 0*SD;
%     F((SD <= truncateKm*1000)) = 1;
    
    pdfHyd =  pd4(2,18000,SD);
    pdfHyd(pdfHyd<0.001)=0;
    
    filt_grid(:,:,ii) = pdfHyd;
    hyd_idx(ii,1) = lat_dy;
    hyd_idx(ii,2) = lon_dx;
    
end


disp('blarg')
end