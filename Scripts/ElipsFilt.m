function AS_propagated =  ElipsFilt(simStruct,averageLklhd_space, time_gaps,...
    grid_v,grid_h)
% Function for creating the time projected similarity matricies using
% eliptical filter (EM Nosal)

array_struct = simStruct.array_struct;

% Grid spacing in m
grid_dy = grid_v/ (length(array_struct.latgrid)-1);
grid_dx = grid_h/ (length(array_struct.longrid)-1);

% Maximum swim speed
s = simStruct.s;

% Maximum swim distance
msd = gather(s*time_gaps);
th = 0:0.01:2*pi;

% Prealocate propagated filter size
try
    AS_propagated = ones([...
    length(simStruct.array_struct.latgrid),...
    length(simStruct.array_struct.longrid),...
    length(time_gaps)], 'gpuArray');
    
    % image dialate only works with uint8 on GPU
    averageLklhd_space = uint8(averageLklhd_space*100);
catch
    AS_propagated = ones([...
    length(simStruct.array_struct.latgrid),...
    length(simStruct.array_struct.longrid),...
    length(time_gaps)]);
end




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
    F(find(SD <= msd(ii))) = 1;
    
    AS_propagated(:,:,ii) = imdilate(averageLklhd_space, F);
    
end

if existsOnGPU(AS_propagated);
    AS_propagated = double(AS_propagated./100);
end

end

