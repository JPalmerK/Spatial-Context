function Sim_mat = SimulationMatrix2D(time_gaps, latgrid, longrid, s);

% This function creates the simulation matrix using the low
% memory approach. This should be used in most cases where not
% exploring the algorithims in depth.


% Grid X/Y space
% Get distance in meters between the lower and upper right
grid_v = vdist(min(latgrid),...
    min(longrid),...
    max(latgrid),...
    min(longrid));


% Get distance in meters between the lower left and lower right
grid_h = vdist(min(latgrid),...
    min(longrid),...
    min(latgrid),...
    max(longrid));




% Grid X/Y space
deltalat_space = gpuArray(grid_v/ (length(latgrid)-1));
deltalon_space = gpuArray(grid_h/ (length(longrid)-1));

% How many grid squares per second can the whale move
lat_persec = s / deltalat_space;
lon_persec = s / deltalon_space;




parfor ii=1:length(time_gaps)
    tic
    aa =  time_gaps(ii).next;
    bb = time_gaps(ii).proj;
    
    cormap= gpuArray();
    sim_score = zeros(size(aa,3),1, 'gpuArray');
    
    for jj = 1:size(aa,3)
        cc =gather(squeeze(aa(:,:,jj)));
        dd =gather(squeeze(bb(:,:,jj)));
        
        cc(cc<.01)=0;
        dd(dd<.01)=0;
        corr2d = xcorr2(cc,dd);
        % Normalize
        cormap =cormap./sum(sum(gather(squeeze(bb(:,:,jj))).^2));
        cormap= cat(3,cormap, corr2d);
        
        
        [ssr,snd] = max(corr2d(:));
        [ij,ji] = ind2sub(size(corr2d),snd);
        
        shiftr = (size(cc,1)-(ij))*deltalat_space;
        shiftc = (size(cc,2)-(ji))*deltalon_space;
        
        dist = sqrt(shiftr^2+shiftc^2);
        speed = nansum(dist/time_gaps(ii).calldiffs(jj));
        speed=nansum(speed+.001)-.001;
        sim_score(jj)=ssr*(speed<s);
        
    end
    time_gaps(ii).cormap=cormap;
    toc
    
    
end
Sim_mat = time_gaps;


end