function [array_struct] = setup_TDOA_grid(hydrophone_struct, array_struct, localize_struct,parm);
    %setup_TDOA_grid create time difference of arrival grids for all master
    % slave pairs
    
    
    % Depth of calling whale and bool for continuous or variable sound  speed
    % profile
    
    grid_depth = localize_struct.parm.grid_depth;
    ssp_choice=localize_struct.parm.ssp_choice;
    
    if isfield(localize_struct.parm,'ssp') ==1
        ssp = localize_struct.parm.ssp;
    end
    
    if isfield(localize_struct.parm,'c_nominal') ==0
        c_nominal=1500;
    else
        c_nominal=localize_struct.parm.c_nominal;
    end
    
    
    % For each hydrophone create TDOA grid
    for(j=1:length(hydrophone_struct))
        
        % Hyerophone depth
        depth(j)=hydrophone_struct(j).depth;
    end
    
    
    % For each array (typically only 1) create for the sound speed
    % and distance between each MARU and lat/lon in the grid spacing
    for(k=1:length(array_struct))
        
        
        local_array=[array_struct(k).master,array_struct(k).slave];
        
        % Create a grid space the size of the lat/long grid
        d=zeros(length(array_struct(k).latgrid),length(array_struct(k).longrid),length(local_array));
        speed=d;
        
        % For each hydrophone pair in the local array
        for(i=1:length(local_array))
            
            % Hydrophone location
            array=hydrophone_struct(local_array(i)).location;
            
            % Distance between the calling animal and the seafloor at
            % hydrohone i
            d_distance = depth(local_array(i))-grid_depth;
            
            % Iterate through the entire lat/long grid space
            for(j=1:length(array_struct(k).latgrid))
                
                % Report the percentage of the aray done
                perccount(j,length(array_struct(k).latgrid));
                for(jj=1:length(array_struct(k).longrid));
                    
                    
                    % Get geographical references for each lat/lon
                    geo(1)=array_struct(k).latgrid(j);
                    geo(2)=array_struct(k).longrid(jj);
                    
                    % Calculate the distance between hydrophone (i)
                    % and the grid location
                    h_distance=vdist(geo(1),geo(2),array(1),array(2));
                    
                    % Total h&v distance between the lat/lon point
                    % and hydrophone (i)
                    d(j,jj,i) = sqrt(h_distance^2 + d_distance^2);
                    
                    % Calculate the sound sepeed
                    if(ssp_choice==2)
                        
                        % Mean soundspeed using SSP profiel
                        speed(j,jj,i) = aute(ssp(:,2),ssp(:,1),d_distance,h_distance);
                    else
                        speed(j,jj,i)=ssp(1,2);
                    end
                    
                end
            end
        end
        
        
        % Master hydrophone location and depth
        master=hydrophone_struct(local_array(1)).location;
        master_depth=depth(local_array(1));
        
        dt=parm.skip/parm.sample_freq;
        
        
        % For each hydrophone in the slave array, calculate the expected
        % tdoa and phase
        for(i=2:length(local_array))
            
            % TDOA- distance divideded by sepped
            array_struct(k).toa_diff{i}= d(:,:,i)./speed(:,:,i) - d(:,:,1)./speed(:,:,1);
            
            array=hydrophone_struct(local_array(i)).location;
            h_distance=vdist(master(1),master(2),array(1),array(2));
            
            array_depth = depth(local_array(i));
            slant_range = sqrt(h_distance^2 + (master_depth-array_depth)^2);
            
            travel_time=slant_range/c_nominal; % should be worst case, mean speed no
            % slower than this
            
            array_struct(k).phase(i-1)=ceil(travel_time/dt);
        end
        
    end
    
    
