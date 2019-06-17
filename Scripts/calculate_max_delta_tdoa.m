function max_tdoa = calculate_max_delta_tdoa(hydrophone_struct, array_struct,...
    c, max_speed) 
% Function for calculating the maximum allowabale change in TDOA between 
% the hydrophone (master) and slave assuming a maximum swimming speed 
 
% Input:  
% c -  sound speed (m/s)
% # note currently only accepts single sound speed 
% array_struct - created by localize_parameter_input  
% hydrophone_struct - structure containing hydrophone information 
% max_speed - maximum animal speed in m/s 
 
% Returns: 
 
 
% check if max speed given assume 4.4m/s if not 
if nargin<5 
    max_speed=16; 
end 
 
 
% Loction of the focal (master) hydrophone 
focal_loc=hydrophone_struct(array_struct.master).location; 
focal_depth = hydrophone_struct(array_struct.master).depth; 
 
% Calculate the maximum TDOA between the master and each slave hydroophone 
for ii=1:length(array_struct.slave) 
 
        % slave/child hydrophone location and depth 
        slave_loc = hydrophone_struct(array_struct.slave(ii)).location; 
        slave_depth = hydrophone_struct(array_struct.slave(ii)).depth; 
         
        % change in depth between two hydrophones 
        depth_dist = abs(slave_depth-focal_depth); 
         
        % horizontal distance between two hydrophones 
        h_distance=vdist(focal_loc(1),focal_loc(2),... 
            slave_loc(1),slave_loc(2)); 
         
        % total distance between master and slave hydrophone 
        ab = sqrt(depth_dist^2+ h_distance^2); 
 
        % maximum allowable change in TDOA allowed based on swim speed 
        max_tdoa(ii) =(1/c)* (sqrt(ab^2+max_speed^2)-ab); 
 
 
end 
 
 
 
 
 
 
 
 
end