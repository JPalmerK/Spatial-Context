%% Draw the Agents
function drawAgents(obj)

% Update the arrival array
if isempty(obj.arrivalArray)
    UpdateArrArray(obj);
end

% Hydrophone locations
hyd_table = struct2table(obj.hydrophone_struct);

TimeColorVals = parula(obj.spaceWhale.param_sim.dur+2);
ColorVals = lines( max([length(obj.spaceWhale.agent), max(obj.Cluster_id)]));
child_hyds = obj.array_struct.slave(obj.child_idx);


figure;

subplot(3,1,2)
hold on
scatter(obj.arrivalArray(:,end-1),...
    obj.arrivalArray(:,end-2),[],...
    ColorVals(obj.arrivalArray(:,end),:), 'f')

scatter(hyd_table.location(:,2),...
    hyd_table.location(:,1), 80,...
    'k', 'filled', 'd')
scatter(...
    hyd_table.location([obj.array_struct.master,...
    obj.array_struct.slave(obj.child_idx)],2),...
    hyd_table.location([obj.array_struct.master,...
    obj.array_struct.slave(obj.child_idx)],1),...
    'r', 'filled', 'd')

title('TOA on Parent')

subplot(3,1,1)

hold on
scatter(obj.arrivalArray(:,end-1),...
    obj.arrivalArray(:,end-2),[],...
    [TimeColorVals(round(obj.arrivalArray(:,1)),:)],...
    'f')

scatter(hyd_table.location(:,2), hyd_table.location(:,1), 80, 'k', 'filled', 'd')
scatter(...
    hyd_table.location([obj.array_struct.master, obj.array_struct.slave(obj.child_idx)],2),...
    hyd_table.location([obj.array_struct.master, obj.array_struct.slave(obj.child_idx)],1),...
    'r', 'filled', 'd')

colorbar
title('True Cluster')

if ~isempty(obj.Cluster_id)
    
    subplot(3,1,3)
    hold on
    scatter(obj.arrivalArray(:,end-1), obj.arrivalArray(:,end-2),[],...
        ColorVals(obj.Cluster_id,:), 'f')
    
    scatter(hyd_table.location(:,2), hyd_table.location(:,1),...
        80, 'k', 'filled', 'd')
    
    scatter(...
        hyd_table.location([obj.array_struct.master, obj.array_struct.slave(obj.child_idx)],2),...
        hyd_table.location([obj.array_struct.master, obj.array_struct.slave(obj.child_idx)],1),...
        'r', 'filled', 'd')
    
    titlestr = [obj.titleStr, ' Expected Clusters'];
    title(titlestr)
    
end



end