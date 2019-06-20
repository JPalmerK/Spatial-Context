% Function for drawing the current simulation agents

function drawWhales(spaceWhale, hyd_arr)
figure



% Make color scale
ColorVals = lines(length(spaceWhale.agent));
TimeColorVals = parula(spaceWhale.param_sim.dur);

for ii=1:length(spaceWhale.agent)
    
    subplot(1,2,2)
    hold on 
    plot(spaceWhale.agent(ii).location(:,2),...
        spaceWhale.agent(ii).location(:,1),...
        'Color', [ColorVals(ii,:)])

    % Add points indicating when the agent produced a call
    call_times = spaceWhale.agent(ii).call_times;
    scatter(spaceWhale.agent(ii).location(call_times,2),...
        spaceWhale.agent(ii).location(call_times,1),40,...
       [ColorVals(ii,:)], 'filled')
   
   
   subplot(1,2,1)
    hold on
    % Add points indicating when the agent produced a call
    call_times = spaceWhale.agent(ii).call_times;
    scatter(spaceWhale.agent(ii).location(call_times,2),...
        spaceWhale.agent(ii).location(call_times,1),40,...
        [TimeColorVals(call_times,:)], 'filled')

   
   

end
subplot(1,2,1)
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')

subplot(1,2,2)
scatter(hyd_arr(:,2), hyd_arr(:,1), 80, 'k', 'filled', 'd')


% Create ylabel
ylabel({'Lat'},'FontWeight','bold');
% Create xlabel
xlabel({'Lon'},'FontWeight','bold');
% Create title
title({'Simulated Movement and Call Times'},'FontSize',14);
end