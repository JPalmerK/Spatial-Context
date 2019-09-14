function truthTable =createTruthTable(obj)
% This function is the classifier function that could use a bit
% of love. It's a bit hacky. Looks at calls and assigns a
% random probability score indicating that each call came from
% species 1.



if ~isfield(obj, 'truthtTable') ||  isempty(obj.Cluster_id)

    % Stick the real id's next to the predicted clusters
    truthTable =gather(obj.arrivalArray(:,end));
  
    
else
    truthTable = obj.truthtTable;
    
end

truthTable= array2table(truthTable, 'VariableNames',...
    {'TrueClust'});
truthTable.TrueSpp = zeros([height(truthTable), 1]);

agents = unique(truthTable.TrueClust);

%Pick one agent to be the target and the rest to be the false positives
spp = randsample(agents,1);
truthTable.TrueSpp(truthTable.TrueClust==spp) =1;

    
end










