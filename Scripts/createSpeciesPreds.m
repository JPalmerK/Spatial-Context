function truthTable =createSpeciesPreds(obj)
% This function is the classifier function that could use a bit
% of love. It's a bit hacky. Looks at calls and assigns a
% random probability score indicating that each call came from
% species 1.




if ~isfield(obj, 'truthTable') 
    truthTable =createTruthTable(obj);
    
else
    truthTable = obj.truthTable;
    
end

score = zeros(height(truthTable), 1);

truthTable.ArrivalTime = gather(obj.arrivalArray(:,1));
spp1_idx = find(truthTable.TrueSpp ==1);
score(spp1_idx) = betarnd(obj.betaParm1, obj.betaParm2, [length(spp1_idx), 1]);


spp0_idx = find(truthTable.TrueSpp ==0);
score(spp0_idx) = betarnd(7, 20, [length(spp0_idx), 1]);
truthTable.Score = score;






end
