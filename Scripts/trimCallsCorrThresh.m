function callsTrimmed = trimCallsCorrThresh(hyd, corrThresh)
% Trim the hydrophone structure based on the correlation threshold for the
% template 


callsTrimmed = hyd;

% trim the scores

% Get the id's of the calls that were detected by two or mor hydrophone and
% with cross correlation scores above
for ii =1:length(callsTrimmed)
    
    if ~isempty(callsTrimmed(ii).detection.calls)
    scores = cat(1,callsTrimmed(ii).detection.calls.calltype_1_score);
    k2 = find(scores>corrThresh);
     
    % Trim calls
    callsTrimmed(ii).detection.calls = callsTrimmed(ii).detection.calls(k2);
    end
    
     
end

end