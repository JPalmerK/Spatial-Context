function TDOA_vals = UpdateTDOA(simStruct)
% Function for extracting/calculating TDOA values from arrival
% times matrix. In progress- allow different array configurations



TDOA_vals =zeros(size(simStruct.arrivalArray,1), length(simStruct.child_idx));

% Time difference of arrivals (can only handle two atm)

for jj =1:length(simStruct.child_idx)
    TDOA_vals(:,jj) = simStruct.arrivalArray(:, 1)-simStruct.arrivalArray(:, jj+1);
end

end
