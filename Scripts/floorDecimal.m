function rounded = floorDecimal(val, decimal)

% Function for rounding up to the nerest decimal
% input - val: value to be rounded
% decimal - integer number of places to round to


rounded=floor(10.^decimal*val)/10^decimal;


end