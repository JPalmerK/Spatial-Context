function rounded = ceilDecimal(val, decimal)

% Function for rounding up to the nerest decimal
% input - val: value to be rounded
% decimal - integer number of places to round to


rounded=ceil(10.^decimal*val)/10^decimal;


end