function Rout = resistor_standard_value(Rin, tol)
% Rout = resistor_standard_value(Rin, tol)
%
% Convert Rin to a 1-percent or 5-percent standard resistor value.
% Rin can be a vector or an array.
%

if (tol == 1)
	Rout = resistor_one_percent_value(Rin);
elseif (tol == 5)
	Rout = resistor_five_percent_value(Rin);
else
	error('Tolerance (tol) must be 1 or 5')
end