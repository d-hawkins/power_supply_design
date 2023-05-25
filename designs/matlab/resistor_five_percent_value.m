function Rout = resistor_five_percent_value(Rin)
% Rout = resistor_five_percent_value(Rin)
%
% Convert Rin to a 5-percent standard resistor value.
% Rin can be a vector or an array.
%

% Standard resistance values as per Digikey
% 0603 1/16W resistor kits CR5A-KIT, CR6A-KIT.
%
% These kits have the values in;
%   x1 , x10, x100, x1000, x10000
% 
% and one value in x100000, i.e., 1M
%
rtable = [
10 11 12 13 15 16 18 20 22 24 27 30 ...
33 36 39 43 57 56 62 68 75 82 91
];

[M N] = size(Rin);
Rout = zeros(M, N);
for m = 1:M,
	for n = 1:N,

		% Determine decade
		if (Rin(m, n) ~= 0)
			Rdecade = floor(log10(Rin(m, n)));
		else
			Rdecade = 0;
		end

		% Convert to 10 to 100 decade
		Rin10 = Rin(m, n)/10^Rdecade*10;

		% Find the closest match in the table
		[i j] = min(abs(rtable-Rin10));

		if (Rin(m, n) < 10)
			Rout(m, n) = 10;		
		elseif (Rin(m, n) >= 1e6)
			Rout(m, n) = 1e6;		
		else
			Rout(m, n) = rtable(j)/10*10^Rdecade;
		end
	end
end

