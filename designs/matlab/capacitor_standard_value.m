function Cout = capacitor_standard_value(Cin)
% Cout = capacitor_standard_value(Cin)
%
% Convert Cin to a standard capacitor value.
% Cin can be a vector or an array.
%
% Standard capacitance values.
%
% I got this table from a mailer. However reviewing
% PCC14-KIT-ND and PCC13-KIT-ND, it seems like only every
% second value is common, i.e, 10, 12, 15, etc. 
%
ctable_common_low = [
10 11 12 13 15 16 18 20 22 24 27 30 ...
33 36 39 43 47 51 56 62 68 75 82 91
];

% Skip every second entry
ctable_common_low = ctable_common_low(1:2:end);

ctable_common_high = [
10 15 22 33 47 68];
ctable = [
ctable_common_low*1e-13 ...
ctable_common_low*1e-12 ...
ctable_common_low*1e-11 ...
ctable_common_low*1e-10 ...
ctable_common_high*1e-9 ...
ctable_common_high*1e-8 ...
ctable_common_high*1e-7 ...
ctable_common_high*1e-6 ...
ctable_common_high*1e-5 ...
ctable_common_high*1e-4 ...
10000e-6
];

[M N] = size(Cin);
Cout = zeros(M, N);
for m = 1:M,
	for n = 1:N,

		% Find the closest match in the table
		[i j] = min(abs(ctable-Cin(m, n)));

		Cout(m, n) = ctable(j);
	end
end

