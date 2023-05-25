function [Lout Rout] = inductor_standard_value(Lin,family)
%
% 9/26/2012 D. W. Hawkins (dwh@ovro.caltech.edu)
%
% Convert Lin to a standard inductor value, where 
% 'family' indicates a manufacturer family identifier
% used in this m-file. Read the source for the available
% options.
%
% The function returns both the inductor value and its
% maximum DC resistance.
%
% eg.,
%
% [Lout Rdcr] = inductor_standard_value(4.6e-6, 'IHLP-2525CZ-01');
%
% will return Lout = 4.7uF, Rdcr = 40mOhm
%

if (nargin ~= 2)
	family = 'IHLP-2525CZ-01'
end

% Vishay IHLP-2525CZ-01 series (table copied from lp25cz01.pdf)
% uH and mOhm
vishay_ihlp2525_cz_01 = [
0.01  1.7
0.15  2.5
0.20  3.0
0.22  2.8
0.33  3.9
0.47  4.2
0.68  5.5
0.82  8
1.0   10
1.5   15
2.2   20
3.3   30
4.7   40
6.8   60
8.2   68
10   105];

% Vishay IHLP-3232DZ-01 series (table copied from lp32dz01.pdf)
% uF and mOhm
vishay_ihlp3232_dz_01 = [
 0.22  1.68
 0.33  2.14
 0.47  2.62
 0.68  3.67
 0.82  4.42
 1.0   5.78
 2.2  13.7
 3.3  17.7
 4.7  32.0
 5.6  35.5
 6.8  47.7
 8.2  50.8
10.0  59.9];

% Vishay IHLP-3232DZ-11 series (lower DCR and temperature rise)
% uF and mOhm
vishay_ihlp3232_dz_11 = [
 0.22  1.35
 0.33  2.15
 0.47  2.38
 0.68  3.22
 0.82  3.88
 1.0   4.63
 2.2   9.41
 3.3  14.9
 4.7  22.6
 5.6  28.6
 6.8  33.4
 8.2  45.0
10.0  51.8
15.0  65.3
22.0  94.2
33.0  144];

if (strcmp(family,'IHLP-3232DZ-01') == 1)
	L_table = vishay_ihlp3232_dz_01(:,1)*1e-6;
	R_table = vishay_ihlp3232_dz_01(:,2)*1e-3;
elseif (strcmp(family,'IHLP-3232DZ-11') == 1)
	L_table = vishay_ihlp3232_dz_11(:,1)*1e-6;
	R_table = vishay_ihlp3232_dz_11(:,2)*1e-3;
else
	% IHLP-2525CZ-01
	L_table = vishay_ihlp2525_cz_01(:,1)*1e-6;
	R_table = vishay_ihlp2525_cz_01(:,2)*1e-3;
end

% Convert to standard value
[M N] = size(Lin);
Lout = zeros(M, N);
Rout = zeros(M, N);
for m = 1:M,
	for n = 1:N,

		% Find the closest match in the table
		[i j] = min(abs(L_table-Lin(m, n)));

		Lout(m, n) = L_table(j);
		Rout(m, n) = R_table(j);
	end
end

