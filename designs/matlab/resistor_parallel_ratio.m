function [R1 R2] = resistor_parallel_ratio(Rin)
% [R1 R2] = resistor_ratio(Rin)
%
% Return R1 and R2 in 1-percent standard resistor values,
% where the two resistors in parallel give the value of Rin
%
%  Rin = R1*R2/(R1+R2)
%
% Search algorithm:
%
%  1/Rin = 1/R1 + 1/R2
%
% 1. Start with R1 = R2 = 2*Rin
% 2. Set R1 to next lowest standard value.
% 3. Calculate R2 = 1/(1/Rin - 1/R1)
% 4. Find the nearest standard value
% 5. Repeat at 2
%    Terminate if exact match found, or R2 too large,
%    or N iterations
%

% Multiple input resistors
if (length(Rin) > 1)
	% Call this procedure on each resistor
	M = length(Rin);
	R1 = zeros(1,M);
	R2 = zeros(1,M);
	for m = 1:M,
		[R1(m) R2(m)] = resistor_parallel_ratio(Rin(m));
	end
	return
end

% Check for a single resistor solution
R1 = resistor_one_percent_value(Rin);
if (R1 == Rin)
    R2 = Inf;
    return;
end

% Standard resistance values for 10 to 100.
rtable = [
10.0 10.2 10.5 10.7 11.0 11.3 11.5 11.8 12.1 12.4 12.7 13.0 ...
13.3 13.7 14.0 14.3 14.7 15.0 15.4 15.8 16.2 16.5 16.9 17.4 ...
17.8 18.2 18.7 19.1 19.6 20.0 20.5 21.0 21.5 22.1 22.6 23.2 ...
23.7 24.3 24.9 25.5 26.1 26.7 27.4 28.0 28.7 29.4 30.1 30.9 ...
31.6 32.4 33.2 34.0 34.8 35.7 36.5 37.4 38.3 39.2 40.2 41.2 ...
42.2 43.2 44.2 45.3 46.4 47.5 48.7 49.9 51.1 52.3 53.6 54.9 ...
56.2 57.6 59.0 60.4 61.9 63.4 64.9 66.5 68.1 69.8 71.5 73.2 ...
75.0 76.8 78.7 80.6 82.5 84.5 86.6 88.7 90.9 93.1 95.3 97.6];

% Starting resistance
R0 = 2*Rin;

% Determine decade (relative to 10 to 100-Ohm)
R0_decade = floor(log10(R0)) - 1;
R0 = R0/10^R0_decade;

% Find the standard value below this
R0_index = max(find(rtable < R0));
R0 = rtable(R0_index);

% Initialize the loop
N = 100;
R = zeros(5, N);
R(1, 1) = R0_index;
R(2, 1) = R0_decade;
R(3, 1) = R0_index;
R(4, 1) = R0_decade;
R(5, 1) = R0/2;
rtable_size = length(rtable);
for n = 2:N,
	% Next value of R1
	R(1, n) = R(1, n-1) - 1;
	R(2, n) = R(2, n-1);
	% Decade wrap
	if (R(1, n) == 0) 
		R(1, n) = rtable_size;
		R(2, n) = R(2, n) - 1;
    end
	R1 = rtable(R(1, n))*10^R(2, n);

	% Calculate R2
	R2 = 1/(1/Rin - 1/R1);
	
	% Skip if (R2 < 0) or (R2 > 1M)
	if (R2 < 0) 
		R(3, n) = 1;
		R(4, n) = -10;
		break
	end
	if (R2 > 1e6) 
		R(3, n) = 1;
		R(4, n) = -10;
		break
	end
	
	R2 = resistor_one_percent_value(R2);
	R2_decade = floor(log10(R2)) - 1;
	R2_10 = R2/10^R2_decade;
	m = find(abs(rtable - R2_10) < 0.01);
	if (isempty(m)) 
		fprintf('Failed to find R2 = %.2f\n', R2_10);
		break;
	end
	R(3, n) = m;
	R(4, n) = R2_decade;

	% Calculate parallel value
	R3 = 1/(1/R1 + 1/R2);
	R(5, n) = R3;

	% Terminate on exact match
	if (R3 == Rin) 
		R1 = rtable(R(1, n))*10^R(2, n);
		R2 = rtable(R(3, n))*10^R(4, n);
		return
	end
end
if (n < N)
	R = R(:,1:n-1);
end

% Return the best match
[m n] = min(abs(R(5, :) - Rin));
R1 = rtable(R(1, n))*10^R(2, n);
R2 = rtable(R(3, n))*10^R(4, n);