function [R1 R2] = resistor_ratio(ratio)
% [R1 R2] = resistor_ratio(ratio)
%
% Return R1 and R2 in 1-percent standard resistor values,
% where
%
%  Vout = Vref * (R1 + R2)/R2
%  
%  R2/R1 = Vref/(Vout - Vref) = ratio
%
% and ratio is between 0 and 1.
% 
%

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

% NxN table of two-decades of resistor ratios
rtable2 = [rtable/10 rtable];
ratios = rtable2'*(1./rtable2);

% Find the closest match in the table
[i1 j1] = min(abs(ratios-ratio));
[i2 j2] = min(i1);

R2 = rtable2(j1(j2));
R1 = rtable2(j2);

