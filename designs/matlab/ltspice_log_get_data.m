function y = ltspice_log_get_data(filename, key, len)
%
% Read the stepped data from the LTSpice log file
%

fd = fopen(filename, 'r');

% Find the line with the key
l = fgetl(fd);
str = ['Measurement: ' key];
while (strcmp(l,str) ~= 1)
	l = fgetl(fd);
end

% Skip the next line
l = fgetl(fd);

% Get the measurement lines
for n = 1:len
	l = fgetl(fd);
	m = sscanf(l,'%f');
	y(n) = m(2);
end
fclose(fd);
