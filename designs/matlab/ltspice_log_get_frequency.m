function y = ltspice_log_get_frequency(filename)
%
% Read the stepped frequency data from the LTSpice log file
%
fd = fopen(filename, 'r');

% Parse the file for the lines starting with '.step freq='
n = 1;
while (feof(fd) ~= 1)
	% Read a line
	l = fgetl(fd);

	% Is it a frequency step
	if (strncmp(l,'.step freq=', 10) == 1)
		y(n) = sscanf(l,'.step freq=%f');
		n = n + 1;
	end
	% else skip it
end
fclose(fd);
