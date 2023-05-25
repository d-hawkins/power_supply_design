function LTC3600_bode_response(example)
% -----------------------------------------------------------------
% LTC3600_bode_response.m
%
% 9/20/2012 D. W. Hawkins (dwh@ovro.caltech.edu)
%
% LTC3600 Bode Response Plots.
%
% This m-file reads and plots the Bode response measured using
% the LTSpice LTC3600_bode.asc circuit. The circuit writes the
% results to its .log file. That log file is then renamed to
% something like LTC3600_bode1.txt so that LTspice does not
% overwrite the results. This m-file parses the .log file and
% extracts the LTspice measurements.
%
% -----------------------------------------------------------------
% Examples
% --------
%
% #1: 1A load, default compensation network R = 100k, C = 50pF.
%
% -----------------------------------------------------------------

%clear all

% Add the LTspice .log parsers to the path
% * assumes the script is called from within designs/LTC3600/
path('../matlab',path);

% -----------------------------------------------------------------
% LTspice log file
% -----------------------------------------------------------------
%
% Select example 
if (nargin == 0)
	% Default to example 1
	example = 1; 
end

if ~isnumeric(example)
	error('Error: please provide an example number')
end

fprintf('Loading example#%d data\n', example)
if (example == 1) 
	% LTC3600_bode1.txt contains a coarse frequency sweep
	% for R = 100kOhm, C = 50pF
	%
	filename = 'LTC3600_bode1.txt';
elseif (example == 2) 
	% LTC3600_bode2.txt contains a fine frequency sweep
	% (excluding 100Hz and 300Hz)
	% for R = 15kOhm, C = 330pF
	%
	filename = 'LTC3600_bode2.txt';
elseif (example == 3) 
	% LTC3600_bode3.txt contains a fine frequency sweep
	% (excluding 100Hz and 300Hz) for a 5V output supply
	% (see the .txt file comments)
	%
	filename = 'LTC3600_bode3.txt';
else
	error('Invalid example number')
end

% -----------------------------------------------------------------
% LTC3600 Parameters
% -----------------------------------------------------------------
%
% These parameters match those used in the LTspice simulation
%
% Switching frequency
fsw  = 1e6;

if ((example == 1) || (example == 2))
	% Feedback reference and output voltage
	Vref = 2.5;
	Vout = 2.5;

	% Inductor
	L    = 3.3e-6;
	Rdcr = 30e-3;

	% Output capacitance
	Cout  = 26e-6;
	Resr  = 1e-3;

	% Compensation
	if (example == 1) 
		R1 = 100e3;
		C1 = 50e-12;
		C2 = 8e-12;   % Added to match LTspice response
	elseif (example == 2)
		R1 = 15e3;
		C1 = 330e-12;
		C2 = 8e-12;   % Added to match LTspice response
	end

elseif (example == 3)
	% Feedback reference and output voltage
	Vref = 5;
	Vout = 5;

	% Inductor (IHLP-2525CZ-01)
	L    = 4.7e-6;
	Rdcr = 40e-3;

	% Output capacitance (TDK C2012X5R1A476M 47uF, 10V, 0603)
	Cout  = 14e-6;  % capacitance under 5V DC bias
	Resr  = 3e-3;   % ESR at 1MHz

	% Compensation
	R1 = 15e3;
	C1 = 1000e-12;
	C2 = 100e-12;

end

% Error amplifier
gmEA = 0.63e-3;

% Error amplifier output resistance
% * adjusted to make the compensator calculation match 
%   the LTspice measurement at low frequency 
Rea = 5e6;

% Output load
Iout = 1;
Rout = Vout/Iout;

% -----------------------------------------------------------------
% Documentation figures
% -----------------------------------------------------------------
%
% Documentation figure notes:
% * boost the line width of the main lines to 1.0
% * use 'k--' rather than 'k:' for the dashed lines
%   and leave their width at the default 0.5
% * use 12pt fonts - as this looks good for mag+phase displayed
%   on a single figure page (mag at top, phase at bottom)
%   (change to 16pt for side-by-side figures)
% * do not use titles on the figures
%
% Generate documentation figures (set to 1)
doc_figs = 0;

% Increase the font for documentation figures
if (doc_figs==1)
	set(0, 'DefaultTextFontSize', 12);
	set(0, 'DefaultAxesFontSize', 12);
else
	set(0, 'DefaultTextFontSize', 10);
	set(0, 'DefaultAxesFontSize', 10);
end

% -----------------------------------------------------------------
% Idealized frequency response calculations
% -----------------------------------------------------------------
%
% Frequency
f = 10.^linspace(2, 6, 1000);
s = j*2*pi*f;

% LC resonance
fLC = 1/(2*pi*sqrt(L*Cout));

% Inductor-capacitor response
ZL    = s*L + Rdcr;
ZCout = 1./( 1./(1./(s*Cout) + Resr));
Kout  = ZCout./(ZCout + ZL);

% Compensator response
Zcomp = 1./(s*C2 + 1./(R1 + 1./(s*C1)) + 1/Rea);

% Output-to-compensator (feedback) gain
Kcomp = -gmEA.*Zcomp;

% -----------------------------------------------------------------
% Idealized frequency response plots
% -----------------------------------------------------------------

figure(1)
hold off
% 0dB line
semilogx([f(1) f(end)], [0 0],'k--')
hold on
% LC response
fha(1) = semilogx(f, 20*log10(abs(Kout)), 'k','LineWidth',1.0);
% LC resonanace
semilogx(fLC*[1 1], [-80 60],'k--')
% Compensator
fha(3) = semilogx(f, 20*log10(abs(Kcomp)), 'r','LineWidth',1.0);

% Switching frequency
plot(fsw*[1 1],[-80 60],'r--')

if (doc_figs==0)
	title('LTC3600 Response')
end
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

figure(2)
hold off
% 0-degree line
semilogx([f(1) f(end)], [0 0],'k--')
hold on
% -90-degree line
semilogx([f(1) f(end)], [-90 -90],'k--')
% LC response
fhb(1) = semilogx(f, angle(Kout)/pi*180, 'k','LineWidth',1.0);
% LC resonanace
semilogx(fLC*[1 1], [-180 180],'k--')
% Compensator
fhb(3) = semilogx(f, angle(Kcomp)/pi*180, 'r','LineWidth',1.0);

% Switching frequency
plot(fsw*[1 1],[-180 180],'r--')

if (doc_figs==0)
	title('LTC3600 Response')
end
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')

% -----------------------------------------------------------------
% LTSpice simulation results
% -----------------------------------------------------------------
%
% Read the frequency data and the number of steps
fs = ltspice_log_get_frequency(filename);
len = length(fs);

% Now read the measurements
gain = [
ltspice_log_get_data(filename, 'gain_db', len);
ltspice_log_get_data(filename, 'phase',   len)];

kco = [
ltspice_log_get_data(filename, 'kco_db',    len);
ltspice_log_get_data(filename, 'kco_phase', len)];

kcomp = [
ltspice_log_get_data(filename, 'kcomp_db',    len);
ltspice_log_get_data(filename, 'kcomp_phase', len)];

kpwr = [
ltspice_log_get_data(filename, 'kpwr_db',    len);
ltspice_log_get_data(filename, 'kpwr_phase', len)];

% Calculate the cross-over
m = find(fs < fsw/2);
fco_est   = interp1(gain(1,m),fs(m),0,'spline');
phase_est = interp1(fs,gain(2,:),fco_est,'spline');

% Calculate the compensation gain at fsw
kcomp_at_fsw = interp1(f,20*log10(abs(Kcomp)),fsw,'linear');

fprintf(' * LC Resonance                = %.2f kHz\n', fLC/1e3);
fprintf(' * Cross-over frequency        = %.2f kHz\n', fco_est/1e3);
fprintf(' * Cross-over phase-margin     = %.2f degrees\n', phase_est);
fprintf(' * Compensation gain at fsw    = %.2f dB\n', kcomp_at_fsw);


% Single-pole response 
% * manually adjust the gain and pole until it overlays the data
if ((example == 1) || (example == 2))
f1 = 2000;
g1 = 4;
else
f1 = 2000;
g1 = 7;
end
K1 = g1./(1+s/(2*pi*f1));

figure(1)
% Measured cross-over
plot(fco_est*[1 1],[-80 60],'b--')
semilogx(fco_est,0,'bo','LineWidth',1.0)

% Power stage single-pole response
if (doc_figs==0)
	semilogx(f, 20*log10(abs(K1)), 'c--');
end

% Measured responses
fha(4) = semilogx(fs,gain(1,:),'b','LineWidth',1.0);
semilogx(fs,gain(1,:),'bx','LineWidth',1.0)
semilogx(fs,kpwr(1,:),'kx','LineWidth',1.0)
semilogx(fs,kcomp(1,:),'rx','LineWidth',1.0)
fha(2) = semilogx(fs,kco(1,:),'g','LineWidth',1.0);
semilogx(fs,kco(1,:),'gx','LineWidth',1.0)
axis([100 1e6 -80 60])
legend(fha,'LC','Control-to-output','Compensation','Open-loop','Location','SouthWest');

figure(2)
% Measured cross-over
plot(fco_est*[1 1],[-180 180],'b--')
semilogx(fco_est,phase_est,'bo','LineWidth',1.0)

% Power stage single-pole response
if (doc_figs==0)
	semilogx(f, angle(K1)/pi*180, 'c--');
end

fhb(4) = semilogx(fs,gain(2,:),'b','LineWidth',1.0);
semilogx(fs,gain(2,:),'bx','LineWidth',1.0)
semilogx(fs,kpwr(2,:),'kx','LineWidth',1.0)
semilogx(fs,kcomp(2,:),'rx','LineWidth',1.0)
fhb(2) = semilogx(fs,kco(2,:),'g','LineWidth',1.0);
semilogx(fs,kco(2,:),'gx','LineWidth',1.0)
axis([100 1e6 -180 180])
legend(fhb,'LC','Control-to-output','Compensation','Open-loop','Location','SouthWest');
