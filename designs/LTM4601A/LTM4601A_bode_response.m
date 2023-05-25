function LTM4601A_bode_response(test)
% -----------------------------------------------------------------
% LTM4601A_bode_response.m
%
% 6/18/2012 D. W. Hawkins (dwh@ovro.caltech.edu)
%
% LTM4601A Bode Response Plots.
%
% This m-file reads and plots the Bode response measured using
% the LTSpice LTM4601A_ex*_bode.asc circuit, where * = 1, 2, etc.
% The circuit writes the results to its .log file. That log file
% is then renamed to something like LTM4601A_ex*_bode1.txt so that
% LTspice does not overwrite the results. This m-file parses the
% .log file and extracts the LTspice measurements.
%
% -----------------------------------------------------------------
% Test Numbers:
% -------------
%
% Each design was used for several tests, so 'test' is used as
% the variable to configure for a specific example and test.
%
% 1) Example #1: 12V to 1.5V@12A single-phase controller design.
%
%    * CFF = 1f (not installed)
%    * LC Resonance                = 11.20 kHz
%    * Cross-over frequency        = 23.36 kHz
%    * Cross-over phase-margin     = 63.77 degrees
%    * Compensation gain at fsw    = 6.88 dB
%
% 2) Example #1: 12V to 1.5V@12A single-phase controller design.
%
%    * CFF = 330pF (between Vout and FB)
%    * LC Resonance                = 11.20 kHz
%    * Cross-over frequency        = 53.44 kHz
%    * Cross-over phase-margin     = 101.45 degrees
%    * Compensation gain at fsw    = 14.85 dB
%
% 3) Example #2: 12V to 1.5V@24A dual-phase controller design.
%
%    * CFF = 1f (not installed)
%    * LC Resonance                = 7.92 kHz
%    * Cross-over frequency        = 23.41 kHz
%    * Cross-over phase-margin     = 64.52 degrees
%    * Compensation gain at fsw    = 8.00 dB
%
% 4) Example #2: 12V to 1.5V@24A dual-phase controller design.
%
%    * CFF = 330pF
%    * LC Resonance                = 7.92 kHz
%    * Cross-over frequency        = 55.32 kHz
%    * Cross-over phase-margin     = 104.33 degrees
%    * Compensation gain at fsw    = 15.96 dB
%
% So the Bode response for the single- and dual-phase controller
% is pretty similar. Now change the output voltage to 0.95V.
%
% 5) Example #3: 12V to 0.95V@12A single-phase controller design.
%
%    * CFF = 1f (not installed)
%    * LC Resonance                = 7.03 kHz
%    * Cross-over frequency        = 16.76 kHz
%    * Cross-over phase-margin     = 54.29 degrees
%    * Compensation gain at fsw    = 11.97 dB
%
% 6) Example #3: 12V to 0.95V@12A single-phase controller design.
%
%    * CFF = 330p
%    * LC Resonance                = 7.03 kHz
%    * Cross-over frequency        = 22.23 kHz
%    * Cross-over phase-margin     = 75.07 degrees
%    * Compensation gain at fsw    = 15.96 dB
%
% I didn't bother with an Example #4: 0.95V dual-phase controller,
% since the single-phase controller could not meet the output
% voltage regulation requirements.
%
% -----------------------------------------------------------------

%clear all

% Add the LTspice .log parsers to the path
% * assumes the script is called from within designs/LTM4601A/
path('../matlab',path);

% -----------------------------------------------------------------
% LTspice log file
% -----------------------------------------------------------------
%
% Select test 
if (nargin == 0)
	% Default to test 1
	test = 1; 
end

fprintf('Loading test#%d data\n', test)
if (test == 1) 
	filename = 'LTM4601A_ex1_bode1.txt';
elseif (test == 2)
	filename = 'LTM4601A_ex1_bode2.txt';
elseif (test == 3)
	filename = 'LTM4601A_ex2_bode1.txt';
elseif (test == 4)
	filename = 'LTM4601A_ex2_bode2.txt';
elseif (test == 5)
	filename = 'LTM4601A_ex3_bode1.txt';
elseif (test == 6)
	filename = 'LTM4601A_ex3_bode2.txt';
else
	error('Invalid test number')
end

% -----------------------------------------------------------------
% LTM4601A Parameters
% -----------------------------------------------------------------
%
% These parameters match those used in the LTspice simulation
%
% Switching frequency
fsw  = 850e3;

% Feedback reference voltage
Vref = 0.6;

% Output voltage
if (test <= 4)
	Vout = 1.5;
else
	Vout = 0.95;
end

% Output voltage set resistors
Ra = 60.4e3;
if (test <= 4)
	Rb = 40.2e3;
else
	Rb = 103.5e3;
end
if (test == 2) || (test == 4) || (test == 6)
	Ca = 330e-12;
end

% Inductor
L    = 0.47e-6;
Rdcr = 2.0e-3; % (guess)

% Output capacitance
if (test == 1) || (test == 2)
	% Single-phase 1.5V output
	Cout  = 330e-6;
	Resr  = 6e-3;
	Cout2 = 100e-6;
	Resr2 = 2e-3;
elseif (test == 3) || (test == 4)
	% Dual-phase 1.5V output
	Cout  = 2*330e-6;
	Resr  = 6e-3/2;
	Cout2 = 2*100e-6;
	Resr2 = 2e-3/2;
elseif (test == 5) || (test == 6)
	% Single-phase 0.95V output
	Cout  = 3*330e-6;
	Resr  = 6e-3/3;
	Cout2 = 100e-6;
	Resr2 = 2e-3;
end

% Error amplifier (guess)
gmEA = 1.5e-3;

% Error amplifier output resistance
% * adjusted to make the compensator calculation match 
%   the LTspice measurement at low frequency 
Rea = 0.7e6;

% Compensation (guess based on Bode measurements - also depends on gmEA)
R1 = 9e3;
C1 = 1000e-12;
C2 = 39e-12;

% Output load (per phase)
% * in LTspice, the load is implemented as a 1A resistive load
%   and a 11A current load step. There is a slight difference in
%   the response for a current sink versus a resistor at low
%   frequencies, but they are the same at high frequencies. This
%   script treats the load as purely resistive.
Iout = 12;
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
f = 10.^linspace(2, log10(1e6), 1000);
s = j*2*pi*f;

% LC resonance
fLC = 1/(2*pi*sqrt(L*(Cout+Cout2)));

% Inductor-capacitor response
ZL    = s*L + Rdcr;
ZCout = 1./( 1./(1./(s*Cout) + Resr) + 1./(1./(s*Cout2) + Resr2) );
Kout  = ZCout./(ZCout + ZL);

% Compensator response
Zcomp = 1./(s*C2 + 1./(R1 + 1./(s*C1)) + 1/Rea);

% Feedback attenuation
% * nominally alpha = Vref/Vout, but it depends on Ra/Rb values
if (test == 1) || (test == 3) || (test == 5)
	alpha = Rb/(Ra+Rb);
else
	Za = 1./(1/Ra + s*Ca);
	alpha = Rb./(Za + Rb);
end

% Output-to-compensator (feedback) gain
Kcomp = -alpha*gmEA.*Zcomp;

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
	title('LTM4601A Response')
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
	title('LTM4601A Response')
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
if (test == 1) || (test == 2) 
	f1 = 300;
	g1 = 12;
elseif (test == 3) || (test == 4)
	f1 = 90;
	g1 = 45;
elseif (test == 5) || (test == 6)
	f1 = 180;
	g1 = 8;
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
semilogx(fs,kcomp(2,:),'rx','LineWidth',1.0)
fhb(2) = semilogx(fs,kco(2,:),'g','LineWidth',1.0);
semilogx(fs,kco(2,:),'gx','LineWidth',1.0)
axis([100 1e6 -180 180])
legend(fhb,'LC','Control-to-output','Compensation','Open-loop','Location','SouthWest');
