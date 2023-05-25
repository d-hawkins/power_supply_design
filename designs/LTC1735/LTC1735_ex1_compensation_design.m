% -----------------------------------------------------------------
% LTC1735_ex1_compensation_design.m
%
% 6/16/2012 D. W. Hawkins (dwh@ovro.caltech.edu)
%
% LTC1735 Compensation Component Calculations.
%
% This m-file reads in the control-to-output response,
% Kco(s), measured using LTSpice. The component values
% for two compensation networks are then calculated
% based on user-entered poles and zeros. The user visually
% determines the locations of the poles and zeros based
% on the measured cross-over frequency and phase margin
% and the compensator gain at the switching frequency
% (which should be no greater than 0dB).
%
% By manually editing the parameters in this file, you can
% analyze the trade-off of creating a loop with good noise
% characteristics, with a desired closed-loop bandwidth.
% If you are willing to allow gain in the compensator loop
% at the switching frequency, then you can move the cross-over
% frequency higher and the low-frequency zeros higher, and
% improve the transient response of the controller.
%
% -----------------------------------------------------------------
% Example:
% --------
%
% 1) User provided inputs:
%    * fco = 50kHz
%    * R1s = 1.76
%    * fz  = 15kHz
%    * fp  = 60kHz
%
% 2) Script generated output:
%
%    * Compensator components:
%      - R1 = 40.20 kOhm
%      - C1 = 270 pF
%      - C2 = 100 pF
%      - fz = 14.66 kHz
%      - fp = 54.25 kHz
%    * Target cross-over frequency = 50.00 kHz
%    * Actual cross-over frequency = 49.75 kHz
%    * Cross-over phase-margin     = 56.84 degrees
%    * Compensation gain at fsw    = 0.02 dB
%
%   Contrast this to the script output for the original
%   components:
%
%    * Compensator components:
%      - R1 = 33.00 kOhm
%      - C1 = 330 pF
%      - C2 = 100 pF
%      - fz = 14.61 kHz
%      - fp = 62.84 kHz
%    * Target cross-over frequency = 50.00 kHz
%    * Actual cross-over frequency = 46.96 kHz
%    * Cross-over phase-margin     = 60.97 degrees
%    * Compensation gain at fsw    = 0.01 dB
%
%    So not a big difference, i.e., the original network is fine.
%
% -----------------------------------------------------------------

clear all

% Add the LTspice .log parsers to the path
% * assumes the script is called from within designs/LTC1735/
path('../matlab',path);

% -----------------------------------------------------------------
% User entered parameters for the compensation design
% -----------------------------------------------------------------

% Desired cross-over frequency
fco = 50e3; 

% R1 scale factor
% * fco determines the nominal value of R1
% * adjust this scale factor until actual fco matches desired value
R1s = 1.76; 

% Zero and pole location
fz = 15e3;  % adjust until phase margin is 60-degrees
fp = 55e3;  % adjust until gain at fsw = 0

% -----------------------------------------------------------------
% LTspice log file
% -----------------------------------------------------------------
%
% LTC1735_ex1_bode3.txt contains the measurements at all
% frequency points (including 100 and 300Hz).
%
filename = 'LTC1735_ex1_bode3.txt';

% -----------------------------------------------------------------
% LTC1735 Parameters
% -----------------------------------------------------------------
%
% These parameters match those used in the LTspice simulation
%
% Switching frequency
fsw  = 500e3;

% Feedback reference voltage
Vref = 0.8;

% Output voltage
Vout = 3.3;

% Output voltage set resistors
Ra = 3.57e3;
Rb = 1.15e3;

% Output-to-compensation feedback gain at DC
alpha_dc = Rb/(Ra+Rb);

% Inductor
L    = 2.2e-6;
Rdcr = 20e-3;

% Sense resistor
Rsns = 7e-3;

% Output capacitance
Cout = 150e-6;
Resr = 18e-3;

% Error amplifier
gmEA = 1.3e-3;

% Error amplifier output resistance
% * adjusted to make the compensator calculation match 
%   the LTspice measurement at low frequency 
Rea = 0.5e6;

% Compensation values used to measure Kco(s)
Rc1 = 33e3;
Cc1 = 330e-12;
Cc2 = 100e-12;

% Output load
% * in LTspice, the load is implemented as a 1A resistive load
%   and a 5A current load step. There is a slight difference in
%   the response for a current sink versus a resistor at low
%   frequencies, but they are the same at high frequencies. This
%   script treats the load as purely resistive.
Iout = 6;
Rout = Vout/Iout;

% =================================================================
% Measured control-to-output response
% =================================================================
%
% -----------------------------------------------------------------
% Read Kco(s)
% -----------------------------------------------------------------

% Read in the control-to-output gain measured from LTSpice
f = ltspice_log_get_frequency(filename);
s = j*2*pi*f;
len = length(f);

Kco_dB  = ltspice_log_get_data(filename, 'kco_db',    len);
Kco_deg = ltspice_log_get_data(filename, 'kco_phase', len);
Kco = 10.^(Kco_dB/20).*exp(j*pi*Kco_deg/180);

% Interpolate to determine the gain and phase at cross-over
fco_gain  = interp1(f,Kco_dB,fco,'spline');
fco_phase = interp1(f,Kco_deg,fco,'spline');

fprintf('\n');
fprintf('Control-to-output response\n');
fprintf('----------------------------\n\n');
fprintf(' * Target  fco  = %.2f kHz\n', fco/1e3);
fprintf(' * Gain at fco  = %.2f dB\n', fco_gain);
fprintf(' * Phase at fco = %.2f degrees\n', fco_phase);
fprintf('\n');

% =================================================================
% Compensation network #1: alpha*Zcomp(s)
% =================================================================
%
% This network uses the error amplifier output compensation 
% network only, and a passive-only attenuation network,
% i.e., alpha = Vref/Vout = Rb/(Ra+Rb), so there are no poles
% and zeros contributed by alpha.
%
% Zcomp is a type II network, with a series R1-C1, in parallel
% with C2.

% -----------------------------------------------------------------
% Calculation
% -----------------------------------------------------------------
%
% The compensator mid-band gain is set by R1. The closed-loop
% gain can be moved up to 0dB at fco by setting the compensator
% gain to match the inverse of the control-to-output gain at fco.
%
% i.e., 10^(-fco_gain/20) = gmEA*alpha*R1
%
% This determines the value of R1.
%
R1 = R1s*10^(-fco_gain/20)/(gmEA*alpha_dc);
R1 = resistor_standard_value(R1, 1);

% Now calculate C1 and C2 based on the target fz and fp
C1 = 1/(2*pi*R1*fz);
C2 = 1/(2*pi*R1*(fp-fz));

% Determine the nearest standard values
C1 = capacitor_standard_value(C1);
C2 = capacitor_standard_value(C2);

if (0)
R1 = Rc1;
C1 = Cc1;
C2 = Cc2;
end

% Recalculate the zero and pole locations
fz = 1/(2*pi*R1*C1);
fp = (C1+C2)/(2*pi*R1*C1*C2);

% Compensator impedance and gain
Zcomp = 1./(s*(C1+C2)).*(1+s*R1*C1)./(1+s*R1*C1*C2/(C1+C2));
Kcomp = -alpha_dc*gmEA.*Zcomp;

% Loop gain
Kloop = Kcomp.*Kco;

% Calculate the cross-over
fco_est   = interp1(20*log10(abs(Kloop)),f,0,'spline');
phase_est = interp1(f,angle(Kloop)/pi*180,fco_est,'spline');

% Calculate the compensation gain at fsw
kcomp_at_fsw = interp1(f,20*log10(abs(Kcomp)),fsw,'spline');

% Print the results
fprintf('Compensation network #1\n');
fprintf('-----------------------\n\n');
fprintf(' * Attenuator components:\n');
fprintf('   - Ra = %.2f kOhm\n', Ra/1e3);
fprintf('   - Rb = %.2f kOhm\n', Rb/1e3);
fprintf(' * Compensator components:\n');
fprintf('   - R1 = %.2f kOhm\n', R1/1e3);
fprintf('   - C1 = %d pF\n', round(C1*1e12));
fprintf('   - C2 = %d pF\n', round(C2*1e12));
fprintf('   - fz = %.2f kHz\n', fz/1e3);
fprintf('   - fp = %.2f kHz\n', fp/1e3);
fprintf(' * Target cross-over frequency = %.2f kHz\n', fco/1e3);
fprintf(' * Actual cross-over frequency = %.2f kHz\n', fco_est/1e3);
fprintf(' * Cross-over phase-margin     = %.2f degrees\n', phase_est);
fprintf(' * Compensation gain at fsw    = %.2f dB\n', kcomp_at_fsw);
fprintf('\n');

figure(1)
hold off
semilogx(f,20*log10(abs(Kco)),'g');
hold on
semilogx([100 1e6],[0 0],'k:');
semilogx(fsw*[1 1],[-60 80],'k:')
semilogx(fco*[1 1],[-60 80],'r:')
semilogx(fco,fco_gain,'gx')
semilogx(fco_est,0,'bx')
semilogx(f,20*log10(abs(Kcomp)),'r');
semilogx(f,20*log10(abs(Kloop)),'b');
axis([100 1e6 -60 80])
title('LTC1735 Compensation')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

figure(2)
hold off
semilogx(f,angle(Kco)/pi*180,'g');
hold on
semilogx([100 1e6],[0 0],'k:');
semilogx([100 1e6],[-90 -90],'k:');
semilogx(fsw*[1 1],[-180 180],'k:')
semilogx(fco*[1 1],[-180 180],'r:')
semilogx(fco,fco_phase,'gx')
semilogx(fco_est,phase_est,'bx')
semilogx(f,angle(Kcomp)/pi*180,'r');
semilogx(f,angle(Kloop)/pi*180,'b');
axis([100 1e6 -180 180])
title('LTC1735 Compensation')
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')

% =================================================================
% Compensation network #2: alpha(s)*Zcomp(s)
% =================================================================
%
% This network adds Ca across Ra to produce a frequency
% dependent attenuation alpha(s). The network creates a zero
% followed by a pole Vout/Vref later in frequency. By placing
% the zero of the feedback on the zero of the compensator
% you should get a Type III response. The pole of the compensator
% will need to be moved down by Vout/Vref to keep the gain at
% fsw at 0dB.
%
% fz_alpha = 1/(2*pi*Ra*Ca)
%
% fp_alpha = fz_alpha*Vout/Vref
%

% Use the same R1 and C1 as for the first scheme.
%
% Calculate Ca using fz_alpha = fz_comp
Ca = 1/(2*pi*Ra*fz);
Ca = capacitor_standard_value(Ca);

% Move fp in the compensator down
% * comment this out to see the effect on the gain and
%   cross-over if this pole is not moved down.
% * scale the result to allow some gain at fsw,
%   eg., scale by 0.5
C2 = 1/(2*pi*R1*(fp-fz))*Vout/Vref;
C2 = capacitor_standard_value(C2);

% Frequency dependent alpha(s)
Za = 1./(1/Ra + s*Ca);
alpha = Rb./(Za+Rb);

% Compensator impedance and gain
Zcomp = 1./(s*(C1+C2)).*(1+s*R1*C1)./(1+s*R1*C1*C2/(C1+C2));
Kcomp = -alpha*gmEA.*Zcomp;

% Loop gain
Kloop = Kcomp.*Kco;

% Calculate the cross-over
fco_est   = interp1(20*log10(abs(Kloop)),f,0,'spline');
phase_est = interp1(f,angle(Kloop)/pi*180,fco_est,'spline');

% Calculate the compensation gain at fsw
kcomp_at_fsw = interp1(f,20*log10(abs(Kcomp)),fsw,'spline');

% Print the results
fprintf('Compensation network #2\n');
fprintf('-----------------------\n\n');
fprintf(' * Attenuator components:\n');
fprintf('   - Ra = %.2f kOhm\n', Ra/1e3);
fprintf('   - Ca = %d pF\n', round(Ca*1e12));
fprintf('   - Rb = %.2f kOhm\n', Rb/1e3);
fprintf(' * Compensator components:\n');
fprintf('   - R1 = %.2f kOhm\n', R1/1e3);
fprintf('   - C1 = %d pF\n', round(C1*1e12));
fprintf('   - C2 = %d pF\n', round(C2*1e12));
fprintf(' * Target cross-over frequency = %.2f kHz\n', fco/1e3);
fprintf(' * Actual cross-over frequency = %.2f kHz\n', fco_est/1e3);
fprintf(' * Cross-over phase-margin     = %.2f degrees\n', phase_est);
fprintf(' * Compensation gain at fsw    = %.2f dB\n', kcomp_at_fsw);
fprintf('\n');

figure(1)
semilogx(fco_est,0,'c+')
semilogx(f,20*log10(abs(Kcomp)),'m--');
semilogx(f,20*log10(abs(Kloop)),'c--');

figure(2)
semilogx(fco_est,phase_est,'c+')
semilogx(f,unwrap(angle(Kcomp))/pi*180,'m--');
semilogx(f,angle(Kloop)/pi*180,'c--');
