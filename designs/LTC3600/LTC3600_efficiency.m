function LTC3600_efficiency(filename)
% -----------------------------------------------------------------
% LTC3600_efficiency.m
%
% 9/20/2012 D. W. Hawkins (dwh@ovro.caltech.edu)
%
% LTC3600 Efficiency Plots.
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
% LTC3600_efficiency.txt contains the measurements at the
% different load current steps.
%
%filename = 'LTC3600_efficiency.txt';

% -----------------------------------------------------------------
% LTSpice simulation results
% -----------------------------------------------------------------
%
% Read the output current and the number of steps
iout = ltspice_log_get_iout(filename);
len = length(iout);

% Efficiency
eff = ltspice_log_get_data(filename, 'eff', len);

% Power measurements
P_in    = ltspice_log_get_data(filename, 'p_in',    len);
P_out   = ltspice_log_get_data(filename, 'p_out',   len);
P_U1    = ltspice_log_get_data(filename, 'p_u1',    len);
P_Cout  = ltspice_log_get_data(filename, 'p_cout',  len);
P_L     = ltspice_log_get_data(filename, 'p_l',     len);

figure(1)
hold off
plot(iout,eff,'LineWidth',1.0)
hold on
plot(iout,eff,'x')
xlabel('Load Current (A)')
ylabel('Efficiency (%)')
title('LTC3600 12V to 2.5V@1.5A Power Supply')
grid on
axis([-inf inf 50 100])

figure(2)
hold off
clear h
h(1)  = plot(iout,P_in-P_out, 'k','LineWidth',1.0);
hold on
% Show that the component powers add close to that of the 
% whole supply (not all powers are written to the log file)
plot(iout,P_Cout+P_L+P_U1,'kx')
%
% Component losses
h(2) = plot(iout,P_U1,    'b','LineWidth',1.0);
h(3) = plot(iout,P_L,     'r','LineWidth',1.0);
h(4) = plot(iout,P_Cout,  'g','LineWidth',1.0);
legend(h,'Total','U1','L','Cout','Location','NorthWest')
xlabel('Load Current (A)')
ylabel('Power Dissipation (W)')
title('LTC3600 12V to 2.5V@1.5A Power Supply')
grid on
axis([-inf inf -0.01 0.5])
