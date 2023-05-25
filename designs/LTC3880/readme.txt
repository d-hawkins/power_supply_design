LTC3880 Designs
---------------

Design notes/comments.

-------------------------------------------------------------------------------
Fairchild MOSFET model edits
----------------------------

The Fairchild MOSFET models have been designed for use with PSpice.
The models need to be edited for use with LTspice.

1) Change the MOSFET model LEVEL from 7 to 8 to select the BSim3 model.

2) Add Rser=0 to each of the inductors. Alternatively set the LTspice
   control panel Hacks! option to set Rser=0.
   
-------------------------------------------------------------------------------
12V to 0.95V@45A 3-phase design
-------------------------------

Simulations are created for 1-, 2-, and 3-phases to compare
the efficiency and simulation times. Since the efficiency
remains approximately the same, with the same per-phase
power dissipation, a 2-phase controller can be used during
simulation tests (since the simulation runs faster).

Component SpiceLines: 
  VIN_ON=9 VIN_OFF=8 Vout_0=0.95 Vout_1=0.95 Ilim0_range=1 Ilim1_range=1 OC_limit0=1 OC_limit1=1 Mode_ll=2 Fault_response=0 Retry_delay=.1m
  Freq=500K PHs_0=0 PHs_1=180 Ton0_delay=.3m Ton0_rise=.5m Ton1_delay=.3m Ton1_rise=.5m Toff0_delay=.2m Toff0_fall=.3m Toff1_delay=.2m Toff1_fall=.3m Vout0_range=1 Vout1_range=1

The Altera Stratix IV GT FPGA core power supply voltage variation
specification is; Vout = 950mV +/- 30mV

How well does each design do?

LTC3880_3phase_design1
 * simulation time = 269s (5 mins)
 * Vout(avg) = 949mV 
 * Dip  to 894mV (26mV violation)
 * Peak to 994mV (14mV violation)

LTC3880_3phase_design2
 * simulation time = 583s (10 mins)
 * Vout(avg) = 949mV 
 * Dip  to 907mV (13mV violation)
 * Peak to 989mV ( 9mV violation)

LTC3880_3phase_design3
 * 3-phase controller phase setup:
   - use configuration resistors to setup the first
     part with 120-degree and 240-degree phasing, and
     the second part to have 0-degree and 180-degree
     phasing (the 180-degree phase is not used)
   - use a common RUN signal between the two parts
   - look at the 4-phase circuit on p109 in the data
     sheet; the VOUT_CFG and VTRIM_CFG pins are not in
     common between parts, i.e., each device needs its
     own set of resistors (otherwise the output voltage
     is not correct), and the figure shows the
     compensation network repeated twice (which is
     pointless, so don't copy that).
 * simulation time = 1426s (24 mins), 1544s (26 mins)
 * Vout(avg) = 949mV
 * Dip  to 907mV (13mV violation for 2.7us)
 * Peak to 990mV (10mV violation for 3.5us)
   So the violations are fairly symmetric. 
 * dVout = 2mV (peak-to-peak ripple)
   dVith = 6mV 
   There is a little too much high-frequency gain, 
   the compensation network zero could be moved
   down to fix that.

-------------------------------------------------------------------------------
12V to 0.95V@40A 2-phase design
-------------------------------

LTC3880_2phase_design1
  * Relative to the 3-phase design
     - change the inductor and current sense
     - change the load step
     - change the Ilim_range to low
     - change the compensation components
     
  * Transient response
     - Vith voltage = 0.7V at  1A load
                    = 1.5V at 20A load
                    
     - Vout ripple voltage 13mV
       Vith ripple voltage 20mV
       
     - Dip to   881mV (39mV violation for 10us)
       Peak to 1007mV (27mV violation for 10us)
       
       The loop bandwidth needs to be increased - but this is fine for now.
  
  * Efficiency analysis
     