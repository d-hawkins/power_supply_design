TDK Ferrite Beads
-----------------

The MPZ series are for power supply ferrite beads

0402 (1005M)
0603 (1608M)
0805 (2012M)

The TDK tool does not include the current rating
of the parts, however, it is listed on the data sheet.
The power rating is about 1/3W for the 0805 parts.

TDK Libraries
-------------

TDK Virtual Component Library

http://www.tdk.com/tvcl_netlist_index.php


mpz1005.zip   0402 (1005M)
mpz1608.zip   0603 (1608M)
mpz2012.zip   0805 (2012M)

c2012-class2.zip contains 

   C2012X5R1A106K  10uF, 10V
   C2012X5R1A475K 4.7uF, 10V
   C2012X5R1A225K 2.2uF, 10V

   The three capacitors used to create the PI filter.

Component Files
---------------

Extract just the files needed from each of the zip files:

1) 10uF, 10V
 
   C2012X5R1A106K_Thickness_1_25_mm (not present)
   C2012X5R1A106K_Thickness_0_85_mm (use this instead)

2) 4.7uF, 10V

   C2012X5R1A445K_Thickness_1_25_mm
   
3) 2.2uF, 10V

   C2012X5R1A225K_Thickness_1_25_mm
  
4) Ferrite bead, 2A, 600Ohm@100MHz, 100mOhm DC resistance

   MPZ2012S601A_s.mod
   
LTspice simulations
-------------------

1) TDK_ferrite_bead.asc

   Look at a couple of beads.
   The 2A bead response looks like a damped LRC resonant circuit,
   the 6A response is a zero then a pole.
   
2) TDK_pi_filter.asc

   PI-filter response for a couple of beads. This shows the
   rejection the PI-filter provides over the 1MHz to 1GHz range.
   
-------------------------------------------------------------------
0603 sized bead tests
---------------------

1A beads
--------

Digikey search results for MPZ1608 (0603) 1A beads:

MPZ1608S601A 600-Ohm@100MHz, DCR=150-mOhm
MPZ1608B471A 470-Ohm@100MHz, DCR=150-mOhm
MPZ1608D101B 100-Ohm@100MHz, DCR=150-mOhm

2A beads
--------

Digikey search results for MPZ1608 (0603) 2A beads:

MPZ1608S181A 180-Ohm@100MHz, DCR=50-mOhm
MPZ1608S121A 120-Ohm@100MHz, DCR=45-mOhm
MPZ1608Y101B 100-Ohm@100MHz, DCR=40-mOhm


The filter response was better for the beads that had
higher-resistance at 100MHz.

The criteria for bead selection should be select a part with the
appropriate current rating, and the highest resistance at 100MHz.
The resistance needs to be better than 100Ohms.