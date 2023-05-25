LTC3600 DC1667A
---------------

The DC1667A is the evaluation board for the LTC3600.

The board has jumper options for output voltages of
0.8V, 1.2V, 1.8V, 3.3V, and a user resistor. These
jumpers select resistors that set the reference
voltage, i.e., Vref = Iset x Rset, where Iset = 50uA.

The evaluation board has a compensation network
consisting of 56k in series with 68pF, or you can
short the ITH pin to INTVCC to use internal compensation.

The LTspice transient response simulation shows this
compensation results in a transient response that meets
the typcail requirements of a 1.2V power supply (5% variation).
However, the compensation voltage is clipped at 0V,
and there is a lot of gain at the switching frequency.
The hardware should be checked to see if it matches the
simulation. Its possible that the compensation pin has
some high-frequency capacitance that slows down the
ITH pin transitions.