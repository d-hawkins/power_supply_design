[Transient Analysis]
{
   Npanes: 3
   Active Pane: 1
   {
      traces: 1 {34603011,0,"I(L1)"}
      X: ('m',1,0,0.0001,0.001)
      Y[0]: (' ',1,-1.2,0.3,2.1)
      Y[1]: (' ',0,1e+308,5,-1e+308)
      Amps: (' ',0,0,1,-1.2,0.3,2.1)
      Log: 0 0 0
      GridStyle: 1
   },
   {
      traces: 1 {524292,0,"V(ith)"}
      X: ('m',1,0,0.0001,0.001)
      Y[0]: (' ',1,0,0.2,2.2)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Volts: (' ',0,0,1,0,0.2,2.2)
      Log: 0 0 0
      GridStyle: 1
   },
   {
      traces: 3 {524290,0,"V(out)"} {589830,0,"1.14V"} {589830,0,"1.26V"}
      X: ('m',1,0,0.0001,0.001)
      Y[0]: (' ',2,1.13,0.01,1.27)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Volts: (' ',0,0,2,1.13,0.01,1.27)
      Log: 0 0 0
      GridStyle: 1
   }
}
