[Transient Analysis]
{
   Npanes: 3
   {
      traces: 1 {34603012,0,"I(L)"}
      X: ('m',1,0,0.0001,0.001)
      Y[0]: (' ',0,-2,1,9)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Amps: (' ',0,0,1,-2,1,9)
      Log: 0 0 0
      GridStyle: 1
   },
   {
      traces: 1 {524291,0,"V(comp)"}
      X: ('m',1,0,0.0001,0.001)
      Y[0]: (' ',1,0.2,0.2,2)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Volts: (' ',0,0,1,0.2,0.2,2)
      Log: 0 0 0
      GridStyle: 1
   },
   {
      traces: 3 {524290,0,"V(out)"} {589830,0,"3.135V"} {589830,0,"3.465V"}
      X: ('m',1,0,0.0001,0.001)
      Y[0]: (' ',2,3.1,0.05,3.5)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Volts: (' ',0,0,2,3.1,0.05,3.5)
      Log: 0 0 0
      GridStyle: 1
   }
}