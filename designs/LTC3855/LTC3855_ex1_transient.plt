[Transient Analysis]
{
   Npanes: 3
   {
      traces: 1 {34603012,0,"I(L)"}
      X: ('m',1,0,0.0001,0.001)
      Y[0]: (' ',0,-6,3,27)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Amps: (' ',0,0,0,-6,3,27)
      Log: 0 0 0
      GridStyle: 1
   },
   {
      traces: 1 {524291,0,"V(comp)"}
      X: ('m',1,0,0.0001,0.001)
      Y[0]: (' ',1,0.3,0.1,1.5)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Volts: (' ',0,0,1,0.3,0.1,1.5)
      Log: 0 0 0
      GridStyle: 1
   },
   {
      traces: 3 {524290,0,"V(out)"} {589830,0,"0.92V"} {589830,0,"0.98V"}
      X: ('m',1,0,0.0001,0.001)
      Y[0]: ('m',0,0.918,0.006,0.984)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Volts: ('m',0,0,0,0.918,0.006,0.984)
      Log: 0 0 0
      GridStyle: 1
   }
}