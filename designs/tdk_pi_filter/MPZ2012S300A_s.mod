*----------------------------------------------------------------------
* SPICE Netlist Generated by TDK Corporation
* Copyright(C) 2012 TDK Corporation.
* All Rights Reserved.
*----------------------------------------------------------------------
* TDK P/N: MPZ2012S300A (Chip Bead)
* Property: |Z|=30ohm at 100MHz
* Size(LxWxT): 2x1.25x0.85mm, 0.079x0.049x0.033inches
* Model Type: Simple Model
* Model Generated on Mar. 26, 2012
*----------------------------------------------------------------------
* Software License Agreement:
* 1)This simulation model is being provided solely for informational
*   purposes. Please refer to the specifications of the products in
*   terms of detailed characteristics of such products.
* 2)In no event shall TDK be liable for any loss or damage arising,
*   directly or indirectly, from the use of any information contained
*   in this simulation model, including, but not limited to loss or
*   damages arising from any inaccuracies, omissions or errors in
*   connection with such information.
* 3)Any and all copyrights on this simulation model is owned by TDK
*   Corporation. Duplication or redistribution of this simulation model
*   without prior written permission from TDK Corporation is prohibited.
* 4)This simulation model is subject to any modification or change
*   without any prior notice.
* 5)The use of this simulation model is governed by the terms of the
*   software license agreement.
*
*   http://www.tdk.co.jp/etvcl/readme.htm
*----------------------------------------------------------------------
* External Node Assignments:
*
*  1 ---@@@--- 2
*
*----------------------------------------------------------------------
* Applicable Conditions:
*   Temperature = 25 degC
*   DC Bias Current = 0 A
*   Small Signal Operation
*----------------------------------------------------------------------
* LTspice modification: Added Rser=0 to L1
*----------------------------------------------------------------------
.SUBCKT MPZ2012S300A_s 1 2
C1 1 11 1.00000000E-13
L1 1 11 1.70000000E-07 Rser=0
R1 1 11 3.00000000E+01
R2 2 11 4.00000000E-03
.ENDS MPZ2012S300A_s
*----------------------------------------------------------------------
