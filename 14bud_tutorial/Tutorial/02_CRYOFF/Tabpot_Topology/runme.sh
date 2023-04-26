#!/bin/bash
offget_charge COU,atm1,0.6645,ln8,ln12,ln16,ln20,ln24,ln28 intra.off  | grep COU | awk '{print $1, $4}' > temp_chgfile
adjust_charge temp_chgfile QQequations > chgfile && rm temp_chgfile
off2top protocol.mol intra.off template_BUU.itp BUU.itp
off2top protocol.nonbonded_list intra.off template_nonbonded_BUU_water.itp temp_nonbonded.itp
off2top protocol.nonbonded.para intra.off temp_nonbonded.itp nonbonded.itp && rm temp_nonbonded.itp

off2tab protocol.tab intra.off nonbonded.itp && cp -a BUU_OW_OW.xvg BUU.xvg
cp -a BLYPSP-4F_b0.xvg BUU_b0.xvg
