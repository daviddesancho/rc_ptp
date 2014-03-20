#!/bin/csh -f

# This script is aimed to set the parameters for rc_ptp.py code

set proot = "1L8C"
set fileroot = "go_${proot}_Zn"
set temp = 290
while ($temp < 330)

@ temp =  ($temp + 10)

# set lower and upper bounds for q0 and q1 for TP (l,u) and binning (L,U)

set b = 10
set l = 0.03
set u = 0.7
set L = 0.20
set U = 1.

# run with options

./rc_ptp.py -b 5 -l ${l} -u ${u} -L ${L} -U ${U} \
--q=data/analysis/contacts/${fileroot}_T${temp}_qab_tmp.dat > \
data/analysis/ptp/${fileroot}_T${temp}_qab_l${l}_u${u}_error.dat
exit
end
