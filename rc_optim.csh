#!/bin/csh -f

# This script is aimed to set the parameters for rc_optim.py code

set proot = "1L8C"
set fileroot = "go_${proot}_Zn"
set temp = 300

# optimization options: max / gauss

set opt = gauss

# initial value for lambda

set lambda0 = 0.001 

# set lower and upper bounds for q0 and q1 for TP (l,u) and binning (L,U)

set l0 = 0.3
set u0 = 0.65
set L0 = 0.2
set U0 = 0.8
set l1 = 5.0
set u1 = 25.0 
set L1 = 2.0
set U1 = 30.0

# run with options

./rc_optim.py --${opt} -l ${lambda0}  -n 40 \
	--l0=${l0} --u0=${u0} \
	--L0=${L0} --U0=${U0} \
	--l1=${l1} --u1=${u1} \
	--L1=${L1} --U1=${U1} \
	--q0=data/xtc/${fileroot}_T300_qab.dat \
	--q1=data/xtc/${fileroot}_T300_drmsab.dat \
	> ${fileroot}_opt_${opt}.dat
