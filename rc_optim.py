#!/usr/bin/env python

# This script is aimed to optimize the reaction coordinate 
# q = q0 + lambda*q1 using the Bayesian approach to analyze 
# the transition paths (Best & Hummer (PNAS, 2005)
# David de Sancho, 2010

import sys,os,math
import getopt
import rc_lib
from numpy import *
from scipy import optimize

Usage="""
Usage:
	
        rc_optim.py [-h,--help][-l lambda]
                    [--l0 lower][--u0 upper][--L0 lower][--U0 upper]
                    [--q0 file][--l1 lower][--u1 upper][--L1 lower]
                    [--U1 upper][--q1 file]

        in which:
                 -h output help
                 -n number of bins
                 --gauss: optimize maximum of a gaussian fit to p(TP|Q)
                 --max: optimize maximum of p(TP|Q)
                 -l lambda to be optimized 
                 -li lower bound for the transition coordinate
                 -ui upper bound for the transition coordinate
                 -Li lower bound for binning the equilibrium coordinate 
                 -Ui upper bound for binning the equilibrium coordinate
                 --qi file with reaction coordinate in column 2

"""

# parse input
if (len(sys.argv) < 7):
	print Usage
	sys.exit(0)

optlist,files = getopt.getopt(sys.argv[1:],'hl:n:',['help','gauss','max','l0=','l1=','u0=','u1=','L0=','L1=','U0=','U1=','q0=','q1='])
# default values
optim = 'max'
maxiter = 100
nbins = 40
lmbd0 = 0.
LB0 = 0.
UB0 = 1.
lb0 = 0.
ub0 = 1.
LB1 = 0.
UB1 = 1.
lb1 = 0.
ub1 = 1.
for opt,val in optlist:
	if opt in [ "-h", "--help" ]:
		print Usage
		sys.exit(0)
	elif opt in ["--gauss"]:
		optim = 'gauss'
	elif opt in ["--max"]:
		optim = 'max'
	elif opt in [ "-l"]:
		lmbd0 = float(val)
	elif opt in [ "-n"]:
		nbins = int(val)
	elif opt in [ "--l0"]:
		lb0 = float(val)
	elif opt in [ "--u0"]:
		ub0 = float(val)
	elif opt in [ "--L0"]:
		LB0 = float(val)
	elif opt in [ "--U0"]:
		UB0 = float(val)
	elif opt in [ "--l1"]:
		lb1 = float(val)
	elif opt in [ "--u1"]:
		ub1 = float(val)
	elif opt in [ "--L1"]:
		LB1 = float(val)
	elif opt in [ "--U1"]:
		UB1 = float(val)
	elif opt in [ "--q0"]:
		fileq0 = val
	elif opt in [ "--q1"]:
		fileq1 = val
	else:
		print "\nUNKNOWN OPTION: %s\n\n" %opt
		print Usage
		sys.exit(1)
sys.stdout.write ("# p(TP|x) variational rc optimization \n")
sys.stdout.write ("# (python/scipy implementation by dd363)\n")
sys.stdout.write ("#    File with reaction coordinate q0: %s \n" %fileq0)
sys.stdout.write ("#    File with reaction coordinate q1: %s \n" %fileq1)

# read values of the order parameters q0 and q1 from file
rawdata = open(fileq0,"r").readlines()
q0 = []
for r in rawdata:
	rs = r.split()
	q0.append(float(rs[1]))
q0 = array(q0)
lq0 = len(q0)
rawdata = open(fileq1,"r").readlines()
q1 = []
for r in rawdata:
	rs = r.split()
	q1.append(float(rs[1]))
q1 = array(q1)
lq1 = len(q1)
if (lq0 != lq1):
	print "ERROR: Files have different lenghts!!!!"
	sys.exit(1)
else:
	lq = lq0
	sys.stdout.write ("#    Total length of the equilibrium trajectory: %i \n" %lq0)

# calculate p(TP|x) for q0 using input lower and upper bounds
nTP0,lTP0,bins0,peq0,pqTP0,qTP0 = rc_lib.bayesian(q0,lb0,ub0,LB0,UB0,nbins)
pTP0 = float(lTP0)/lq
pTPq0 = zeros((nbins),float)
for i in arange(nbins):
	if (peq0[i] > 0):
		pTPq0[i] = pqTP0[i]*pTP0/peq0[i]
pTPq0_max = max(pTPq0)
sys.stdout.write ("# TP analysis for coordinate q0:\n")
sys.stdout.write ("#    Total transition path length (q0): %i; p(TP) = %6.4e\n" %(lTP0,pTP0))
sys.stdout.write ("#    Number of transition paths (q0): %i \n" %nTP0)
sys.stdout.write ("#    lb0 = %8.6f; ub0 = %8.6f;LB0 = %8.6f;UB0 = %8.6f \n" %(lb0,ub0,LB0,UB0))

if (optim =='gauss'):
	# fit p(TP|q0) to gaussian
	mu = mean(qTP0)
	sigma = std(qTP0)
	height = 0.2
	p0 = []
	p0.append(mu)
	p0.append(sigma)
	p0.append(height)
	p0 = array(p0)
	opt_results = optimize.fmin(rc_lib.gauss_sls,p0,\
		args=(bins0,pTPq0),\
		xtol=1e-6, \
		ftol=1e-6, \
		full_output=1, \
		disp=0)
	p0opt = opt_results[0]
	gauss0 = rc_lib.gaussian(p0opt,bins0)
	pTPq0fit_max = max(gauss0)

# optimize p(TP|q) for q = q0 + lambda*q1
#help(optimize.fmin)
if (optim == 'max'):
	# optimize maximum of p(TP|q) 
	#mpTPqmax,q,lb,ub,LB,UB = rc_lib.rc_pTPmax(lmbd_opt,q0,lb0,ub0,LB0,UB0,q1,lb1,ub1,LB1,UB1,nbins)
	opt_results = optimize.fmin(rc_lib.rc_pTPmax,lmbd0, \
		args=(q0,lb0,ub0,LB0,UB0,q1,lb1,ub1,LB1,UB1,nbins), \
		xtol=1e-10, \
		ftol=1e-10, \
		maxiter=None, \
		maxfun=None, \
		full_output=1, \
		disp=0, \
		retall=0, \
		callback=None)
	lmbd_opt = opt_results[0]
	pTPq_opt = -1*opt_results[1]
	iter = opt_results[2]
	evals = opt_results[3]
	sys.stdout.write ("# RC optimization results: q = q0 + lambda*q1:\n")
	sys.stdout.write ("# Criterion: maximization of p(TP|q) \n")
	sys.stdout.write ("#    p(TP|q)0 = %8.6e; lambda0 = %8.6e \n" %(pTPq0_max,lmbd0))
	sys.stdout.write ("#    p(TP|q)max = %8.6e; lambda = %8.6e \n" %(pTPq_opt,lmbd_opt))
	sys.stdout.write ("#    iterations = %i \n" %iter)
	sys.stdout.write ("#    evaluations = %i \n" %evals)

elif (optim == 'gauss'):
	# optimize maximum of gaussian fit to p(TP|q)
	#a = rc_lib.rc_pTPgauss(lmbd0,q0,lb0,ub0,LB0,UB0,q1,lb1,ub1,LB1,UB1,nbins)
	# optimize p(TP|x)
	popt = opt_results[0]
        opt_results = optimize.fmin(rc_lib.rc_pTPgauss,lmbd0, \
                args=(q0,lb0,ub0,LB0,UB0,q1,lb1,ub1,LB1,UB1,nbins), \
                xtol=1e-10, \
                ftol=1e-10, \
                maxiter=maxiter, \
                maxfun=None, \
                full_output=1, \
                disp=0, \
                retall=0, \
                callback=None)
	lmbd_opt = opt_results[0]
	pTPq_opt = -1*opt_results[1]
	iter = opt_results[2]
	evals = opt_results[3]
	sys.stdout.write ("# RC optimization results: q = q0 + lambda*q1:\n")
	sys.stdout.write ("# Criterion: maximum of gaussian fit to p(TP|q) \n")
	sys.stdout.write ("#    p(TP|q)0 = %8.6e; gauss0 = %8.6e; lambda0 = %8.6e \n" %(max(pTPq0),pTPq0fit_max,lmbd0))
	sys.stdout.write ("#    gauss_max = %8.6e; lambda = %8.6e \n" %(pTPq_opt,lmbd_opt))
	sys.stdout.write ("#    iterations = %i \n" %iter)
	sys.stdout.write ("#    evaluations = %i \n" %evals)
qopt = q0 + lmbd_opt*q1
lb = lb0 + lmbd_opt*lb1
ub = ub0 + lmbd_opt*ub1
LB = LB0 + lmbd_opt*LB1
UB = UB0 + lmbd_opt*UB1

# recalculate pTP for optimized reaction coordinate
nTP,lTP,bins,peq,pqTP,qTP = rc_lib.bayesian(qopt,lb,ub,LB,UB,nbins)
pTP = float(lTP)/lq
pTPq = zeros((nbins),float)
for i in arange(nbins):
	if (peq[i] > 0):
		pTPq[i] = pqTP[i]*pTP/peq[i]
sys.stdout.write ("# TP analysis for optimized coordinate q = q0 + lambda*q1:\n")
sys.stdout.write ("#    Total transition path length (qopt): %i; p(TP) = %6.4e\n" %(lTP,pTP))
sys.stdout.write ("#    Number of transition paths (qopt): %i \n" %nTP)
sys.stdout.write ("#    lb = %8.6f; ub = %8.6f;LB = %8.6f;UB = %8.6f \n" %(lb,ub,LB,UB))

if optim == 'max':
	# output results for q0 and q
	for i in arange(nbins):
		sys.stdout.write ("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" %(bins0[i],peq0[i],pqTP0[i],pTPq0[i],bins[i],peq[i],pqTP[i],pTPq[i]))
elif optim == 'gauss':
	# fit to gaussian
	mu = mean(qTP)
	sigma = std(qTP)
	height = 0.2
	p0 = []
	p0.append(mu)
	p0.append(sigma)
	p0.append(height)
	p0 = array(p0)
	opt_results = optimize.fmin(rc_lib.gauss_sls,p0,\
		args=(bins,pTPq),\
		xtol=1e-6, \
		ftol=1e-6, \
		full_output=1, \
		disp=0)
	popt = opt_results[0]
        pTPq_opt = -1*opt_results[1]
	gauss = rc_lib.gaussian(popt,bins)
	# output results for q0, q and gaussian fit
	sys.stdout.write ("#  bins       peq0       pqTP0      pTPq0     gauss     bins       peq        pqTP       pTPq      gauss\n")
	for i in arange(nbins):
		sys.stdout.write ("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" %(bins0[i],peq0[i],pqTP0[i],pTPq0[i],gauss0[i],bins[i],peq[i],pqTP[i],pTPq[i],gauss[i]))

sys.exit(0)
