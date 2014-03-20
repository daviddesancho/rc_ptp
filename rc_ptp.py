#!/usr/bin/env python

# This script is aimed the do the Bayesian analysis of 
# transition paths between two states given an equilibrium
# trajectory (Best & Hummer (PNAS, 2005)
# David de Sancho, 2010

import sys,os,math
import getopt
import rc_lib
from numpy import *

Usage="""
Usage:
        rc_ptp.py [-h,--help][-b blocks][-n bins][-l lower][-u upper]
                  [-L lower][-U upper][--q file]
        in which:
                 -h output help
                 -b number of blocks
                 -n number of bins
                 -l lower bound for the transition coordinate
                 -u upper bound for the transition coordinate
                 -L lower bound for binning the equilibrium coordinate 
                 -U upper bound for binning the equilibrium coordinate
                 --q file with reaction coordinate in column 2

"""

# parse input
if (len(sys.argv) < 4):
	print Usage
	sys.exit(0)
optlist,files = getopt.getopt(sys.argv[1:],'hb:l:u:L:U:n:',['help','q='])
# default values for parameters
error = False
nbins = 40
L = 0.
U = 1.
l = 0.
u = 1.
# assign values for the parameters
for opt,val in optlist:
	if opt in [ "-h", "--help" ]:
		print Usage
		sys.exit(0)
	elif opt in [ "-b"]:
		error = True
		nblocks = int(val)
	elif opt in [ "-n"]:
		nbins = int(val)
	elif opt in [ "-l"]:
		l = float(val)
	elif opt in [ "-u"]:
		u = float(val)
	elif opt in [ "-L"]:
		L = float(val)
	elif opt in [ "-U"]:
		U = float(val)
	elif opt in [ "--q"]:
		fileq = val
	else:
		print "\nUNKNOWN OPTION: %s\n\n" %opt
		print Usage
		sys.exit(1)
sys.stdout.write ("# p(TP|x) analysis (python/scipy implementation by dd363) \n")
sys.stdout.write ("# File with reaction coordinate: %s \n" %fileq)
sys.stdout.write ("#      q_low: %12.4f;   q_high: %12.4f\n" %(l,u))
sys.stdout.write ("#      bin_low: %12.4f;   bin_high: %12.4f\n" %(L,U))

# read data from file
rawdata = open(fileq,"r").readlines()
q = []
for r in rawdata:
	rs = r.split()
	q.append(float(rs[1]))
q = array(q)
lq = len(q)

# calculate p(TP|Q) for each coordinate using lower and upper bounds
nTP,lTP,bins,peq,pqTP,qTP = rc_lib.bayesian(q,l,u,L,U,nbins)
pTP = float(lTP)/lq
pTPq = zeros((nbins),float)
for i in arange(nbins):
	if (peq[i] > 0):
		pTPq[i] = pqTP[i]*pTP/peq[i]

if (error):
	# error analysis from block averages
	nTP_block,lTP_block,bins_block,peq_block,pqTP_block = \
		rc_lib.block_errors(q,l,u,L,U,nbins,nblocks)
	# output results
	sys.stdout.write ("# Total length of the equilibrium trajectory: %i \n" %lq)
	sys.stdout.write ("# Total transition path length: %i; p(TP) = %6.4e\n" %(lTP,pTP))
	sys.stdout.write ("# Number of transition paths: %i \n" %nTP)
	sys.stdout.write ("# Transition paths in each block:")
	for i in arange(nblocks):
		sys.stdout.write ("%6i" %nTP_block[i])
	sys.stdout.write ("\n")
	sys.stdout.write ("# Transition path length per block: ")
	for i in arange(nblocks):
		sys.stdout.write ("%6i" %(lTP_block[i]))
	sys.stdout.write ("\n")
	pTPq_block = zeros((nblocks,nbins),float)
	sys.stdout.write ("# p(TP) per block:")
	for i in arange(nblocks):
		pTP_block = float(lTP_block[i])/lq*nblocks
		sys.stdout.write ("%6.4f" %(pTP_block))
		for j in arange(nbins):
			if (peq_block[i,j] > 0):
				pTPq_block[i,j] = pqTP_block[i,j]/peq_block[i,j]*pTP_block
	sys.stdout.write ("\n")
	# compute averages
	peq_mean = mean(peq_block,axis=0)
	peq_std = std(peq_block,axis=0)
	pqTP_mean = mean(pqTP_block,axis=0)
	pqTP_std = std(pqTP_block,axis=0)
	pTPq_mean = mean(pTPq_block,axis=0)
	pTPq_std = std(pTPq_block,axis=0)
	sys.stdout.write ("# Bins, Peq, <Peq>, std(Peq), p(x|TP), <p(x|TP)>, std(p(x|TP)), p(TP|x), <p(TP|x)>, std(p(TP|x)) \n")
	for i in arange(nbins):
		sys.stdout.write ("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" %(bins[i],peq[i],peq_mean[i],peq_std[i],pqTP[i],pqTP_mean[i],pqTP_std[i],pTPq[i],pTPq_mean[i],pTPq_std[i]))
else:
	# output results
	sys.stdout.write ("# Total length of the equilibrium trajectory: %i \n" %lq)
	sys.stdout.write ("# Total transition path length: %i; p(TP) = %6.4e\n" %(lTP,pTP))
	sys.stdout.write ("# Number of transition paths: %i \n" %nTP)

	sys.stdout.write ("# Bins, Peq, p(x|TP), p(TP|x): \n")
	for i in arange(nbins):
		sys.stdout.write ("%10.4f %10.4f %10.4f %10.4f\n" %(bins[i],peq[i],pqTP[i],pTPq[i]))

sys.exit()
