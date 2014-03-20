#!/usr/bin/env python

# Useful functions for reaction coordinate optimization

import sys,os,math
from numpy import *
from scipy import optimize

def gaussian(param,x):
	# gaussian function of x with parameters mu, sigma and height
	mu = param[0]
	sigma = param[1]
	height = param[2]
	gauss = height*exp(-((x-mu)/sigma)**2)
	return gauss

def gauss_sls(param,x,y):
	# sls function for gaussian optimization purposes
	gauss = gaussian(param,x)
	diff = sum((gauss - y)**2)
	return diff

def bayesian(q,lb,ub,LB,UB,nbins):
	# Bayesian analysis of equilibrium trajectories 
	# to extract properties of transition path ensembles
	# Best & Hummer (PNAS, 2005)
	# analyze trajectory considering boundaries
	isTP = False
	lq = len(q)
	qeq = []
	iqTP = []
	qTP = [] # q for transition paths
	iqmaybeTP = []
	qmaybeTP = [] # q for possible transition paths
	state = '-1'
	ntp = 0
	for i in arange(lq):
		# check if coordinate is between binning bounds
		if (q[i] <= UB) and (q[i] >= LB):
			qeq.append(q[i])
			# assign state when beyond boundaries
			if (q[i] > ub): # state "1"
				if (state == 0):
					ntp +=1
					qTP.extend(qmaybeTP)	# transition path?
					iqTP.extend(iqmaybeTP)
				state = 1
				qmaybeTP = []
				iqmaybeTP = []
			elif (q[i] < lb): # state "0"
				if (state == 1):
					ntp +=1
					qTP.extend(qmaybeTP) # transition path?
					iqTP.extend(iqmaybeTP)
				state = 0
				qmaybeTP = []
				iqmaybeTP = []
			else:
			# save values of q that may be in a TP
				if (state == 1) and (q[i] < ub):
					qmaybeTP.append(q[i])
					iqmaybeTP.append(i+1)
				elif (state == 0) and (q[i] > lb):
					qmaybeTP.append(q[i])
					iqmaybeTP.append(i+1)
	# convert lists to numpy arrays
	qeq = array(qeq)
	qTP = array(qTP)
	iqTP = array(iqTP)
	for i in arange(len(qTP)):
		print iqTP[i],qTP[i]
	sys.exit()
	ltp = len(qTP)
	if (ltp == 0):
		print " No transitions in range"
		sys.exit()
	hist_eq,bin_edges = histogram(qeq,bins=nbins,range=(LB,UB),normed=True)
	hist_tp,bin_edges = histogram(qTP,bins=nbins,range=(LB,UB),normed=True)
	bin_center = zeros(nbins)
	for i in arange(nbins):
		bin_center[i] = (bin_edges[i] + bin_edges[i+1])/2.
	return ntp,ltp,bin_center,hist_eq,hist_tp,qTP

def block_errors(q,l,u,L,U,nbins,nblocks):
	# Carry out error analysis using block averages
	lq = len(q)
	block = lq/nblocks
	nTP_block = []
	lTP_block  = []
	bins_block = []
	peq_block = []
	pqTP_block = []
	qTP_block = []
 	for i in arange(nblocks):
		ib = i*block
		ie = (i+1)*block - 1
		qblock = q[ib:ie]
		nTP,lTP,bins,peq,pqTP,qTP = bayesian(qblock,l,u,L,U,nbins)
		nTP_block.append(nTP)
		lTP_block.append(lTP)
		bins_block.append(bins)
		peq_block.append(peq)
		pqTP_block.append(pqTP)
	nTP = array(nTP_block)
	lTP = array(lTP_block)
	bins = array(bins_block[0])
	peq = array(peq_block)
	pqTP = array(pqTP_block)
	return nTP,lTP,bins,peq,pqTP

def rc_comb(lmbd,q0,lb0,ub0,LB0,UB0,q1,lb1,ub1,LB1,UB1,nbins):
	# Obtain new reaction coordinate variationally 
	q = q0 + lmbd*q1
	# define upper and lower bounds for new coordinate
	lb = lb0 + lmbd*lb1
	ub = ub0 + lmbd*ub1
	LB = LB0 + lmbd*LB1
	UB = UB0 + lmbd*UB1
	nTP,lTP,bins,peq,pqTP,qTP = bayesian(q,lb,ub,LB,UB,nbins)
	lq = len(q)
	pTP = float(lTP)/lq
	pTPq = zeros((nbins),float)
	for i in arange(nbins):
		if (peq[i] > 0):
			pTPq[i] = pqTP[i]*pTP/peq[i]
	pTPqmax = max(pTPq)
	mpTPqmax = -1*pTPqmax
	return mpTPqmax

def rc_pTPmax(lmbd,q0,lb0,ub0,LB0,UB0,q1,lb1,ub1,LB1,UB1,nbins):
	# Optimization of the maximum of p(TP|q)
	# Obtain new reaction coordinate variationally 
	q = q0 + lmbd*q1
	# define upper and lower bounds for new coordinate
	lb = lb0 + lmbd*lb1
	ub = ub0 + lmbd*ub1
	LB = LB0 + lmbd*LB1
	UB = UB0 + lmbd*UB1
	nTP,lTP,bins,peq,pqTP,qTP = bayesian(q,lb,ub,LB,UB,nbins)
	lq = len(q)
	pTP = float(lTP)/lq
	pTPq = zeros((nbins),float)
	for i in arange(nbins):
		if (peq[i] > 0):
			pTPq[i] = pqTP[i]*pTP/peq[i]
	pTPqmax = max(pTPq)
	mpTPqmax = -1*pTPqmax
	return mpTPqmax

def rc_pTPgauss(lmbd,q0,lb0,ub0,LB0,UB0,q1,lb1,ub1,LB1,UB1,nbins):
	# Optimization of the maximum of a gaussian fit to p(TP|q)
	# Obtain new reaction coordinate variationally 
	q = q0 + lmbd*q1
#	for i in arange(len(q)):
#		print q[i]
	# define upper and lower bounds for new coordinate
	lb = lb0 + lmbd*lb1
	ub = ub0 + lmbd*ub1
	LB = LB0 + lmbd*LB1
	UB = UB0 + lmbd*UB1
	nTP,lTP,bins,peq,pqTP,qTP = bayesian(q,lb,ub,LB,UB,nbins)
	lq = len(q)
	pTP = float(lTP)/lq
	pTPq = zeros((nbins),float)
	for i in arange(nbins):
		if (peq[i] > 0):
			pTPq[i] = pqTP[i]*pTP/peq[i]
	# giving initial parameters for gaussian fit
	mu = mean(qTP)
	sigma = std(qTP)
	height = max(pTPq) 
	p0 = []
	p0.append(mu)
	p0.append(sigma)
	p0.append(height)
	p0 = array(p0)
	# optimize gaussian
	opt_results = optimize.fmin(gauss_sls,p0,\
		args=(bins,pTPq), \
                xtol=1e-8, \
                ftol=1e-8, \
                full_output=1, \
                disp=0)
	popt = opt_results[0]
	mheight_opt  = -1*popt[2]
	return mheight_opt
