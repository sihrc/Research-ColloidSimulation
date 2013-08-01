"""g(r) program to validate array_main.py
This was adapted from gofrone.m from Rebecca J. Christianson"""
import os, sys, csv, random
import numpy as np
import matplotlib.pyplot as plt


"""Importing the data"""
datafile = open(sys.argv[1],'r')

data = csv.reader(datafile)
"""Reading the parameter data"""
info = data.next()
num_part = int(info[0])
total_lines = float(info[1])/float(info[2]) * num_part
radius = float(info[3])
side = float(info[4])

current = np.zeros(num_part,3)
steps = 100
dx = .1;
part = 0
while data.line_num < 61000#total_lines + 1:
	part = part%num_part
	print data.line_num
	coord = data.next()
	current[part,0] = float(coord[0])
	current[part,1] = float(coord[1])
	current[part,2] = float(coord[2])
	if part == num_part - 1:
		minmax = np.ptp(current,axis = 0)
		volume = minmax.prod()
		density = float(num_part/volume)

		r = np.zeros(num_part,1)
		gofr = np.zeros(steps,4)

		for i in range(steps):
			outer_radius = i*dx
			inner_radius = radius - dx
			dist = np.sum(current**2, axis = 1)**.5
			count = len(np.where[(dist<outer_radius)*(dist>inner_radius)][0])
			gofr[i,0] = count
			gofr[i,1] = dx
			gofr[i,2] = inner_radius
			gofr[i,3] = outer_radius

	part += 1
