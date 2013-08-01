"""g(r) program to validate array_main.py
This was adapted from gofrone.m from Rebecca J. Christianson"""
import os, sys, csv, random, time
import numpy as np
import numexpr as ne
import matplotlib.pyplot as plt


def distance(current):
	x1 = current[:,0]
	y1 = current[:,1]
	z1 = current[:,2]

	x = x1 - x1[:,np.newaxis]
	y = y1 - y1[:,np.newaxis]
	z = z1 - z1[:,np.newaxis]
	return ne.evaluate("(x**2 + y**2 + z**2)**.5")


pi = np.pi
"""Importing the data"""
datafile = open(sys.argv[1],'r')
result = open('gofr_results.csv','wb')
data = csv.reader(datafile)
"""Reading the parameter data"""
coord = data.next()
num_part = int(coord[0])
frames = float(coord[5])
total_lines = float(coord[1])/float(coord[2]) * num_part/frames
radius = float(coord[3])
side = float(coord[4])


part = 0
steps = 100
dr = .001

current = np.zeros(shape = (num_part,3))
gofrset = np.zeros(shape = (steps,4))
#results = csv.writer(result)
while 1:
	coord = data.next()
	if len(coord) == 5:
		break
	current[part,0] = float(coord[0])
	current[part,1] = float(coord[1])
	current[part,2] = float(coord[2])
	part += 1
	if part == num_part:
		part = 0
		mini = np.amin(current,axis = 0)
		maxi = np.amax(current,axis = 0)
		volume = (maxi - mini).prod()
		density = float(num_part/volume)
		r = distance(current)
		for i in range(steps):
			radius = i*dr
			inn = radius - dr
			count = 0
			vshell = 0
			minchange = current < mini + radius 
			maxchange = current > maxi - radius
			mini_h = (current - mini)[minchange]
			maxi_h = (maxi - current)[maxchange]
			vshell += np.sum(ne.evaluate("1/3.0*pi*(2*(radius**3 - inn**3) + 3*mini_h*(radius**2 - inn**2))"))
			vshell += np.sum(ne.evaluate("1/3.0*pi*(2*(radius**3 - inn**3) + 3*maxi_h*(radius**2 - inn**2))"))
			vshell += np.sum(((-minchange + -maxchange))*ne.evaluate("4/3.0*pi*(radius**3 - inn**3)"))
			count += np.sum((r<radius) * (r>inn))
			gofrset[i,0] = radius
			gofrset[i,1] = count/(vshell*density)
			gofrset[i,2] = count
			gofrset[i,3] = vshell
			#print gofrset[i,1]
	#results.writerow(gofrset[:,0:1])
	
for i in range(len(gofrset[:])):
	plt(gofrset[i,0:1])
	raw_input()
	

		

print "I'm Done!"
"""
			for m in range(num_part):
				if (current[m,0] - radius < mini[0]) | current[m,0] + radius > maxi[0]:
					continue
				elif (current[m,1] - radius < mini[1]) | current[m,1] + radius > maxi[1]:
					continue
				elif (current[m,2] - radius < mini[2]):
					h = current(m,2) - mini[2]
					vshell += ne.evaluate("1/3.0*np.pi*(2*(radius**3 - inn**3) + 3*h*(radius**2 - inn**2))")
				elif (current[m,2] + radius > maxi[2]):
					h = maxi[2] - current[m,2]
					vshell += ne.evaluate("1/3.0*np.pi*(2*(radius**3 - inn**3) + 3*h*(radius**2 - inn**2))")
				else
					vshell += 4/3.0*np.pi*(radius**3 - inn**3)
"""


