"""
Brownian Dynamics Simulation: simulating N suspensed particles
Chris Lee
Created: October 12, 2012
Modified: December 14

Usage: python msd.py [datafile.txt]
"""

import matplotlib.pyplot as plt
import sys
import csv
import numpy as np

"""Importing the data"""
datafile = open(sys.argv[1],'r')

data = csv.reader(datafile)
"""Reading the parameter data"""
info = data.next()
num_part = int(info[0])
total_lines = float(info[1])/float(info[2])*num_part

"""Extracting Data for MSD calculations"""
#Initializing data structures
n = 0
particles = np.zeros(shape=(num_part,3))
averages= []
values = [0,0,0]
value = np.zeros(num_part)
#Running through the file
while data.line_num < total_lines + 1:
        n += 1
        a = data.next()
       	values = np.array([float(a[0]),float(a[1]), float(a[2])])
        particles[n-1] += values
        value[n-1] = (particles[n-1].sum())**2
        if n == num_part:
            averages.append(value.sum()/num_part)
        n = n%num_part

"""Plots the Mean Square Displacement for all the particles"""
plt.plot(averages)
plt.xlabel('Step Iteration')
plt.ylabel('MSD (mme )')
plt.show()