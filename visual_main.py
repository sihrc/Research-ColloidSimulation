"""
Brownian Dynamics Simulation: simulating N suspensed particles
Chris Lee
October 12, 2012

Usage: python visual_main.py [datafile.txt][max loops per sec]
"""

from visual import *
import os
import sys
import csv
import random


"""Importing the data"""
datafile = open(sys.argv[1],'r')

data = csv.reader(datafile)
"""Reading the parameter data"""
info = data.next()
num_part = int(info[0])
total_lines = float(info[1])/float(info[2]) * num_part
radius = float(info[3])
side = float(info[4])
"""Initializaing dictionary and iterator"""
particles = []
n = 0
"""Setting up window"""
window = display(title = "Colloidal Particles Simulation" , height = 900, width = 1800, center=(0,0,0))
window.autoscale = True
"""Making the container"""
#side = side * 10
print side
a = side/2.0
for i in [1,2,3]:
        print a
        curve(pos = [(a,a,a),(a+(i==1)*-side,a+(i==2)*-side,a+(i==3)*-side)], color = color.blue)
        curve(pos = [(-a,-a,a),(-a+(i==1)*side,-a+(i==2)*side,a+(i==3)*-side)], color = color.blue)
        curve(pos = [(-a,a,-a),(-a+(i==1)*side,a+(i==2)*-side,-a+(i==3)*side)], color = color.blue)
        curve(pos = [(a,-a,-a),(a+(i==1)*-side,-a+(i==2)*side,-a+(i==3)*side)], color = color.blue)
        
"""Extracting Data into sphere objects"""
#Initializing spheres
for i in range(num_part):
        #RGB = (random.random(),random.random(),random.random())
        RGB = (1,1,1)
        box = sphere(pos=(0,0,0), radius = radius, color = RGB)
        particles.append(box)
        
#Running through the datafile
while 1:
        coord = data.next()
        if len(coord) == 5:
                break
        ball = particles[n]     
        ball.pos = vector(float(coord[0]),float(coord[1]),float(coord[2]))
        if float(coord[3]) == 1:
                ball.color = (1,0,0)
        else:
                ball.color = RGB
        n += 1
        print ball.pos
        if n == num_part:
                n = 0
        #Controls the # of times the loop runs per second
        rate(float(sys.argv[2]))

