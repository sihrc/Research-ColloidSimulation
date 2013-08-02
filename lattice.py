#Generating Lattice
#Bravais Vector for FCC 

import numpy as np 
import itertools as it
import csv,sys,math

#"""
#Getting Parameters
num_part = int(sys.argv[1])
sigma = float(sys.argv[2])
volume_fraction = float(sys.argv[3])
print sys.argv
#Opening the lattice file
datafile = open('lattice'+sys.argv[1]+"_"+str(volume_fraction*100)+'.txt','wb')
data = csv.writer(datafile)

#Uni-axial Distance between particles
a = (14*4.0/3.0*math.pi*sigma**3/volume_fraction)**(1.0/3.0)/2

#Calculating size of combinations
size = 0
while 1:
	size += 1
	if (size - 2)*(size-1)*size > 6*num_part:
		break
#Setting up Bravais Vectors
v1 = np.array([0,a,a ])
v2 = np.array([a,0,a])
v3 = np.array([a,a,0])

#Setting up Combinations
ints = np.linspace(-size/2,size/2,size+1)
combinations = list(it.permutations(ints,3))
#print combinations

#Initializing particles
particles = np.zeros(shape = (len(combinations),4))

#Iterating through combinations
i = 0
for comb in combinations:
	particles[i,0:3] = comb[0]*v1 + comb[1] * v2 + comb[2] * v3
	particles[i,3] = np.sum(particles[i,0:3]**2)**.5
	i += 1

#Sorting particles based on distance from origin
particles.view('f8,f8,f8,f8').sort(order=['f3'],axis = 0)

#Writing positions to file
for line in particles:
	data.writerow(line)

