"""
Brownian Dynamics Simulation: Simulating Crystallization of almost hard-sphere colloids
Chris Lee
May 6, 2013

Citation: "Simulation of nucleation in almost hard-sphere colloids" - Filion, Ni, Frenkel, Dijkstra
Usage: python main.py [number of particles][duration][time step][volume fraction]
"""

#Dependencies - Located in the "Install" folder
import numpy as np
import numexpr as ne
import itertools as it
import csv, time, sys, os, math
from subprocess import call

def generate(n):
    """Generates n number of particles data in a matrix with n rows"""
    #Initiate the Matrix of particles
    particles = np.zeros(shape = (n,10)) 
    #Calculates the size of the container based on volume fraction
    side = ((n*4.0/3.0*math.pi*(sigma/2.0)**3/volume_fraction))**(1.0/3.0)

    #""" Add/Remove the # to (un)comment entire section
    #Filling in random positions for the particles
    particles[:,0:3] = np.random.uniform(-.5*side, .5*side,size = (n,3))
    return particles
    #"""

    """ Add/Remove the # to (un)comment entire section
    #Generate Determined Lattice
    read = lattice()
    particles[:,0:3] = read[0:n,0:3]
    print "Done importing lattice"
    return particles
    #"""
    
def no_overlap(positions):
    """Runs the initial particles in a simulation to sort out overlaps :: Returns new positions"""
    #Setting up distance matrix
    x = positions[:,0] - positions[:,0][:,np.newaxis]
    y = positions[:,1] - positions[:,1][:,np.newaxis]
    z = positions[:,2] - positions[:,2][:,np.newaxis] 
    distance = (x**2 + y**2 + z**2)**.5
    #Filling in for (self - self) particles
    np.fill_diagonal(distance, 2*sigma)

    #Searching for overlap cases
    nearenough = distance<sigma
    np.fill_diagonal(nearenough,False)

    #Moving the particles from overlapping
    positions[:,0] += np.sum(nearenough*-sigma/10.0*x/distance ,axis = 1)
    positions[:,1] += np.sum(nearenough*-sigma/10.0*y/distance ,axis = 1)
    positions[:,2] += np.sum(nearenough*-sigma/10.0*z/distance ,axis = 1)
    return positions

def overlaps(positions):
    """Checks if the particles are still overlapping :: Returns False when no longer overlapping"""
    distance = np.sum((positions - positions[:,np.newaxis])**2,axis = 2)**.5 < sigma
    np.fill_diagonal(distance,False)
    if distance.any():
        return True
    else:
        return False


def lattice():
    """Generates a Bravais Lattice for positions of particles :: Returns a matrix of positions"""
    #Uni-axial Distance between particles
    a = .5*(16.0/3.0*math.pi*(sigma/2.0)**3/volume_fraction)**(1.0/3.0)
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
    combinations = [p for p in it.product(ints,repeat = 3)]
    #Initializing particles
    bravais_lattice = np.zeros(shape = (len(combinations),4))

    #Iterating through combinations
    i = 0
    for comb in combinations:
        bravais_lattice[i,0:3] = comb[0]*v1 + comb[1] * v2 + comb[2] * v3
        bravais_lattice[i,3] = np.sum(bravais_lattice[i,0:3]**2)**.5
        i += 1
    #Sorting particles based on distance from origin
    bravais_lattice.view('f8,f8,f8,f8').sort(order=['f3'],axis = 0)
    return bravais_lattice

def random_force():
    """Calculates the random force on particles in each component [x,y,z]
    using a Gaussian Random Distribution"""

    #Initializing matrix
    total = np.zeros(shape=(num_part,3))

    #Averaging to get rid of extremeties
    for i in range(10):
        total += np.random.normal((6.*mass*gamma*kbT)**.5/2, size=(num_part,3))

    particles[:,3:6] +=  random_coeff * ne.evaluate("total/10.0")

def wall_check():
    """Implementing Periodic Boundary Conditions. Interaction with particles beyond the boundaries"""
    #Finding particles near the boundaries
    under = np.where(particles[:,0:3] < -side/2 + inter_dist)[0]
    over = np.where(particles[:,0:3] > side/2 - inter_dist)[0]

    #Finds particles near enough to boundaries
    set1 = particles[:,0:3][under]

    #Offset over particles to set them close enough to interact with under
    set2 = particles[:,0:3][over] - side 

    #Calculate inter-distances between over and under
    distance = set2 - set1[:,np.newaxis]
    #Finds distances that are close enough
    near = distance < inter_dist
    if -near.all():
        return
    #Calculating the forces
    forces = ne.evaluate("near*beta_epsil/beta * -4*(-12*(sigma**12)/distance**13 + 6*(sigma**6)/distance**7)")

    #Updating
    particles[:,3:6][under] += -np.sum(forces,axis = 1)
    particles[:,3:6][over] += -np.sum(forces,axis = 0)

def inter_force():
    """Calculates all interaction potentials within the volume and returns the velocity it causes
    Equations adopted from simulatingon of Nucleation in almost hard-sphere colloids. Journal of Chemical Physics 134, 134901 (2011)"""
    x = particles[:,0] - particles[:,0][:,np.newaxis]
    y = particles[:,1] - particles[:,1][:,np.newaxis]
    z = particles[:,2] - particles[:,2][:,np.newaxis] 

    #Calculates squared distance matrix
    distance = ne.evaluate("x**2+y**2+z**2")

    #Getting rid of 0 distances
    np.fill_diagonal(distance, 1)

    #Indicator
    particles[:,9] = 0

    #Find particles close enough for particle interaction
    near = (distance < inter_dist**2)
    
    #Ignoring self-self distance
    np.fill_diagonal(near,False)

    #Returning if no particles are close enough
    if near.any():
        forces = ne.evaluate("inter_coeff*near*beta_epsil/beta * 4*(12*(sigma**12)/distance**6.5 - 6*(sigma**6)/distance**3.5)")
        particles[np.where(near)[0],9] = 1
    else:
        forces = 0
    
    #Assign velocities to particles
    particles[:,3] = np.sum(ne.evaluate("x/distance**.5*forces"),axis = 1)
    particles[:,4] = np.sum(ne.evaluate("y/distance**.5*forces"),axis = 1)
    particles[:,5] = np.sum(ne.evaluate("z/distance**.5*forces"),axis = 1)
    #print "inter",particles[:,3:5]
    #raw_input()


def Eulers():
    """Calculating new positions using Forward Euler Method: Updates Particles"""
    #Calculating Interaction force
    inter_force()

    #Checks acceleration due to periodic wall boundaries
    wall_check()

    #Random Force - Van der Waal
    random_force()

    #Euler's Method 
    particles[:,3:6] += particles[:,6:9]*dt/(mass*gamma)
    particles[:,0:3] += particles[:,3:6]*dt

    #Periodic Boundary Conditions
    over = particles[:,0:3]>side/2
    under = particles[:,0:3]<-side/2
    particles[:,0:3] +=  over*-side + under*side


        
##################################################################################
##########################     Main Simulation    ################################
##################################################################################
"""Parameters"""
global V,side,sigma,inter_str,num_part, mass, gamma, dt, particles, data, kbT, inter_coeff, inter_dist, orig, k, volume_fraction

"""Experiment Inputs"""
num_part=int(sys.argv[1])                       #Number of Particles
dt= float(sys.argv[3])                          #Time Step
duration=float(sys.argv[2])                     #Duration of Simulation in minute-s
volume_fraction = float(sys.argv[4])            #Starting Volume Fraction

#Container and particles
sigma = 1                                       #Diameter of Particles (m)
mass = 1                                        #Mass of each particle (kg?)
gamma = 1                                       #Friction coefficient
pi = np.pi                                      #pi!

#Parameters provided by cited paper
boltz = 1.38e-34                       #Boltzman Constant
temp = 298                             #Temperature (K)
kbT = boltz * temp                     #Energy
beta = kbT**-1                         #Energy constant in potential equation
beta_epsil = 40.0                      #Energy constant with potential well
epsil = beta_epsil/beta                #Potential well
inter_dist = sigma*(2.0**(1.0/6.0))    #Interaction potential cut-off

#Calibration Constants
random_coeff = 1                       #Calibration constant for random force
inter_coeff = 1e22                     #Calibration constant for interaction force

"""Writing Data to File"""
#Setting up: Write data to file
filename = str(num_part) + "_" + str(volume_fraction) + "_" + str(duration) 
if (filename + ".txt") in os.listdir("."):
    ext = raw_input("Do you want to add an extension to your filename?")
    if ext == "n":
        filename = filename + ".txt"
        os.remove(filename)
    else:
        filename = filename + ext + ".txt"
else:
    filename = filename + ".txt"

frames = 1/dt                        #Framerate recorded
datafile = open(filename,'wb')
data = csv.writer(datafile)

#MSD Validation
#msdfile = open("msdfile.txt","wb")
#msd = csv.writer(msdfile)
#msd.writerow([num_part, duration, dt, sigma/2., side])
                                            
"""Simulatation"""
#Timing the simulation
start_time = time.clock()

#Creating the Particles
particles = generate(num_part)
#Sizes the Container Based on particles
side = 2.5*float(np.amax(abs(particles[:,0:3])))
#Settling overlaps in the initialization
while overlaps(particles[:,0:3]):
    particles[:,0:3] = no_overlap(particles[:,0:3])
#Writing initial States
data.writerow([num_part, duration, dt, sigma/2.0, side, frames])

#Printing the simulation properties
print "Number of Particles", num_part
print "duration", duration
print "time step", dt
print "Side of container", side

#Recording initial positions
orig = np.array(particles[0,0:3])

#Simulation timer/counter
elapsed_time = 0
i = 0

#Running the Simulation
while elapsed_time<duration:
    elapsed_time += dt
    start = time.clock()
    Eulers()
    for part in range(num_part):
        if (i%(frames)) == 1:
            data.writerow(particles[part,[0,1,2,9]])
        #Outputs [dx,dy,dz]
        #msd.writerow(particles[part,6:9]*dt)
    if i%(duration/(100*dt)) == 1:
        print 100*elapsed_time/duration, "%", " completed"
        print time.clock()/(elapsed_time/duration) - time.clock(), "seconds remaining"
    i += 1

#Runtime
print time.clock() - start_time
data.writerow('Done1')

#Waits for User Input
raw_input()

#Begins the Animation in visual_main.py <- must be in the same directory.
command = "python visual_main.py " + filename + " " + str(num_part)
call(command)
print "Done"