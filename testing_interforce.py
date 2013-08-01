import numpy as np 
import matplotlib.pyplot as matplot

a = np.linspace(.1,2**(1.0/6.0)+5, 100)
i=0
forces = np.zeros(shape = (100))
for distance in a:
	beta_epsil = 40
	sigma = 1
	beta = (1.38e-23 * 298.15)**-1
	forces[i] = beta_epsil/beta * -4*(-12*(sigma**12)/distance**6.5 + 6*(sigma**6)/distance**3.5)
	i+=1

print forces

matplot(a,forces)
