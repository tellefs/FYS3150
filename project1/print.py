import matplotlib.pyplot as plt
import numpy as np


infile = open('data.dat', 'r')
u_arr = [];
n = 0
for line in infile:
	values= line.split()
	u_arr.append(float(values[0]))
	n=n+1
infile.close()

x_axis = np.linspace(0, 1, n)

x = np.linspace(0, 1, n)
u_exact = 1 - (1-np.exp(-10))*x - np.exp(-10*x)
	
plt.plot(x, u_exact, label="u_exact")
plt.plot(x_axis, u_arr, label="u_arr")
plt.legend()
plt.show()


	

