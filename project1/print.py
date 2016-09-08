import matplotlib.pyplot as plt
import numpy as np

n = 10
x_axis = np.linspace(0, 1, n)

infile = open('data.dat', 'r')
u_arr = [];
for line in infile:
	values= line.split()
	u_arr.append(float(values[0]))
infile.close()

x = np.linspace(0, 1, 1000)
u_exact = 1 - (1-np.exp(-10))*x - np.exp(-10*x)

plt.plot(x, u_exact, label="u_exact")
plt.plot(x_axis, u_arr, label="u_arr")
plt.legend()
plt.show()


	

