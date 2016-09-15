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


infile = open('special_data.dat', 'r')
u_arr_special = [];
for line in infile:
	values= line.split()
	u_arr_special.append(float(values[0]))
infile.close()

infile = open('error_data.dat', 'r')
error = [];
for line in infile:
	values= line.split()
	error.append(float(values[0]))
infile.close()
	


plt.plot(x, u_exact, label="u_exact")
plt.plot(x_axis, u_arr, label="u_arr")
plt.plot(x_axis, u_arr_special, label='u_arr_special')
plt.xlabel('x')
plt.ylabel('U(x)')
plt.legend()
plt.title('Plot for n = %i' %(n-2))
plt.show()


