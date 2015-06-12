import numpy as np
import matplotlib.pyplot as plt
"""Note that the above implements the binary model synchronously in it egf and hrg are both set 1
   The indexes 0 to 6 represent EGF,HRG,EGFR,RAF,Pi3k,ERK,AKT"""

def rules( initial_condition ):
	y = np.copy(initial_condition)
	x = initial_condition

	if( x[0] == 1 or x[1] == 1):  #EGFR
		y[2] = 1
	else:
		y[2] = 0

	if( x[2] == 1 or x[6] == 1): #Raf
		y[3] = 1
	else:
		y[3] = 0

	if( x[2] == 1 and x[5] != 1):  #Pi3k
		y[4] = 1
	else:
		y[4] = 0

	if(x[3] == 1):  #ERK
		y[5] = 1
	else:
		y[5] = 0

	if( x[4] == 1):  #Akt
		y[6] = 1
	else:
		y[6] = 0
	return y

def main():
	initial_condition = np.array([1, 1, 0, 0, 0, 0, 0])  #initial value of the species
	t_stop = 10  

	####Setting the initial condition
	y = np.copy(initial_condition)
	step = 1

	##this calculates the value of each species during the time course
	while (step < t_stop + 1):
		if( step == 1):
			temp = rules( y )	
		else:
			temp = rules( y[step-1])
		y = np.vstack((y, temp))
		step += 1
	print y
	###Plot graph
	plt.title("Synch")
	lines = plt.plot(np.arange(t_stop + 1), y[:,5], np.arange(t_stop + 1), y[:,4])	
	plt.setp(lines[0],antialiased = False,color ='#000000',linewidth = 2, marker = "o", markeredgecolor = 'green', label = "erk")
	plt.setp(lines[1],antialiased = False,color ='red',linewidth = 2,linestyle = '--', marker = 'D', markeredgecolor = 'blue', markerfacecolor = 'none', label = "akt")
	plt.legend(loc='upper right')
	plt.xlabel('Time')
	plt.ylabel('Species')
	plt.axis([0,11,-0.05,1.2])
	plt.grid(True)
	plt.show()

if __name__ == "__main__":
	main()


	