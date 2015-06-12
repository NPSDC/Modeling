import numpy as np
import itertools
import matplotlib.pyplot as plt
from random import randint

"""Since egfr is more sensitive to egf than hrg we start with taking egf as 0 and then plot
   indexes 0 to 6 are for egf,hrg,egfr,raf,pi3k,erk,akt"""



def main():
	initial_condition = np.array([0,1,0,0,0,0,0])  # stores the initial values of the various components at t=0
	t_stop = 20
	lis = [1,2,3,4,5]
	a =  list(itertools.permutations(lis))
	poss_perms = np.asarray(a)

	y = np.copy(initial_condition)
	step =  1

	while( step < t_stop + 1):
		size,breadth = poss_perms.shape
		order = poss_perms[randint(0,size - 1),:]

		if(step == 1):
			new_order = rules(order, y, breadth)
		else:
			new_order = rules(order,y[step-1 ], breadth)

		
		y = np.vstack((y,new_order))
		print y[4],times[4]
		step += 1



	plt.title("ASynch")
	lines = plt.plot(np.arange(t_stop + 1), y[:,5],np.arange(t_stop + 1), y[:,6])	
	plt.setp(lines[0],antialiased = False,color ='#000000',linewidth = 2, marker = "o", markeredgecolor = 'green', label ="erg")
	plt.setp(lines[1],antialiased = False,color ='red',linewidth = 2,linestyle = '--', marker = "D", markerfacecolor = "none",label ="akt" )
	plt.legend(loc='upper right')
	plt.axis([0,21,-0.05,1.2])
	plt.xlabel('Time')
	plt.ylabel('Species')	
	plt.show()	
	
def rules(new_order,y_old,size):
	y = np.copy(y_old)
	for i in xrange(size):
		new = changeValues(y, new_order[i])
		y = new
	return y

def changeValues(y_old,choice):
	y = np.copy(y_old)

	if(choice == 1):  #egfr
		if((y_old[0] + 0.5*y_old[1]) >= 1):
			y[2] = 2
		elif((y_old[0] + 0.5*y_old[1]) > 0):
			y[2] = 1
		else:
			y[2] = 0

	elif(choice == 2): #raf
		if((y_old[2] + 2*y_old[6]) > 1):
			y[3] = 1
		else:
			y[3] = 0

	elif(choice == 3): #pi3k
		if((y_old[2] - y_old[5]) > 0):
			y[4] = 1
		else:
			y[4] = 0

	elif(choice == 4): #erk
		if( y_old[3] == 1):
			y[5] = 1
		else:
			y[5] = 0

	elif(choice == 5):
		if(y_old[4] == 1):
			y[6] = 1
		else:
			y[6] = 0

	return y

if __name__ == "__main__":
	main()