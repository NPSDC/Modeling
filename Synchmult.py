import numpy as np
import matplotlib.pyplot as plt

"""Since egfr is more sensitive to egf than hrg we start with taking egf as 0 and then plot
   indexes 0 to 6 are for egf,hrg,egfr,raf,pi3k,erk,akt"""		

def main():
	initial_condition = np.array([0,1,0,0,0,0,0])  # stores the initial values of the various components at t=0
	t_stop = 20

	y = np.copy(initial_condition)
	step =  1

	while( step < t_stop + 1):
		if(step == 1):
			new_order = rules(y)
		else:
			new_order = rules(y[step-1])
		
		y = np.vstack((y,new_order))
		step += 1

	print y

	plt.title("SynchMulti")
	lines = plt.plot(np.arange(t_stop + 1), y[:,5], np.arange(t_stop + 1), y[:,6])	
	plt.setp(lines[0],antialiased = False,color ='#000000',linewidth = 2, marker = "o", markeredgecolor = 'green', label ="erg")
	plt.setp(lines[1],antialiased = False,color ='red',linewidth = 2,linestyle = '--',marker = "D", markerfacecolor = "none", label = "akt" )
	plt.legend(loc='upper right')
	ax = plt.gca()
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	plt.axis([0,22,-0.1,1.2])
	plt.xlabel('Time')
	plt.ylabel('Species')	
	plt.show()	

def rules(y_old):
	y = np.copy(y_old)

	#egfr
	if((y_old[0] + 0.5*y_old[1]) >= 1):
		y[2] = 2
	elif((y_old[0] + 0.5*y_old[1]) > 0):
		y[2] = 1
	else:
		y[2] = 0

	#raf
	if((y_old[2] + 2*y_old[6]) > 1):
		y[3] = 1
	else:
		y[3] = 0

	#pi3k
	if((y_old[2] - y_old[5]) > 0):
		y[4] = 1
	else:
		y[4] = 0

	#akt
	if( y_old[3] == 1):
		y[5] = 1
	else:
		y[5] = 0

	#erk
	if(y_old[4] == 1):
		y[6] = 1
	else:
		y[6] = 0

	return y

if __name__ == "__main__":
	main()