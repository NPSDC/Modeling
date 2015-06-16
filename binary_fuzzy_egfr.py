import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz
from mfs import *

def find_index(arr, value):
	for i in xrange(arr.size):
		if( arr[i] == value):
			return i
	return -1

def find_closest(arr, value):
	dif = abs(arr[0] - value)
	min_dif = dif
	index = 0
	for i in xrange(arr.size):
		dif = abs(arr[i] - value)
		if dif < min_dif:
			min_dif = dif
			index = i
	return index

def compute_egfr_change(egf_value, hrg_value, time_value, initial_values, mfs ):	

	'''Rules---
	If egf is high or hrg is high and time is high then positive_change_egfr is high
	If egf is high or hrg is high and time is low then positive_change_egfr is low
	if egf is low and hrg is low and time is high then negative_change_egfr is high
	if egf is low and hrg is low and time is low then negative_change_egfr is low
	'''
	####Positive CHange#######
	####Antecedent1
	a1_1 = mfs[0][1][initial_values[0] == egf_value] #egf_high[egf == egf_value]
	a1_2 = mfs[1][1][ initial_values[1] == hrg_value ] #hrg_high[hrg == hrg_value]
	a1_3 = mfs[3][1][ initial_values[3] == time_value] #time_high[time == time_value]

	if( a1_1.size == 0):
		a1_1 = mfs[0][1][ find_closest(initial_values[0], egf_value)]
	if( a1_2.size == 0):
		a1_2 = mfs[1][1][ find_closest(initial_values[1], hrg_value)]
	if( a1_3.size == 0):
		a1_3 = mfs[3][1][ find_closest(initial_values[3], time_value)]


	a1 = max( a1_1, a1_2) 
	a1 = min(a1, a1_3 )
	#Consequent 1
	c1 = np.fmin(np.linspace(a1, a1, mfs[2][3].size), mfs[2][3])  #mfs[2][3] is positve_change_egfr_high

	####Antecedent2
	'''Note a2_1 and a2_2 are the same as a1_1 and a1_2 and hence not calculated'''
	a2_3 = mfs[3][0][ initial_values[3] == time_value] #time_low[time == time_value]

	a2 = max( a1_1, a1_2) 
	a2 = min(a2, a2_3 )
	#Consequent 2
	c2 = np.fmin(np.linspace(a2, a2, mfs[2][2].size), mfs[2][2])	#mfs[2][2] is positve_change_egfr_low

	c_com_positive = np.fmax(c1, c2)
	pos_change = fuzz.defuzz(initial_values[2][1], c_com_positive, 'centroid' ) #initial_values[2][1] is positive_change_egfr


	#### Negative Change ######
	####Antecedent 3
	a3_1 = mfs[0][0][initial_values[0] == egf_value] #egf_low[egf == egf_value]
	a3_2 = mfs[1][0][initial_values[1] == hrg_value] #hrg_low[hrg == hrg_value]
	a3_3 = mfs[3][1][initial_values[3] == time_value] #time_high[time == time_value]

	if( a3_1.size == 0):
		a3_1 = mfs[0][0][ find_closest(initial_values[0], egf_value)]
	if( a3_2.size == 0):
		a3_2 = mfs[1][0][ find_closest(initial_values[1], hrg_value)]
	if( a3_3.size == 0):
		a3_3 = mfs[3][0][ find_closest(initial_values[3], time_value)]

	a3 = min(a3_1, a3_2, a3_3)
	#Consequent 3
	c3 = np.fmin(np.linspace(a3, a3, 100), mfs[2][5]) #mfs[2][5] is negative_change_egfr_high

	'''Note a4_1 and a4_2 are the same as a3_1 and a3_2 and hence not calculated'''
	a4_3 = mfs[3][0][initial_values[3] == time_value] #time_low[time == time_value]
	a4 = min(a3_1, a3_2, a4_3)

	#Consequent 4
	c4 = np.fmin(np.linspace(a4, a4, mfs[2][4].size), mfs[2][4]) #mfs[2][4] is negative_change_egfr_low

	c_com_negative = np.fmax(c3, c4)
	neg_change = fuzz.defuzz(initial_values[2][2], c_com_negative, 'centroid') #initial_values[2][2] is negative_change_egfr

	return pos_change + neg_change

def compute_raf_change(egfr_value, akt_value, time_value, initial_values, mfs):
	"""Rules---
	If egfr is high or akt is high and time is high then positive_change_raf is high
	If egfr is high or akt is high and time is low then positive_change_raf is low
	if egfr is low and akt is low  and time is low then negative_change_raf is low
	if egfr is low and akt is low  and time is high then negative_change_raf is high
	"""

	a1_1 = mfs[0][1][initial_values[0] == egfr_value] #egfr_high[egfr == egfr_value]
	a1_2 = mfs[1][1][initial_values[1] == akt_value] #akt_high[akt == akt_value]
	a1_3 = mfs[3][1][initial_values[3] == time_value]

	if( a1_1.size == 0):
		a1_1 = mfs[0][1][ find_closest(initial_values[0], egfr_value)]
	if( a1_2.size == 0):
		a1_2 = mfs[1][1][ find_closest(initial_values[1], akt_value)]
	if(a1_3.size == 0):
		a1_3 = mfs[3][1][ find_closest(initial_values[3], time_value)]

	a1 = max( a1_1 , a1_2 )
	a1 = min(a1_3, a1)
	#print egfr_value
	c1 = np.fmin( np.linspace(a1, a1, 100), mfs[2][1])

	a2_1 = mfs[0][0][initial_values[0] == egfr_value] #egfr_low[egfr == egfr_value]
	a2_2 = mfs[1][0][initial_values[1] == akt_value] #akt_low[akt == akt_value]
	a2_3 = mfs[3][0][initial_values[3] == time_value]

	if( a2_1.size == 0):
		a2_1 = mfs[0][0][ find_closest(initial_values[0], egfr_value)]
	if( a2_2.size == 0):
		a2_2 = mfs[1][0][ find_closest(initial_values[1], akt_value)]
	if( a2_3.size == 0):
		a2_3 = mfs[3][0][ find_closest(initial_values[3], time_value)]

	a2 = min(a2_1 ,a2_2 )
	a2 = max(a2_3, a2)
	c2 = np.fmin( np.linspace(a2, a2, 100), mfs[2][0])

	c_com = np.fmax(c1, c2)
	return fuzz.defuzz(initial_values[2], c_com, 'centroid')

def compute_pi3k_change(egfr_value, erk_value, time_value, initial_values, mfs):
	"""Rules ---
	if egfr is high and erk is low and time is high then pi3k is high
	if egfr is low or erk is high or time is low then pi3k is low"""

	a1_1 = mfs[0][1][initial_values[0] == egfr_value] #egfr_high[egfr == egfr_value]
	a1_2 = mfs[1][0][ initial_values[1] == erk_value] #erk_low[erk == erk_value]
	a1_3 = mfs[3][1][ initial_values[3] == time_value]

	if( a1_1.size == 0):
		a1_1 = mfs[0][1][ find_closest(initial_values[0], egfr_value)]
	if( a1_2.size == 0):
		a1_2 = mfs[1][0][ find_closest(initial_values[1], erk_value)]
	if( a1_3.size == 0):
		a1_3 = mfs[3][1][ find_closest(initial_values[3], time_value)]

	a1 = min(a1_1 , a1_2, a1_3)
	c1 = np.fmin( np.linspace(a1, a1, 100), mfs[2][1])

	a2_1 = mfs[0][0][ initial_values[0] == egfr_value] #egfr_low[egfr == egfr_value]
	a2_2 = mfs[1][1][ initial_values[1] == erk_value] #erk_high[erk == erk_value]
	a2_3 = mfs[3][0][ initial_values[3] == time_value]

	if( a2_1.size == 0):
		a2_1 = mfs[0][0][ find_closest(initial_values[0], egfr_value)]
	if( a2_2.size == 0):
		a2_2 = mfs[1][1][ find_closest(initial_values[1], erk_value)]
	if( a2_3.size == 0):
		a2_3 = mfs[3][0][ find_closest(initial_values[3], time_value)]

	a2 = max(a2_1 , a2_2, a2_3)
	c2 = np.fmin( np.linspace(a2, a2, 100), mfs[2][0] )

	c_com = np.fmax(c1, c2)
	return fuzz.defuzz(initial_values[2], c_com, 'centroid')

def compute_erk_change(raf_value, time_value, initial_values, mfs):
	"""Rules-
		If raf is high and time is high erk is high
		If raf is low or time is low then erk is low"""
	a1_1 = mfs[0][1][initial_values[0] == raf_value]
	a1_2 = mfs[2][1][initial_values[2] == time_value]

	if( a1_1.size == 0):
		a1_1 = mfs[0][1][ find_closest(initial_values[0], raf_value)]
	if( a1_2.size == 0):
		a1_2 = mfs[2][1][ find_closest(initial_values[2], time_value)]

	a1 = min(a1_1, a1_2)
	c1 = np.fmin( np.linspace(a1, a1, 100), mfs[1][1])

	a2_1 = mfs[0][0][initial_values[0] == raf_value]
	a2_2 = mfs[2][0][initial_values[2] == time_value]

	if( a2_1.size == 0):
		a2_1 = mfs[0][0][ find_closest(initial_values[0], raf_value)]
	if( a2_2.size == 0):
		a2_2 = mfs[2][0][ find_closest(initial_values[2], time_value)]

	a2 = max(a2_1, a2_2)
	c2 = np.fmin( np.linspace(a2, a2, 100), mfs[1][0])

	c_com = np.fmax(c1,c2)
	return fuzz.defuzz( initial_values[1], c_com, 'centroid')

def compute_akt_change(pi3k_value, time_value, initial_values, mfs):
	"""Rules-
		If pi3k is high and time is high akt is high
		If pi3k is low or time is low then akt is low"""
	a1_1 = mfs[0][1][initial_values[0] == pi3k_value]
	a1_2 = mfs[2][1][initial_values[2] == time_value]

	if( a1_1.size == 0):
		a1_1 = mfs[0][1][ find_closest(initial_values[0], pi3k_value)]
	if( a1_2.size == 0):
		a1_2 = mfs[2][1][ find_closest(initial_values[2], time_value)]

	a1 = min(a1_1, a1_2)
	c1 = np.fmin( np.linspace(a1, a1, 100), mfs[1][1])

	a2_1 = mfs[0][0][initial_values[0] == pi3k_value]
	a2_2 = mfs[2][0][initial_values[2] == time_value]

	if( a2_1.size == 0):
		a2_1 = mfs[0][0][ find_closest(initial_values[0], pi3k_value)]
	if( a2_2.size == 0):
		a2_2 = mfs[2][0][ find_closest(initial_values[2], time_value)]

	a2 = max(a2_1, a2_2)
	c2 = np.fmin( np.linspace(a2, a2, 100), mfs[1][0])

	c_com = np.fmax(c1,c2)
	return fuzz.defuzz( initial_values[1], c_com, 'centroid')

def compute_egfr(change_egfr_reflected, not_updated, initial_cond, time_egfr_index, initial_values, mfs):

	flag = 0
	egfr = initial_cond[2]

	if(time_egfr_index > 100):
		flag = 1
		change_egfr_reflected = egfr
	
	if 0 in not_updated or 1 in not_updated:		
		time_egfr_index = time_egfr_index + 1		
	else:
		flag = 0
		time_egfr_index = 2

	if flag == 0:
		egfr = change_egfr_reflected + compute_egfr_change(initial_cond[0], initial_cond[1], initial_values[7][ time_egfr_index - 1], \
			   (initial_values[0], initial_values[1], initial_values[2], initial_values[7]), (mfs[0], mfs[1], mfs[2], mfs[7]))
	
	return (egfr, time_egfr_index, change_egfr_reflected)

def check_raf(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs):

	if 2 in not_updated and 6 in not_updated:
		time_indexes[1] = time_indexes[1] + 1	
	elif 2 in not_updated and initial_cond[2] > 0.8:
		time_indexes[1] = time_indexes[1] + 1
	elif 6 in not_updated and initial_cond[6] > 0.8:
		time_indexes[1] = time_indexes[1] + 1
	else:
		time_indexes[1] = 2

	y[3] = compute_raf(initial_cond[2], initial_cond[6], initial_values[7][ time_indexes[1] - 1], \
	 (initial_values[2], initial_values[6], initial_values[3], initial_values[7]), (mfs[2], mfs[6], mfs[3], mfs[7]))

	return (y, time_indexes)

def check_pi3k(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs ):

	if 2 in not_updated and 5 in not_updated:		
		time_indexes[2] = time_indexes[2] + 1
		#if time_indexes[0] < 300:
		#	print time_indexes[0]
	else:
		time_indexes[2] =  2

	y[4] = compute_pi3k(initial_cond[2], initial_cond[5], initial_values[7][ time_indexes[2] - 1], \
		 (initial_values[2], initial_values[5], initial_values[4], initial_values[7]), (mfs[2], mfs[5], mfs[4], mfs[7]))
			
	return (y, time_indexes)

def check_erk(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs ):

	if 3 in not_updated:		
		time_indexes[3] = time_indexes[3] + 1
	else:
		time_indexes[3] = 2


	y[5] = compute_erk(initial_cond[3], initial_values[7][ time_indexes[3] - 1], \
		 (initial_values[3], initial_values[5], initial_values[7]), (mfs[3], mfs[5], mfs[7]))	

	
	return (y, time_indexes)

def check_akt(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs ):

	if 4 in not_updated:
		y[6] = compute_akt(initial_cond[4], initial_values[7][ time_indexes[4]], \
		 (initial_values[4], initial_values[6], initial_values[7]), (mfs[4], mfs[6], mfs[7]))
		time_indexes[4] = time_indexes[4] + 1

	else:
		time_indexes[4] = 2	
	return (y, time_indexes)


def rules(change_reflected, prev_cond, initial_cond, time_indexes, (initial_values, mfs)):	
	y = np.copy(initial_cond)
	not_updated = []  #contains the indexes not updated
	for i in xrange(len(prev_cond)):
		if str(prev_cond[i]) == str(initial_cond[i]):
			not_updated.append(i)	
		
	y[2], time_indexes[0], change_reflected[2] = compute_egfr(change_reflected[2], not_updated, initial_cond, time_indexes[0], initial_values, mfs )
	#y, time_indexes = check_raf(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs )
	#y, time_indexes = check_pi3k(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs )
	#y, time_indexes = check_erk(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs )
	#y, time_indexes = check_akt(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs )
	#y[5] = initial_cond[3]
	#y[6] = initial_cond[4]
	#y = np.around(y, 8)
	return (y,time_indexes, change_reflected)

def main():
	
	initial_cond = np.array([1, 1, 0, 0, 0, 0, 0], dtype = "float64")
	time_stop = 10
	y = np.copy(initial_cond)
	change_reflected = np.copy(initial_cond) # This is the array to which we would adding the changing values
	y.resize(1, 7)
	step = 1

	egf = np.linspace(0, 1, 100)
	hrg = egf
	egfr = egf
	akt = egf
	raf = egf
	pi3k = egf
	erk = egf

	positive_change_egfr = egf
	positive_change_raf = egf
	positive_change_pi3k = egf
	positive_change_erk = egf
	positive_change_akt = egf

	negative_change_egfr = np.linspace(-1, 0, 100)
	negative_change_raf = negative_change_egfr
	negative_change_pi3k = negative_change_egfr
	negative_change_erk = negative_change_egfr
	negative_change_akt = negative_change_egfr

	time = np.linspace(0, 10, 1001)
	vals = (egf, hrg, (egfr, positive_change_egfr, negative_change_egfr),\
		   (raf, positive_change_raf, negative_change_raf),\
		   (pi3k, positive_change_pi3k, negative_change_pi3k),\
		   (erk, positive_change_erk, negative_change_erk),\
		   (akt, positive_change_akt, negative_change_akt), time)

	mfs = eval_membership_functions(vals)

	times = [ 1, 1, 1, 1, 1] #provides the initial time inputs (indexes 0 for egfr, 1 for raf, 2 for pi3k and so on)
	for i in xrange(1, time.size ):
			if(i > 100 and i %100 == 1):
				y[i - 1][0] = int(y[i - 1][0]) ^ 1
				y[i - 1][1] = int(y[i - 1][1]) ^ 1
			temp, times, change_reflected = rules(change_reflected, y[i-2], y[i - 1], times , (vals, mfs))
			y = np.vstack((y, temp))
			if i < 201:
				print change_reflected, times[0], y[i]
	#print y[98],y[99],y[100],y[101]
	
	plt.title("Synch")
	
	lines = plt.plot(time, y[:,2]	)	
	
	plt.legend(loc='upper right')
	plt.xlabel('Time')
	plt.ylabel('Species')
	plt.axis([-0.2,10.1,-0.05,1.2])
	plt.grid(True)
	plt.show()

if __name__ == "__main__":
	main()