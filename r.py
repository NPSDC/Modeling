import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz

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

def eval_membership_functions(initial_values):
	'''Technically they should have been calcultated using fuzz.someMf but since they are same to save computation we do this'''
	egf_high = fuzz.gaussmf(initial_values[0], 1, 0.1)
	egf_low = fuzz.gaussmf(initial_values[0], 0, 0.1)
	hrg_high = egf_high
	hrg_low = egf_low
	egfr_high = egf_high
	egfr_low =  egf_low
	erk_high = egfr_high  
	erk_low = egfr_low 
	pi3k_high = egfr_high
	pi3k_low = egfr_low
	akt_high = egfr_high  
	akt_low = egfr_low 
	raf_high = egfr_high
	raf_low = egfr_low
	time_high = fuzz.smf(initial_values[7], 0, 1)
	time_low = fuzz.zmf(initial_values[7], 0, 1)
	return ((egf_low, egf_high), (hrg_low, hrg_high), (egfr_low, egfr_high), (raf_low, raf_high), (pi3k_low, pi3k_high), (erk_low, erk_high), \
			(akt_low, akt_high), (time_low, time_high))


def compute_egfr(egf_value, hrg_value, time_value, initial_values, mfs ):	

	"""Rules---
	If egf is high or hrg is high and time is high then egfr is high
	if egf is low and hrg is low or time is low then egfr is low"""

	a1_1 = mfs[0][1][initial_values[0] == egf_value] #egf_high[egf == egf_value]
	a1_2 = mfs[1][1][ initial_values[1] == hrg_value ] #hrg_high[hrg == hrg_value]
	a1_3 = mfs[3][1][ initial_values[3] == time_value]

	if( a1_1.size == 0):
		a1_1 = mfs[0][1][ find_closest(initial_values[0], egf_value)]
	if( a1_2.size == 0):
		a1_2 = mfs[1][1][ find_closest(initial_values[1], hrg_value)]
	if( a1_3.size == 0):
		a1_3 = mfs[3][1][ find_closest(initial_values[3], time_value)]

	a1 = max( a1_1, a1_2) 
	a1 = min(a1, a1_3 )
	c1 = np.fmin(np.linspace(a1, a1, 100), mfs[2][1])

	a2_1 = mfs[0][0][initial_values[0] == egf_value] #egf_low[egf == egf_value]
	a2_2 = mfs[1][0][initial_values[1] == hrg_value] #hrg_low[hrg == hrg_value]
	a2_3 = mfs[3][0][initial_values[3] == time_value]

	if( a2_1.size == 0):
		a2_1 = mfs[0][0][ find_closest(initial_values[0], egf_value)]
	if( a2_2.size == 0):
		a2_2 = mfs[1][0][ find_closest(initial_values[1], hrg_value)]
	if( a2_3.size == 0):
		a2_3 = mfs[3][0][ find_closest(initial_values[3], hrg_value)]

	a2 = min(a2_1, a2_2)
	a2 = max(a2, a2_3)
	c2 = np.fmin(np.linspace(a2,a2,100), mfs[2][0])

	c_com = np.fmax(c1, c2)
	
	a =  fuzz.defuzz(initial_values[2], c_com, 'centroid')
	return a

def compute_raf(egfr_value, akt_value, time_value, initial_values, mfs):
	"""Rules---
	If egfr is high or akt is high and time is high then raf is high
	if egfr is low and akt is low  or time is low then raf is low"""	

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

def compute_pi3k(egfr_value, erk_value, time_value, initial_values, mfs):
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

def compute_erk(raf_value, time_value, initial_values, mfs):
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

def compute_akt(pi3k_value, time_value, initial_values, mfs):
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

def check_egfr(not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs):
	y = np.copy(initial_cond)
	
	if 0 in not_updated and 1 in not_updated:
		
		y[2] = compute_egfr(initial_cond[0], initial_cond[1], initial_values[7][ time_indexes[0] ], (initial_values[0], initial_values[1], initial_values[2],\
		 initial_values[7]), (mfs[0], mfs[1], mfs[2], mfs[7]))
		time_indexes[0] = time_indexes[0] + 1

	else:
		time_indexes[0] = 1


	

	return (y, time_indexes)

def check_raf(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs):

	if 2 in not_updated or 6 in not_updated:
		time_indexes[1] = time_indexes[1] + 1
		
	else:
		time_indexes[1] = 2

	y[3] = compute_raf(initial_cond[2], initial_cond[6], initial_values[7][ time_indexes[1] - 1], \
		 (initial_values[2], initial_values[6], initial_values[3], initial_values[7]), (mfs[2], mfs[6], mfs[3], mfs[7]))
	

	
	return (y, time_indexes)

def check_pi3k(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs ):

	if 2 in not_updated and 5 in not_updated:
		time_indexes[2] = time_indexes[2] + 1

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


	y[5] = compute_erk(initial_cond[3], initial_values[7][ time_indexes[3]-1], \
		 (initial_values[3], initial_values[5], initial_values[7]), (mfs[3], mfs[5], mfs[7]))	

	
	return (y, time_indexes)

def check_akt(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs ):

	if 4 in not_updated:
		y[6] = compute_akt(initial_cond[4], initial_values[7][ time_indexes[4]], \
		 (initial_values[4], initial_values[6], initial_values[7]), (mfs[4], mfs[6], mfs[7]))
		time_indexes[4] = time_indexes[4] + 1

	else:
		time_indexes[4] = 1


	

	
	return (y, time_indexes)

def rules(prev_cond, initial_cond, time_indexes, (initial_values, mfs)):	
	not_updated = []
	for i in xrange(len(prev_cond)):
		if prev_cond[i] == initial_cond[i]:
			not_updated.append(i)		

	y, time_indexes = check_egfr(not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs )
	y, time_indexes = check_raf(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs )
	y, time_indexes = check_pi3k(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs )
	y, time_indexes = check_erk(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs )
	y, time_indexes = check_akt(y, not_updated, prev_cond, initial_cond, time_indexes, initial_values, mfs )
	#y[5] = initial_cond[3]
	#y[6] = initial_cond[4]
	return (y,time_indexes)

def main():
	#Universal sets
	#c_egfr_aggregated = compute_initial_egfr(egf, hrg, egfr)
	#c_raf_aggregated = compute_initial_raf(egfr, akt, raf)
	#c_pi3k_aggregated = compute_initial_pi3k(egfr, erk, pi3k)

	initial_cond = np.array([1, 1, 0, 0, 0, 0, 0], dtype = "float64")
	time_stop = 10
	y = np.copy(initial_cond)
	y.resize(1, 7)
	step = 1

	egf = np.linspace(0, 1, 100)
	hrg = egf
	egfr = egf
	akt = egf
	raf = egf
	pi3k = egf
	erk = egf
	time = np.linspace(0, 10, 1000)
	vals = (egf, hrg, egfr, raf, pi3k, erk, akt, time)
	mfs = eval_membership_functions(vals)
	
	#print rules(y[0], 0, (vals,mfs))
	times = [ 1, 1, 1, 1, 1]
	for i in xrange(1, time.size):
			temp, times = rules(y[i-2], y[i - 1], times , (vals, mfs))
			#print times	
			y = np.vstack((y, temp))
	#for i in xrange(290,400):
	#	print y[i,6]
		#	if i > 290 and i < 350:
		#		print times
			#print time[i ], y[i]
	#for i in xrange(200,400):
	#	print y[i]
	#for i in xrange(400):
	#	print y[i]
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