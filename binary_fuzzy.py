"""This module tries to implement the boolean version of the earlier model in fuzzy logic"""

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

"""def compute_initial_egfr(egf, hrg, egfr):
	'''Rules---
	If egf is high or hrg is high then egfr is high
	if egf is low and hrg is low then egfr is low'''

	#Membership function for egf
	egf_high = fuzz.gaussmf(egf, 1, 0.25)
	egf_low = fuzz.gaussmf(egf, 0, 0.25)
	'''plt.plot(egf,egf_low)
	plt.show()'''
	#Membership function for hrg and egfr
	'''Technically they should have been calcultated using fuzz.someMf but since they are same to save computation we do this'''
	hrg_high = egf_high  
	hrg_low = egf_low 
	egfr_high = egf_high
	egfr_low = egf_low

	#Antecedents
	a1_egfr = egf_high  # Again here np.fmax should have been used but we save computation
	a2_egfr = egf_low

	#Consequents
	c1_egfr = fuzz.relation_min(a1_egfr, egfr_high)
	c2_egfr = fuzz.relation_min(a2_egfr, egfr_low)
	c_egfr_combined = np.fmax(c1_egfr, c2_egfr)
	return c_egfr_combined

def compute_initial_raf(egfr, akt, raf):
	'''Rules---
	If egfr is high or akt is high then raf is high
	if egfr is low and akt is low then raf is low'''
	#Membership function for egf
	egfr_high = fuzz.gaussmf(egfr, 1, 0.25)
	egfr_low = fuzz.gaussmf(egfr, 0, 0.25)
	
	#Membership function for akt and raf
	'''Technically they should have been calcultated using fuzz.someMf but since they are same to save computation we do this'''
	akt_high = egfr_high  
	akt_low = egfr_low 
	raf_high = egfr_high
	raf_low = egfr_low

	#Antecedents
	a1_raf = egfr_high  # Again here np.fmax should have been used but we save computation
	a2_raf = egfr_low

	#Consequents
	c1_raf = fuzz.relation_min(a1_raf, raf_high)
	c2_raf = fuzz.relation_min(a2_raf, raf_low)
	c_raf_combined = np.fmax(c1_raf, c2_raf)
	return c_raf_combined

def compute_initial_pi3k(egfr, erk, pi3k):
	'''Rules ---
	if egfr is high and erk is low then pi3k is high
	if egfr is low or erk is high then pi3k is low'''

	#Memebership function for egfr
	egfr_high = fuzz.gaussmf(egfr, 1, 0.25)
	egfr_low = fuzz.gaussmf(egfr, 0, 0.25)

	#Membership function for akt and raf
	'''Technically they should have been calcultated using fuzz.someMf but since they are same to save computation we do this'''
	erk_high = egfr_high  
	erk_low = egfr_low 
	pi3k_high = egfr_high
	pi3k_low = egfr_low

	#Antecedents
	a1_pi3k = np.fmin(egfr_high, erk_low)
	a2_pi3k = np.fmax(egfr_low, erk_high)

	#Consequents
	c1_pi3k = fuzz.relation_min(a1_pi3k, pi3k_high)
	c2_pi3k = fuzz.relation_min(a2_pi3k, pi3k_low)
	c_pi3k_combined = np.fmin(c1_pi3k, c2_pi3k)
	return c_pi3k_combined"""

def eval_membership_functions(initial_values):
	'''Technically they should have been calcultated using fuzz.someMf but since they are same to save computation we do this'''
	egf_high = fuzz.gaussmf(initial_values[0], 1, 0.25)
	egf_low = fuzz.gaussmf(initial_values[1], 0, 0.25)
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
	return ((egf_low, egf_high), (hrg_low, hrg_high), (egfr_low, egfr_high), (raf_low, raf_high), (pi3k_low, pi3k_high), (erk_low, erk_high), (akt_low, akt_high))


def compute_egfr(egf_value, hrg_value, initial_values, mfs ):	

	"""Rules---
	If egf is high or hrg is high then egfr is high
	if egf is low and hrg is low then egfr is low"""

	a1_1 = mfs[0][1][initial_values[0] == egf_value] #egf_high[egf == egf_value]
	a1_2 = mfs[1][1][ initial_values[1] == hrg_value ] #hrg_high[hrg == hrg_value]

	if( a1_1.size == 0):
		a1_1 = mfs[0][1][ find_closest(initial_values[0], egf_value)]
	if( a1_2.size == 0):
		a1_2 = mfs[1][1][ find_closest(initial_values[1], hrg_value)]

	a1 = max( a1_1, a1_2 )
	c1 = np.fmin(np.linspace(a1,a1,100), mfs[2][1])

	a2_1 = mfs[0][0][initial_values[0] == egf_value] #egf_low[egf == egf_value]
	a2_2 = mfs[1][0][initial_values[1] == hrg_value] #hrg_low[hrg == hrg_value]

	if( a2_1.size == 0):
		a2_1 = mfs[0][0][ find_closest(initial_values[0], egf_value)]
	if( a2_2.size == 0):
		a2_2 = mfs[1][0][ find_closest(initial_values[1], hrg_value)]

	a2 = min(a2_1, a2_2)
	c2 = np.fmin(np.linspace(a2,a2,100), mfs[2][0])

	c_com = np.fmax(c1, c2)
	
	a =  fuzz.defuzz(initial_values[2], c_com, 'lom')
	return a

def compute_pi3k(egfr_value, erk_value, initial_values, mfs):
	"""Rules ---
	if egfr is high and erk is low then pi3k is high
	if egfr is low or erk is high then pi3k is low"""

	a1_1 = mfs[0][1][initial_values[0] == egfr_value] #egfr_high[egfr == egfr_value]
	a1_2 = mfs[1][0][ initial_values[1] == erk_value] #erk_low[erk == erk_value]

	if( a1_1.size == 0):
		a1_1 = mfs[0][1][ find_closest(initial_values[0], egfr_value)]
	if( a1_2.size == 0):
		a1_2 = mfs[1][0][ find_closest(initial_values[1], erk_value)]

	a1 = min(a1_1 , a1_2)
	c1 = np.fmin( np.linspace(a1, a1, 100), mfs[2][1])

	a2_1 = mfs[0][0][ initial_values[0] == egfr_value] #egfr_low[egfr == egfr_value]
	a2_2 = mfs[1][1][ initial_values[1] == erk_value] #erk_high[erk == erk_value]

	if( a2_1.size == 0):
		a2_1 = mfs[0][0][ find_closest(initial_values[0], egfr_value)]
	if( a2_2.size == 0):
		a2_2 = mfs[1][1][ find_closest(initial_values[1], erk_value)]

	a2 = max(a2_1 , a2_2)
	c2 = np.fmin( np.linspace(a2, a2, 100), mfs[2][0] )

	c_com = np.fmax(c1, c2)
	return fuzz.defuzz(initial_values[2], c_com, 'lom')

def compute_raf(egfr_value, akt_value, initial_values, mfs):
	"""Rules---
	If egfr is high or akt is high then raf is high
	if egfr is low and akt is low then raf is low"""	

	a1_1 = mfs[0][1][initial_values[0] == egfr_value] #egfr_high[egfr == egfr_value]
	a1_2 = mfs[1][1][initial_values[1] == akt_value] #akt_high[akt == akt_value]

	if( a1_1.size == 0):
		a1_1 = mfs[0][1][ find_closest(initial_values[0], egfr_value)]
	if( a1_2.size == 0):
		a1_2 = mfs[1][1][ find_closest(initial_values[1], akt_value)]

	a1 = max( a1_1 , a1_2 )
	#print egfr_value
	c1 = np.fmin( np.linspace(a1, a1, 100), mfs[2][1])

	a2_1 = mfs[0][0][initial_values[0] == egfr_value] #egfr_low[egfr == egfr_value]
	a2_2 = mfs[1][0][initial_values[1] == akt_value] #akt_low[akt == akt_value]

	if( a2_1.size == 0):
		a2_1 = mfs[0][0][ find_closest(initial_values[0], egfr_value)]
	if( a2_2.size == 0):
		a2_2 = mfs[1][0][ find_closest(initial_values[1], akt_value)]

	a2 = min(a2_1 ,a2_2 )
	c2 = np.fmin( np.linspace(a2, a2, 100), mfs[2][0])

	c_com = np.fmax(c1, c2)
	return fuzz.defuzz(initial_values[2], c_com, 'lom')

def rules(initial_cond, (initial_values, mfs)):
	y = np.copy(initial_cond)
	y[2] = compute_egfr(initial_cond[0], initial_cond[1], (initial_values[0], initial_values[1], initial_values[2]), (mfs[0], mfs[1], mfs[2]))
	y[3] = compute_raf(initial_cond[2], initial_cond[6], (initial_values[2], initial_values[6], initial_values[3]), (mfs[2], mfs[6], mfs[3]))
	y[4] = compute_pi3k(initial_cond[2], initial_cond[5], (initial_values[2], initial_values[5], initial_values[4]), (mfs[2], mfs[5], mfs[4]))
	y[5] = initial_cond[3]
	y[6] = initial_cond[4]
	return y

def main():
	#Universal sets
	#c_egfr_aggregated = compute_initial_egfr(egf, hrg, egfr)
	#c_raf_aggregated = compute_initial_raf(egfr, akt, raf)
	#c_pi3k_aggregated = compute_initial_pi3k(egfr, erk, pi3k)

	initial_cond = np.array([1, 1, 0, 0, 0, 0, 0], dtype = "float64")
	time_stop = 10
	y = np.copy(initial_cond)
	step = 1

	egf = np.linspace(0,1,100)
	hrg = egf
	egfr = egf
	akt = egf
	raf = egf
	pi3k = egf
	erk = egf

	vals = (egf, hrg, egfr, raf, pi3k, erk, akt)
	mfs = eval_membership_functions(vals)
	
	while(step < time_stop + 1):
		#global vals
		if(step == 1):
			temp = rules(y, (vals, mfs))
		else:
			temp = rules( y[ step - 1, : ],  (vals, mfs))
		y = np.vstack((y, temp))
		step += 1

	print y
	plt.title("Synch")
	lines = plt.plot(np.arange(time_stop + 1), y[:,5], np.arange(time_stop + 1), y[:,6])	
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