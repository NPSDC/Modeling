import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz
from scipy.interpolate import interp1d
from scipy import interp
import math
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

def centroid(x, mfx):
    """
    Defuzzification using centroid (`center of gravity`) method.
    Parameters
    ----------
    x : 1d array, length M
        Independent variable
    mfx : 1d array, length M
        Fuzzy membership function
    Returns
    -------
    u : 1d array, length M
        Defuzzified result
    See also
    --------
    skfuzzy.defuzzify.defuzz, skfuzzy.defuzzify.dcentroid
    """
    return (x * mfx).sum() / mfx.sum()
                                    


def find_if_middle(c_combined, antecedents, component):
	diction = {}
	for i in c_combined:
		if i in diction.keys():
			diction[i] += 1
		else:
			diction[i] = 1
	
	'''for i in diction.keys():
		if diction[i] == component.size:
			for j in antecedents:
				for k in xrange(len(j)):
					if j[k] == i:
						print j[k],k,i
						return k
		else:
			print -100, i'''
		
	print len(diction.keys())

def ret(val, antecedents):
	for i in xrange(len(antecedents)):
		for j in xrange(len(antecedents[i])):
			if antecedents[i][j] == val :
				return i
							
def compute_egfr_change(egf_value, hrg_value, time_values, initial_values, mfs ):	

	'''Rules---
	If egf is high and time_egf is high or hrg is high and time_hrg  is high then positive_change_egfr is high
	If egf is high and  hrg is low and time_egf  is low then positive_change_egfr is low
	If egf is low and hrg is high and time_hrg is low  then positive_change_egfr is low
	If egf is high and hrg is high and time_hrg is low and time_egf is low then positive_change_egfr is low
	if egf is low  and time_egf is high and hrg is low and time_hrg is high then negative_change_egfr is high
	if egf is low  and and hrg is low and (time_egf is low  or time_hrg is low) then negative_change_egfr is low
	'''
	####Positive CHange#######
	####Antecedent1
	
	a1_1 = mfs[0][1][initial_values[0] == egf_value] #egf_high[egf == egf_value]
	a1_2 = mfs[1][1][ initial_values[1] == hrg_value ] #hrg_high[hrg == hrg_value]
	a1_3 = mfs[3][1][ initial_values[3] == time_values[0]] #time_high[time == time_egf_value]
	a1_4 = mfs[3][1][ initial_values[3] == time_values[1]] #time_high[time == time_hrg_value ]

	if( a1_1.size == 0):
		f = interp1d(initial_values[0], mfs[0][1])
		a1_1 = f(egf_value)
	if( a1_2.size == 0):
		a1_2 = interp1d( initial_values[1], mfs[1][1])(egf_value)

	a1 = max( min(a1_1, a1_3) , min(a1_2, a1_4))
	#Consequent 1
	c1 = np.fmin(a1, mfs[2][3])  #mfs[2][3] is positve_change_egfr_high

	####Antecedent2
	a2_1 = a1_1 #egf_high[egf == egf_value]
	a2_2 = mfs[1][0][initial_values[1] == hrg_value] #hrg_low[hrg == hrg_value]
	a2_3 = mfs[3][0][ initial_values[3] == time_values[0]] #time_low[time == time_egf_value]
	
	
	if( a2_2.size == 0):
		a2_2 = interp(hrg_value, initial_values[1],  mfs[1][0])

	a2 = min(a2_1, a2_2, a2_3 )
	#Consequent 2
	c2 = np.fmin(a2, mfs[2][2])	#mfs[2][2] is positve_change_egfr_low

	####Antecedent3
	a3_1 = mfs[0][0][initial_values[0] == egf_value] #egf_low[egf == egf_value]
	a3_2 = a1_2 #hrg_high[hrg == hrg_value]
	a3_3 = mfs[3][0][ initial_values[3] == time_values[1]] #time_low[time == time_hrg_value]

	if( a3_1.size == 0):
		a3_1 = interp( egf_value, initial_values[0], mfs[0][0])
	
	
	a3 = min(a3_1, a3_2, a3_3 )

	#Consequent 3
	c3 = np.fmin(a3, mfs[2][2])	#mfs[2][2] is positve_change_egfr_low

	#Antecedent 4
	a4_1 = a1_1 #egf_high[egf == egf_value]
	a4_2 = a1_2 #hrg_high[hrg == hrg_value]
	a4_3 = a2_3 #time_low[time == time_egf_value]
	a4_4 = a3_3 #time_low[time == time_hrg_value]
	a4 = min(a4_1, a4_2, max(a4_3, a4_4))

	#Consequent 4
	c4 = np.fmin(a4, mfs[2][2])  #mfs[2][2] is positve_change_egfr_low
	c_com_positive = np.fmax(np.fmax(c1, c4), np.fmax(c2, c3))

	try:
		pos_change = fuzz.defuzz(initial_values[2][1], c_com_positive, 'centroid' ) #initial_values[2][1] is positive_change_egfr
	except AssertionError as e:
		pos_change = 0

	

	#plt.plot(mfs[2][1], c_com_positive)
	#plt.show()
	#print a1,a2,a3,a4s
	#print c_com_positive
	#### Negative Change ######

	####Antecedent 6
	a5_1 = a3_1 #egf_low[egf == egf_value]
	a5_2 = a2_2 #hrg_low[hrg == hrg_value]
	a5_3 = a1_3 #time_high[time == time_egf_value]
	a5_4 = a1_4 #time_high[time == time_hrg_value]
		
	a5 = min(a5_1, a5_2, a5_3, a5_4)
	#Consequent 5
	c5 = np.fmin(a5, mfs[2][5]) #mfs[2][5] is negative_change_egfr_high

	#Antecedent 6
	a6_1 = a5_1 #egf_low[egf == egf_value]
	a6_2 = a5_2 #hrg_low[hrg == hrg_value]
	a6_3 = a4_3 #time_low[time == time_egf_value]
	a6_4 = a4_4 #time_low[time == time_hrg_value]
	a6 = min(a6_1, a6_2, max(a6_3, a6_4))

	#Consequent 6
	c6 = np.fmin(a6, mfs[2][4]) #mfs[2][4] is negative_change_egfr_low

	c_com_negative = np.fmax(c5, c6)
	try:
		neg_change = fuzz.defuzz(initial_values[2][2], c_com_negative, 'centroid') #initial_values[2][2] is negative_change_egfr
	except AssertionError as e:
		neg_change = 0
	#print a1,a2,a3,a4,a5,a6,pos_change,neg_change
	#print neg_change,pos_change, neg_change + pos_change
	return pos_change + neg_change

def compute_raf_change(egfr_value, akt_value, time_values, initial_values, mfs):
	"""Rules---
	1)If egfr is high and time_egfr is high or akt is high and time_akt is high then positive_change_raf is high and negative_change_raf is low
	2)If egfr is high and akt is high and time_egfr is low  and time_akt is low then positive_change_raf is low and negative_change_raf is low
	3)If egfr is high and akt is low and time_egfr is low  then positve_change_raf is low and negative_change_raf is low
	4)If egfr is low and akt is high and time_akt is low then positive_change_raf is low and negative_change_raf is low
	5)If egfr is low and akt is low  and time_egfr is high and time_akt is high then negative_change_raf is high and positive_change_raf is low
	6)If egfr is low and akt is low  and (time_egfr is low or time_akt is low) then negative_change_raf is low and positive_change_raf is low
	7)If egfr is mid and akt is mid then positive_change_raf is low and negative_change_raf is low
	8)If egfr is mid and akt is low then positive_change_raf is low and negative_change_raf is low
	9)If egfr is low and akt is mid then positive_change_raf is low and negative_change_raf is low
	10)If egfr is high and akt is mid and time_egfr is low  then positve_change_raf is low and negative_change_raf is low
	11)If egfr is mid and akt is high and time_akt is low then positive_change_raf is low and negative_change_raf is low
	"""
	###Positive Change
	#Antecedent 1
	a1_1 = mfs[0][1][initial_values[0][0] == egfr_value] #egfr_high[egfr == egfr_value]
	a1_2 = mfs[1][1][initial_values[1][0] == akt_value]  #akt_high[akt == akt_value]
	a1_3 = mfs[3][1][initial_values[3] == time_values[0]] #time_high[time == time_egfr_value]
	a1_4 = mfs[3][1][initial_values[3] == time_values[1]] #time_high[time == time_akt_value]

	if( a1_1.size == 0):
		a1_1 = interp( egfr_value, initial_values[0][0], mfs[0][1])
	if( a1_2.size == 0):
		a1_2 = interp(akt_value, initial_values[1][0],  mfs[1][1])
	
	a1 = max( min(a1_1 , a1_3), min(a1_2, a1_4) )
	
	
	#Consequent 1
	c1_pos = np.fmin( a1, mfs[2][3])  #mfs[2][3] is positive_change_raf_high
	c1_neg = np.fmin(a1, mfs[2][4])
	
	##Antecedent 2
	a2_1 = a1_1 #egfr_high[egfr == egfr_value]
	a2_2 = a1_2 #akt_high[akt == akt_value]
	a2_3 = mfs[3][0][initial_values[3] == time_values[0]] #time_low[time == time_egfr_value]
	a2_4 = mfs[3][0][initial_values[3] == time_values[1]] #time_low[time == time_akt_value]
	a2 = min(a2_1, a2_2, a2_3, a2_4)

	#Consequent 2
	c2_neg = np.fmin(a2, mfs[2][2])  #mfs[2][2] is positive_change_raf_low
	c2_pos = np.fmin(a2, mfs[2][4])
	c_com_negative = np.fmax(c1_neg, c2_neg)
	c_com_positive = np.fmax(c1_pos, c2_pos)
	#Antecedent 3
	a3_1 = a1_1 #egfr_high[ egfr == egfr_value]
	a3_2 = mfs[1][0][ initial_values[1][0] == akt_value] #akt_low[akt == akt_value]
	a3_3 = a2_3 #time_low[ time == time_egfr_value]

	if(a3_2.size == 0):
		a3_2 = interp(akt_value, initial_values[1][0], mfs[1][0])

	a3 = min(a3_1, a3_2, a3_3)

	#Consequent 3
	c3_pos = np.fmin(a3, mfs[2][2]) #mfs[2][2] is positive_change_raf_low
	c3_neg = np.fmin(a3, mfs[2][4])
	c_com_negative = np.fmax(c_com_negative, c3_neg)
	c_com_positive = np.fmax(c_com_positive, c3_pos)

	#Antecedent 4
	a4_1 = mfs[0][0][ initial_values[0][0] == egfr_value ] #egfr_low[egfr == egfr_value]
	a4_2 = a1_2 #akt_high[akt == akt_value]
	a4_3 = a2_4 #time_low[time == time_akt_value]

	if(a4_1.size == 0):
		a4_1 = interp(egfr_value, initial_values[0][0], mfs[0][0] )

	a4 = min(a4_1, a4_2, a4_3)	

	#Consequent 4
	c4_pos = np.fmin(a4, mfs[2][2]) #mfs[2][2] is positive_change_raf_low
	c4_neg = np.fmin(a4, mfs[2][4])
	c_com_negative = np.fmax(c_com_negative, c4_neg)
	c_com_positive = np.fmax(c_com_positive, c4_pos)
	#A7
	f = interp1d(initial_values[0][0], mfs[0][6])
	a7_1 = f(egfr_value)
	f = interp1d(initial_values[1][0], mfs[1][6])
	a7_2 = f(akt_value)
	a7 = min(a7_1, a7_2)
	c7_pos = np.fmin(a7, mfs[2][2])
	c7_neg = np.fmin(a7, mfs[2][4])
	c_com_negative = np.fmax(c_com_negative, c7_neg)
	c_com_positive = np.fmax(c_com_positive, c7_pos)
	#A8
	a8_1 = a7_1
	a8_2 = a3_2
	a8 = min(a8_1, a8_2)
	c8_pos = np.fmin(a8, mfs[2][2])
	c8_neg = np.fmin(a8, mfs[2][4])
	c_com_negative = np.fmax(c_com_negative, c8_neg)
	c_com_positive = np.fmax(c_com_positive, c8_pos)
	#A9
	a9_1 = a4_1
	a9_2 = a7_2
	a9 = min(a9_1, a9_2)
	c9_pos = np.fmin(a9, mfs[2][2])
	c9_neg = np.fmin(a9, mfs[2][4])
	c_com_negative = np.fmax(c_com_negative, c9_neg)
	c_com_positive = np.fmax(c_com_positive, c9_pos)
	#A10
	a10_1 = a1_1
	a10_2 = a7_2
	a10_3 = a2_3
	a10 = min(a10_1, a10_2, a10_3)
	c10_pos = np.fmin(a10, mfs[2][2])
	c10_neg = np.fmin(a10, mfs[2][4])
	c_com_negative = np.fmax(c_com_negative, c10_neg)
	c_com_positive = np.fmax(c_com_positive, c10_pos)
	#A11
	a11_1 = a7_1
	a11_2 = a1_2
	a11_3 = a2_4
	a11= min(a11_1, a11_2, a11_3)
	c11_pos = np.fmin(a11, mfs[2][2])
	c11_neg = np.fmin(a11, mfs[2][4])
	c_com_negative = np.fmax(c_com_negative, c11_neg)
	c_com_positive = np.fmax(c_com_positive, c11_pos)
	'''c_pos = np.fmax(np.fmax(np.fmax(c7, c8), np.fmax(c9,c10)), c11)
	c_com_positive = np.fmax(np.fmax(c1, c2), np.fmax(c3, c4))
	c_com_positive = np.fmax(c_com_positive, c_pos)'''
	#find_if_middle(c_com_positive, ((a1_1, a1_2, a1_3,a1_4), (a2_1, a2_2, a2_3, a2_4), (a3_1, a3_2, a3_3), (a4_1, a4_2, a4_3)), initial_values[2][1])
	'''if(math.log10(a1) >=-6 and math.log10(a1) <= -5):
		plt.plot(initial_values[2][1], mfs[2][1])
		plt.show()
		plt.plot(initial_values[2][1], c1)
		plt.show()
		plt.plot(initial_values[2][1], c3)
		plt.show()
		plt.plot(initial_values[2][1], c_com_positive)
		plt.show()
		print c_com_positive'''
	'''plt.plot(initial_values[2][1], c2)
	plt.show()'''
	#if(initial_values[3][5] == time_values[1]):
	'''	plt.plot(initial_values[2][1], c1)
		plt.show()
		plt.plot(initial_values[2][1], c2)
		plt.show()
		plt.plot(initial_values[2][1], c3)
		plt.show()
		plt.plot(initial_values[2][1], c4)
		plt.show()'''
	#plt.plot(initial_values[2][1], c_com_positive)
	#plt.show()
	#	print c_com_positive
	#print a1_1, a1_2, a4_1, a3_2
	#print a1,a2,a3,a4
	'''plt.plot(initial_values[2][0], mfs[0][0])
	plt.show()'''
	#plt.plot(initial_values[2][1], c3)
	#plt.show()
	'''plt.plot(initial_values[2][1], c4)
	print a4
	plt.show()'''
	
	'''plt.plot(initial_values[2][1], c_com_positive)
	plt.show()'''
	

 	###Negative
	#Antecedent 5
	a5_1 = a4_1 #egfr_low[egfr == egfr_value]
	a5_2 = a3_2 #akt_low[akt == akt_value]
	a5_3 = a1_3 #time_high[time == time_egfr_value]
	a5_4 = a1_4 #time_high[time == time_akt_value]

	a5 = min(a5_1 ,a5_2, a5_3, a5_4 )

	#Consequent 5
	c5_neg = np.fmin( a5, mfs[2][5]) #mfs[2][5] is raf_negative_change_high
	c5_pos = np.fmin(a5, mfs[2][2])
	c_com_negative = np.fmax(c_com_negative, c5_neg)
	c_com_positive = np.fmax(c_com_positive, c5_pos)
	#Antecedent 6
	a6_1 = a4_1  #egfr_low[egfr == egfr_value]
	a6_2 = a3_2  #akt_low[akt == akt_value]
	a6_3 = a2_3  #time_low[time == time_egfr_value]
	a6_4 = a2_4  #time_low[time == time_akt_value]

	a6 = min(a6_1, a6_2,max( a6_3, a6_4))
#	print a5_1,a5_2,a5_3,a5_4,a6_3,a6_4	
	#Consequent 6
	c6_neg = np.fmin(a6, mfs[2][4]) #mfs[2][4] is raf_negative_change_low
	c6_pos = np.fmin(a6, mfs[2][2])
	#print a5,a6
	c_com_negative = np.fmax(c_com_negative, c6_neg)
	c_com_positive = np.fmax(c_com_positive, c6_pos)
	try:
		'''	plt.plot(initial_values[2][2], c5)
		plt.plot(initial_values[2][2], c6)
		plt.show()
		plt.plot(initial_values[2][2], c_com_negative)
		plt.show()'''
		neg_change = centroid(initial_values[2][2], c_com_negative) #initial_values[2][2] is raf_negative_change
	except AssertionError as e:
		neg_change = 0
	try :
		pos_change = centroid(initial_values[2][1], c_com_positive) #initial_values[2][1] is positve_change_raf
	except AssertionError as e:
		pos_change = 0
	#print a6_1,a6_2,a6_3,a6_4
	#print a7, a8, a9, a10, a11
	#print a1,a2,a3,a4,a5,a6
	#find_if_middle(c_com_negative, ((a5_1, a5_2, a5_3,a5_4), (a6_1, a6_2, a6_3, a6_4)), initial_values[2][1])
	#if(initial_values[3][39] == time_values[0]):
	print pos_change,neg_change, pos_change + neg_change
#plt.plot(initial_values[2][2], c_com_negative)	print neg_change,pos_change, neg_change + pos_change
	return pos_change + neg_change

def compute_pi3k_change(egfr_value, erk_value, time_values, initial_values, mfs):
	"""Rules ---
	if egfr is high and time_egfr is high and erk is low and time_erk is high then postive_change_pi3k is high
	if egfr is high and erk is low and (time_egfr is low or time_erk is low) then postive_change_pi3k is low
	if egfr is low and time_egfr is high then  negative_change_pi3k is high
	if erk is high and time_erk is high then negative_change_pi3k is high
	if egfr is low and erk is high and time_egfr is low and time_erk is low then negative_change_pi3k is low
	if egfr is low and erk is low and time_egfr is low then negative_change_pi3k is low
	if egfr is high and erk is high and time_erk is low then negative_change_pi3k is low """

	######Positive Change
	#Antecedent 1
	f = interp1d(initial_values[0][0], mfs[0][1])  
	a1_1 = f(egfr_value) #egfr_high[egfr == egfr_value]
	f = interp1d(initial_values[1][0], mfs[1][0])
	a1_2 = f (erk_value) #erk_low[erk == erk_value
	f = interp1d(initial_values[3], mfs[3][1])
	a1_3 = f(time_values[0]) #time_high[time == time_egfr_value]
	a1_4 = f(time_values[1]) #time_high[time == time_erk_value ]

	a1 = min(a1_1 , a1_2, a1_3, a1_4)
	c1 = np.fmin( a1, mfs[2][3]) #positive_change_pi3k is high

	#Antecedent 2
	a2_1 = a1_1 #egfr_high[egfr == egfr_value]
	a2_2 = a1_2 #erk_low[erk == erk_value]
	f = interp1d(initial_values[3], mfs[3][0]) #time_low
	a2_3 = f(time_values[0]) #time_low[time == time_egfr_value]
	a2_4 = f(time_values[1]) #time_low[time == time_erk_value]

	a2 = min(a2_1 , a2_2, max(a2_3, a2_4))
	c2 = np.fmin( a2, mfs[2][2] ) #positive_change_pi3k is low

	c_com_positive = np.fmax(c1, c2)
	pos_change = fuzz.defuzz(initial_values[2][1], c_com_positive, 'centroid')
	#print a1_1,a1_2,a1_3,a1_4,a2_3,a2_4
	#print egfr_value
	#####Negative Change
	##Antecedent3
	f = interp1d(initial_values[0][0], mfs[0][0]) #
	a3_1 = f(egfr_value) #egfr_low[egfr == egfr_value]
	a3_2 = a1_3 #time_high[time == time_egfr_value]
	a3 = min(a3_1, a3_2)
	c3 = np.fmin(a3, mfs[2][5]) #negative_change_pi3k is high

	#Antecedent 4
	f = interp1d(initial_values[1][0], mfs[1][1])
	a4_1 = f(erk_value) #erk_high[erk == erk_value]
	a4_2 = a1_4 #time_high[time == time_erk_value]
	a4 = min(a4_1, a4_2)
	c4 = np.fmin(a4, mfs[2][5])  #negative_change_pi3k is high

	#Antecedent 5
	a5_1 = a3_1 #egfr_low[egfr == egfr_value]
	a5_2 = a4_1 #erk_high[erk == erk_value]
	a5_3 = a2_3 #time_low[time == time_egfr_value]
	a5_4 = a2_4 #time_low[time == time_erk_value]
	a5 = min(a5_1, a5_2, a5_3, a5_4)
	c5 = np.fmin(a5, mfs[2][4])  #negative_change_pi3k is low

	#Antecedent 6
	a6_1 = a3_1 
	a6_2 = a1_2
	a6_3 = a2_3
	a6 = min(a6_1, a6_2, a6_3)
	c6 = np.fmin(a6, mfs[2][4])

	#Antecedent 7
	a7_1 = a1_1
	a7_2 = a4_1
	a7_3 = a2_4
	a7 = min(a7_1, a7_2, a7_3)
	c7 = np.fmin(a7, mfs[2][4])
	c_com_negative = np.fmax(c3, np.fmax(np.fmax(c4, c5), np.fmax(c6, c7)))
	
	neg_change = fuzz.defuzz(initial_values[2][2], c_com_negative, 'centroid')
	#print pos_change,neg_change,time_values
	return pos_change + neg_change

def compute_erk_change(raf_value, time_value, initial_values, mfs):
	"""Rules-
		If raf is high and time is high then positive_change_erk is high
		If raf is high and time is low then positive_change_erk is low
		If raf is low and time is high then negative_change_erk is high
		If raf is low and time is low then negative_change_erk is low"""

	###Positive Change

	#Antecedent 1
	f = interp1d(initial_values[0][0], mfs[0][1])
	a1_1 = f(raf_value) #raf_high[raf == raf_value]
	f = interp1d(initial_values[2], mfs[2][1])
	a1_2 = f(time_value) #time_high[time == time_value]

	a1 = min(a1_1, a1_2)
	c1 = np.fmin( a1, mfs[1][3]) #mfs[1][3] is positive_change_erk_high

	#Antecedent 2
	a2_1 = a1_1
	f = interp1d(initial_values[2], mfs[2][0]) #time_low[time == time_value]
	a2_2 = f(time_value)

	a2 = min(a2_1, a2_2)
	c2 = np.fmin( a2, mfs[1][2]) #mfs[1][2] is positive_change_raf_low

	c_com_positive = np.fmax(c1,c2)
	pos_change = fuzz.defuzz( initial_values[1][1], c_com_positive, 'centroid') #initial_values[1][1] is positive_change_erk

	###Negative Change

	#Antecedent 3
	f = interp1d(initial_values[0][0], mfs[0][0])
	a3_1 = f(raf_value) #raf_low[raf == raf_value]
	a3_2 = a1_2 #time_high[time == time_value]

	a3 = min(a3_1,a3_2)
	c3 = np.fmin(a3, mfs[1][5]) #mfs[1][3] is negative_change_erk_high

	#Antecedent 4
	a4_1 = a3_1 #raf_low[raf == raf_value]
	a4_2 = a2_2 #time_low[time == time_value]

	a4 = min(a4_1, a4_2)
	c4 = np.fmin(a4, mfs[1][4]) #mfs[1][4] is negative_change_erk_low

	c_com_negative = np.fmax(c3, c4)
	neg_change = fuzz.defuzz(initial_values[1][2], c_com_negative, 'centroid') #initial_values[1][2] is negative_change_erk
	
	#print pos_change, neg_change
	return pos_change + neg_change

def compute_akt_change(pi3k_value, time_value, initial_values, mfs):
	"""Rules-
		If pi3k is high and time is high then positive_change_akt is high
		If pi3k is high and time is low then positive_change_akt is low
		If pi3k is low and time is high then negative_change_akt is high
		If pi3k is low and time is low then negative_change_akt is low"""

	###Positive Change
	#Antecedent 1
	f = interp1d(initial_values[0][0], mfs[0][1])	
	a1_1 = f(pi3k_value) #pi3k_high[pi3k == pi3k_value]
	f = interp1d(initial_values[2], mfs[2][1])
	a1_2 = f(time_value) #time_high[time == time_value]

	a1 = min(a1_1, a1_2)
	c1 = np.fmin(a1, mfs[1][3]) #positive_change_akt is high

	#Antecedent 2
	a2_1 = a1_1 #pi3k_high[pi3k == pi3k_value]
	f = interp1d(initial_values[2], mfs[2][0])
	a2_2 = f(time_value) #time_low[time == time_value]

	a2 = min(a2_1, a2_2)
	c2 = np.fmin( a2, mfs[1][2]) #positive_change_akt is low

	c_com_positive = np.fmax(c1,c2)
	pos_change = fuzz.defuzz( initial_values[1][1], c_com_positive, 'centroid') #initial_values[1][1] is positive_change_akt

	###Negative Change

	#Antecedent 3
	f = interp1d(initial_values[0][0], mfs[0][0])
	a3_1 = f(pi3k_value) #pi3k_low[pi3k == pi3k_value]
	a3_2 = a1_2 #time_high[time == time_value]

	a3 = min(a3_1, a3_2)
	c3 = np.fmin(a3, mfs[1][5]) #mfs[1][5] is negative_change_akt_high

	#Antecedent 4
	a4_1 = a3_1 #pi3k_low[pi3k == pi3k_value]
	a4_2 = a2_2 #time_low[time == time_value]

	a4 = min(a4_1, a4_2)
	c4 = np.fmin(a4, mfs[1][4]) #mfs[1][4] is negative_change_akt_low

	c_com_negative = np.fmax(c3, c4)
	neg_change = fuzz.defuzz(initial_values[1][2], c_com_negative, 'centroid') #initial_values[1][2] is negative_change_akt

	return pos_change + neg_change

def compute_egfr(change_egfr_reflected, not_updated, initial_cond, time_egfr_index, initial_values, mfs):

	egfr = initial_cond[2]
	time_size = mfs[7][0].size
	if(time_egfr_index[0] > time_size or time_egfr_index[1] > time_size):
		change_egfr_reflected = egfr
	
	if 0 in not_updated :		
		time_egfr_index[0] = time_egfr_index[0] + 1		
	else:
		time_egfr_index[0] = 2

	if 1 in not_updated:
		time_egfr_index[1] = time_egfr_index[1] + 1		
	else:
		time_egfr_index[1] = 2

	time_egf_index = time_egfr_index[0] - 1
	time_hrg_index = time_egfr_index[1] - 1

	if(time_egf_index >= initial_values[7].size):
		time_egf_index = initial_values[7].size - 1

	if(time_hrg_index >= initial_values[7].size):	
		time_hrg_index = initial_values[7].size - 1

	temp =  change_egfr_reflected + compute_egfr_change(initial_cond[0], initial_cond[1],\
    		(initial_values[7][time_egf_index], initial_values[7][time_hrg_index]), \
			(initial_values[0], initial_values[1], initial_values[2], initial_values[7]), (mfs[0], mfs[1], mfs[2], mfs[7]))

	if(temp >= 0 and temp <= 1 ):
		egfr = temp
	
	return (egfr, time_egfr_index, change_egfr_reflected)

def compute_raf(change_raf_reflected, not_updated, initial_cond, time_raf_index, initial_values, mfs):

	
	raf = initial_cond[3]
	time_size = mfs[7][0].size
	if((time_raf_index[0] > time_size and time_raf_index[0]%time_size == 1) or (time_raf_index[1] > time_size and time_raf_index[1] % time_size ==1 )):
		change_raf_reflected = raf
		
	if 2 in not_updated :
		time_raf_index[0] = time_raf_index[0] + 1
	else:
		time_raf_index[0] = 2

	if  6 in not_updated:
		time_raf_index[1] = time_raf_index[1] + 1
	else:
		time_raf_index[1] = 2

	time_egfr_index = time_raf_index[0] - 1
	time_akt_index = time_raf_index[1] - 1

	if(time_egfr_index >= initial_values[7].size):
		time_egfr_index = initial_values[7].size - 1

	if(time_akt_index >= initial_values[7].size):	
		time_akt_index = initial_values[7].size - 1
	
	temp = change_raf_reflected + compute_raf_change(initial_cond[2], initial_cond[6],\
		   (initial_values[7][ time_egfr_index ], initial_values[7][ time_akt_index ] ),\
	 	   (initial_values[2], initial_values[6], initial_values[3], initial_values[7]), (mfs[2], mfs[6], mfs[3], mfs[7]))

	if(temp >= 0 and temp <= 1):
		raf = temp

	return (raf, time_raf_index, change_raf_reflected)

def compute_pi3k(change_pi3k_reflected, not_updated, initial_cond, time_pi3k_index, initial_values, mfs ):

	pi3k = initial_cond[4]
	time_size = mfs[7][0].size

	if((time_pi3k_index[0] > time_size and time_pi3k_index[0] % time_size == 1) or \
		(time_pi3k_index[1] > time_size and time_pi3k_index[1] % time_size == 1)):
		#print "yes",time_pi3k_index
		change_pi3k_reflected = pi3k

	if 2 in not_updated :		
		time_pi3k_index[0] += 1		
	else:
		time_pi3k_index[0] =  2

	if 5 in not_updated:
		time_pi3k_index[1] += 1
	else:
		time_pi3k_index[1] = 2

	time_egfr_index = time_pi3k_index[0] - 1
	time_erk_index = time_pi3k_index[1] - 1

	if(time_egfr_index >= time_size):
		time_egfr_index = time_size - 1

	if(time_erk_index >= time_size):
		time_erk_index = time_size - 1

	temp = change_pi3k_reflected + compute_pi3k_change(initial_cond[2], initial_cond[5], \
		   (initial_values[7][ time_egfr_index], initial_values[7][time_erk_index]), \
		   (initial_values[2], initial_values[5], initial_values[4], initial_values[7]), (mfs[2], mfs[5], mfs[4], mfs[7]))

	if(temp >= 0 and temp <= 1):
		pi3k = temp

	return (pi3k, time_pi3k_index, change_pi3k_reflected)

def compute_erk(change_erk_reflected, not_updated, initial_cond, time_erk_index, initial_values, mfs ):

	erk = initial_cond[5]
	time_size = mfs[7][0].size

	if(time_erk_index > time_size and time_erk_index % time_size == 1):
		change_erk_reflected = erk

	if 3 in not_updated:		
		time_erk_index = time_erk_index + 1
	else:
		time_erk_index = 2

	time_raf_index = time_erk_index - 1

	if(time_raf_index >= time_size):
		time_raf_index = time_size - 1

	temp = change_erk_reflected + compute_erk_change(initial_cond[3], initial_values[7][ time_raf_index], \
		 (initial_values[3], initial_values[5], initial_values[7]), (mfs[3], mfs[5], mfs[7]))	

	if (temp >= 0 and temp <= 1):
		erk = temp

	return (erk, time_erk_index, change_erk_reflected)

def compute_akt(change_akt_reflected, not_updated, initial_cond, time_akt_index, initial_values, mfs ):

	akt = initial_cond[6]
	time_size = mfs[7][0].size

	if(time_akt_index > time_size and time_akt_index % time_size == 1):
		change_akt_reflected = akt

	if 4 in not_updated:
		time_akt_index += 1		
	else:
		time_akt_index = 2	

	time_pi3k_index = time_akt_index - 1

	if(time_pi3k_index >= time_size):
		time_pi3k_index = time_size - 1

	temp = change_akt_reflected + compute_akt_change(initial_cond[4], initial_values[7][ time_pi3k_index], \
		 (initial_values[4], initial_values[6], initial_values[7]), (mfs[4], mfs[6], mfs[7]))
		
	if(temp >= 0 and temp <= 1):
		akt = temp

	return (akt, time_akt_index, change_akt_reflected)


def rules(change_reflected, prev_cond, initial_cond, time_indexes, (initial_values, mfs)):	
	y = np.copy(initial_cond)
	not_updated = []  #contains the indexes not updated
	for i in xrange(len(prev_cond)):
		if str(prev_cond[i]) == str(initial_cond[i]):
			not_updated.append(i)	
		
	y[2], time_indexes[0], change_reflected[2] = compute_egfr(change_reflected[2], not_updated, initial_cond, time_indexes[0], initial_values, mfs )
	y[3], time_indexes[1], change_reflected[3] = compute_raf(change_reflected[3], not_updated, initial_cond, time_indexes[1], initial_values, mfs )
	#y[4], time_indexes[2], change_reflected[4] = compute_pi3k(change_reflected[4], not_updated,  initial_cond, time_indexes[2], initial_values, mfs )
	#y[5], time_indexes[3], change_reflected[5] = compute_erk(change_reflected[5], not_updated, initial_cond, time_indexes[3], initial_values, mfs )
	#y[6], time_indexes[4], change_reflected[6] = compute_akt(change_reflected[6], not_updated, initial_cond, time_indexes[4], initial_values, mfs )
	return (y,time_indexes, change_reflected)

def main():
	
	initial_cond = np.array([1, 1, 0, 0, 0, 0, 0], dtype = "float64")
	time_stop = 10
	y = np.copy(initial_cond)
	change_reflected = np.copy(initial_cond) # This is the array to which we would adding the changing values
	y.resize(1, 7)
	step = 1

	egf = np.linspace(0, 1, 10)
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

	negative_change_egfr = np.linspace(-1, 0, 10)
	negative_change_raf = negative_change_egfr
	negative_change_pi3k = negative_change_egfr
	negative_change_erk = negative_change_egfr
	negative_change_akt = negative_change_egfr

	time = np.linspace(0, 10, 101)
	time_1 = np.linspace(0, 1, 10)
	vals = (egf, hrg, (egfr, positive_change_egfr, negative_change_egfr),\
		   (raf, positive_change_raf, negative_change_raf),\
		   (pi3k, positive_change_pi3k, negative_change_pi3k),\
		   (erk, positive_change_erk, negative_change_erk),\
		   (akt, positive_change_akt, negative_change_akt), time_1)

	mfs = eval_membership_functions(vals)

	times = [ [1, 1], [1, 1], [1, 1], 1, 1] #provides the initial time inputs (indexes 0 for egfr, 1 for raf, 2 for pi3k and so on)
	for i in xrange(1, time.size):	#Also the for egfr 0 is for egf and 1 for hrg; for raf 0 for egfr and 1 for akt; for pi3k 0 for egfr and 1 for erk
		temp, times, change_reflected = rules(change_reflected, y[i-2], y[i - 1], times , (vals, mfs))
		y = np.vstack((y, temp))		
		#print y[i][4], times[1], change_reflected[4]
		#if i < 31:
			#print y[i][2], times[1]
			#print y[i][6], times[4], change_reflected[6] 
	
	
	plt.title("Synch")
	
	lines = plt.plot(time, y[:,3])	
	
	plt.legend(loc='upper right')
	plt.xlabel('Time')
	plt.ylabel('Species')
	plt.axis([-0.2,10.1,-0.05,1.2])
	plt.grid(True)
	plt.show()	

if __name__ == "__main__":
	main()