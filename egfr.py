import numpy as np
import skfuzzy as fuzz
import matplotlib.pyplot as plt

def transform_gaussmf(x, param, val):
	point = param[0]
	y = fuzz.gaussmf(x, param[0], param[1])
	idx = np.searchsorted(x, point, 'left')
	if(idx < x.size - 1):
		if(x[idx] == x[idx + 1]):
			idx += 1

	if(val == 0):
		for i in xrange(idx):
			y[i] = 0
	else:
		for i in xrange(idx + 1, x.size):
			y[i] = 0

	return y

def calculate_egfr_mfs(time, egf, egfr, egf_index, t):
	size = egf.size
	#egfr = np.linspace(0, egf[egf_index], size)
	time_low = fuzz.zmf(time, 0, time[egf.size*t - t] - time[egf_index*t] + time[t] )
	time_high = fuzz.smf(time, 0, time[egf.size*t - t] - time[egf_index*t] + time[t] )
	egfr_low = fuzz.gaussmf(egfr, 0, egf[egf_index]/20)
	egfr_high = transform_gaussmf(egfr, (egf[egf_index], egf[egf_index]/20), 1)
	return ((time_low, time_high, egfr_low, egfr_high))

def calculate_raf_mfs( egfr, raf):
	egfr_low = fuzz.zmf(egfr, 0, 0.95)
	egfr_high = fuzz.smf(egfr, 0.1, 0.9)
	raf_low = fuzz.gaussmf(raf, 0, egfr[egfr.size - 1]/20)
	raf_high = fuzz.gaussmf(raf, egfr[egfr.size - 1], egfr[egfr.size - 1]/20)
	return (egfr_low, egfr_high, raf_low, raf_high)

def calculate_erk_mfs(raf, erk):
	raf_low = fuzz.zmf(raf, 0, 0.95)
	raf_high = fuzz.smf(erk, 0.1, 0.9)
	erk_low = fuzz.gaussmf(erk, 0, raf[raf.size - 1]/20)
	erk_high = fuzz.gaussmf(erk, raf[raf.size - 1], raf[raf.size - 1]/20)
	return (raf_low, raf_high, erk_low, erk_high)

def calculate_pi3k_mfs(egfr, erk, pi3k):
	egfr_low = fuzz.zmf(egfr, 0, 0.95)
	egfr_high = fuzz.smf(egfr, 0.1, 0.9)
	erk_low = fuzz.zmf(erk, 0, 0.9)
	erk_high = fuzz.smf(erk, 0, 0.9)
	pi3k_low = fuzz.gaussmf(pi3k, 0, egfr[egfr.size - 1]/20)
	pi3k_high = fuzz.gaussmf(pi3k, egfr[egfr.size - 1], egfr[egfr.size - 1]/20)
	return ((egfr_low, egfr_high), (erk_low, erk_high), (pi3k_low, pi3k_high))

def calculate_akt_mfs(pi3k, akt):
	pi3k_low = fuzz.zmf(pi3k, 0, 0.95)
	pi3k_high = fuzz.smf(pi3k, 0, 0.9)
	akt_low = fuzz.gaussmf(akt, 0, akt[akt.size -1]/20)
	akt_high = fuzz.gaussmf(akt, akt[akt.size -1], akt[akt.size -1]/20)
	return ((pi3k_low, pi3k_high), (akt_low, akt_high))

def calculate_egfr(time_mfs, egfr_mfs, egfr, time_index):
	a1 = time_mfs[0][time_index]
	c1 = np.fmin(a1, egfr_mfs[0])
	a2 = time_mfs[1][time_index]
	c2 = np.fmin(a2, egfr_mfs[1])
	c_com = np.fmax(c1, c2)
	try:
		egfr_val = fuzz.defuzz(egfr, c_com, 'centroid')
	except AssertionError as e:
		egfr_val = 0
	return egfr_val

def calculate_raf(egfr_mfs, raf_mfs, raf, egfr_index):
	a1 =  egfr_mfs[0][egfr_index]
	c1 = np.fmin(a1, raf_mfs[0])
	a2 = egfr_mfs[1][egfr_index]
	c2 = np.fmin(a2, raf_mfs[1])
	c_com = np.fmax(c1, c2)
	try:
		raf_val = fuzz.defuzz(raf, c_com, 'centroid')
	except AssertionError as e:
		raf_val = 0
	return raf_val

def calculate_erk(raf_mfs, erk_mfs, erk, raf_index):
	a1 = raf_mfs[0][raf_index]
	c1 = np.fmin(a1, erk_mfs[0])
	a2 = raf_mfs[1][raf_index]
	c2 = np.fmin(a2, erk_mfs[1])
	c_com = np.fmax(c1, c2)
	try:
		erk_val = fuzz.defuzz(erk, c_com, 'centroid')
	except AssertionError as e:
		erk_val = 0
	return erk_val

def calculate_pi3k(egfr_mfs, erk_mfs, pi3k_mfs, pi3k, egfr_index, erk_index):
	a1 =  egfr_mfs[0][egfr_index]
	c1 = np.fmin(a1, pi3k_mfs[0])
	a2 = erk_mfs[1][erk_index]
	c2 = np.fmin(a2, pi3k_mfs[0])
	c_com = np.fmax(c1, c2)
	a3_1 = egfr_mfs[1][egfr_index]
	a3_2 = erk_mfs[0][erk_index]
	a3 = min(a3_1, a3_2)
	c3 = np.fmin(a3, pi3k_mfs[1])
	c_com = np.fmax(c_com, c3)
	try:
		pi3k_val = fuzz.defuzz(pi3k, c_com, 'centroid')
	except AssertionError as e:
		pi3k_val = 0
	return pi3k_val

def calaculate_akt(pi3k_mfs, akt_mfs, akt, pi3k_index):
	a1 = pi3k_mfs[0][pi3k_index]
	c1 = np.fmin(a1, akt_mfs[0])
	a2 = pi3k_mfs[1][pi3k_index]
	c2 = np.fmin(a2, akt_mfs[1])
	c_com = np.fmax(c1, c2)
	try:
		akt_val = fuzz.defuzz(akt, c_com, 'centroid')
	except AssertionError as e:
		akt_val = 0
	return akt_val

def rules(present_values, initial_values, time_indexes, time_length):
	
	y = np.copy(present_values)
	egf_index = np.searchsorted(initial_values[0], present_values[0], 'left')
	egfr_index = np.searchsorted(initial_values[2], present_values[1], 'left')
	raf_index = np.searchsorted(initial_values[3], present_values[2], 'left')
	erk_index = np.searchsorted(initial_values[4], present_values[3], 'left')
	pi3k_index = np.searchsorted(initial_values[5], present_values[4], 'left')

	egfr_mfs = ()
	raf_mfs = ()
	erk_mfs = ()
	pi3k_mfs = ()

	if(egf_index > 0):
		egfr_mfs = calculate_egfr_mfs(initial_values[-1], initial_values[0], initial_values[2], egf_index, time_length[0])

	raf_mfs = calculate_raf_mfs(initial_values[2], initial_values[3])
	pi3k_mfs = calculate_pi3k_mfs(initial_values[2], initial_values[4], initial_values[5])
	erk_mfs = calculate_erk_mfs(initial_values[3], initial_values[4])
	akt_mfs = calculate_akt_mfs(initial_values[5], initial_values[6])
	
	if(present_values[0] == 0):
		y[1] = 0		
	else:
		y[1] = calculate_egfr((egfr_mfs[0], egfr_mfs[1]), (egfr_mfs[2], egfr_mfs[3]), initial_values[2], time_indexes[0])
		
	if(present_values[1] == 0):
		y[2] = 0
		y[4] = 0
	else:
		y[2] = calculate_raf((raf_mfs[0], raf_mfs[1]), (raf_mfs[2], raf_mfs[3]), initial_values[3], egfr_index)
		y[4] = calculate_pi3k((pi3k_mfs[0][0], pi3k_mfs[0][1]), (pi3k_mfs[1][0], pi3k_mfs[1][1]), (pi3k_mfs[2][0], pi3k_mfs[2][1]),  initial_values[5], egfr_index, erk_index)

	if(present_values[2] == 0):
		y[3] = 0
	else:
		y[3] = calculate_erk((erk_mfs[0], erk_mfs[1]), (erk_mfs[2], erk_mfs[3]), initial_values[4], raf_index)

	if(present_values[4] == 0):
		y[5] = 0
	else:
		y[5] = calaculate_akt((akt_mfs[0][0], akt_mfs[0][1]), (akt_mfs[1][0], akt_mfs[1][1]), initial_values[6], pi3k_index)

	not_updated = []
	for i in xrange(len(y)):
		if y[i] == present_values[i]:
			not_updated.append(i)

	if 0 in not_updated:
		time_indexes[0] += 1
	else:
		time_indexes[0] = 1

	return y , time_indexes

def main():
	y = np.array([0.9, 0, 0, 0, 0, 0])
	size = 100
	time_length = [10, 1.0, 1.0]
	time = np.linspace(0, time_length[0], size*time_length[0] + 1)
	egf = np.linspace(0, 1, size + 1)
	egfr = np.linspace(0, y[0], size + 1)
	raf = egfr
	pi3k = egfr
	akt = egfr
	erk = np.linspace(0, y[0], size + 1)
	egfr_out = np.zeros(time.size)
	time_indexes = [1, 1]
	new_raf_out = np.zeros(time_length[0]*size + 1)
	new_erk_out = np.zeros(time_length[0]*size + 1)
	new_pi3k_out = np.zeros(time_length[0]*size + 1)
	new_akt_out = np.zeros(time_length[0]*size + 1)
	'''for i in xrange(1,size+1):
		if i == 0:
			for t in xrange(1, time.size + 1):
				egfr_out[t] = 0
			continue
		egfr_mfs = calculate_egfr_mfs(time, egf, egfr, i, time_length)
		#raf_mfs = calculate_raf_mfs(time_raf, egfr, raf, i, 1)
		for t in xrange(1,time.size):
			if t == 0:
				egfr_out[t] = 0
				raf_out[t] = 0
				continue
			egfr_out[t] = calculate_egfr((egfr_mfs[0], egfr_mfs[1]), (egfr_mfs[2], egfr_mfs[3]), egfr, t)
			#raf_out[t] = calculate_egfr((raf_mfs[0], raf_mfs[1]), (raf_mfs[2], raf_mfs[3]), raf, t)
		#time = time_raf
		if i%10 == 0:
			#pass
			plt.plot(time, egfr_out)
			plt.xlabel('Time')
			plt.ylabel('egfr')
			#time_length = 1
			plt.axis([-0.01,time_length + .1,-0.01,1.1])
	plt.show()'''
	initial_vals = (egf, egf, egfr, raf, erk, pi3k, akt, time)
	y.resize(1, y.size)
	time_raf_delay = time_length[1]*size
	time_erk_delay = time_length[2]*size
	time_pi3k_delay = time_raf_delay
	time_akt_delay = time_erk_delay

	if((time_raf_delay%2 == 0 or (time_raf_delay + 1)%2 == 0) and time_erk_delay%2 == 0 or (time_erk_delay+1)%2 == 0):
		for i in xrange(1, time.size):
			temp, time_indexes = rules(y[i - 1], initial_vals, time_indexes, time_length)
			y = np.vstack((y, temp))
					
			new_raf_out[i] = y[i, 2]
			new_erk_out[i] = y[i, 3]
			new_pi3k_out[i] = y[i, 4]
			new_akt_out[i] = y[i, 5]

			if( i > time_raf_delay ):
				y[i, 2] = new_raf_out[i - time_raf_delay]
			else:
				y[i, 2] = 0

			if(i > time_pi3k_delay):
				y[i ,4] = new_pi3k_out[i - time_pi3k_delay]
			else:
				y[i, 4] = 0

			if( i > time_erk_delay):
				y[i, 3] = new_erk_out[i - time_erk_delay]
			else:
				y[i, 3] = 0

			if( i > time_akt_delay):
				y[i, 5] = new_akt_out[i - time_akt_delay]
			else:
				y[i, 5] = 0

			#print y[i]
		plt.plot(time, y[:,1], time, y[:,5], time, y[:,4])
		plt.xlabel('Time')
		plt.ylabel('egfr')
		plt.axis([-0.1,time_length[0] + 0.01, -0.01, 1.01])
		plt.show()

	else:
		print "Not a valid time_raf_delay or time_erk_delay"

if(__name__ == '__main__'):
	main()