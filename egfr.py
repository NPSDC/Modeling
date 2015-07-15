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

def rules(present_values, initial_values, time_indexes, time_length):
	
	y = np.copy(present_values)
	egf_index = np.searchsorted(initial_values[0], present_values[0], 'left')
	egfr_index = np.searchsorted(initial_values[2], present_values[1], 'left')

	egfr_mfs = ()
	raf_mfs = ()

	if(egf_index > 0):
		egfr_mfs = calculate_egfr_mfs(initial_values[-1], initial_values[0], initial_values[2], egf_index, time_length[0])
	if(egfr_index > 0):
		raf_mfs = calculate_raf_mfs(initial_values[2], initial_values[3])

	
	if(present_values[0] == 0):
		y[1] = 0
	else:
		y[1] = calculate_egfr((egfr_mfs[0], egfr_mfs[1]), (egfr_mfs[2], egfr_mfs[3]), initial_values[2], time_indexes[0])

	if(present_values[1] == 0):
		y[2] = 0
	else:
		y[2] = calculate_raf((raf_mfs[0], raf_mfs[1]), (raf_mfs[2], raf_mfs[3]), initial_values[3], egfr_index)

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
	y = np.array([0.9, 0, 0])
	size = 100
	time_length = [10, 0.2]
	time = np.linspace(0, time_length[0], size*time_length[0] + 1)
	egf = np.linspace(0, 1, size + 1)
	egfr = np.linspace(0, y[0], size + 1)
	raf = egfr
	egfr_out = np.zeros(time.size )
	time_indexes = [1, 1]
	new_raf_out = np.zeros(time_length[0]*size + 1)
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
	initial_vals = (egf, egf, egfr, raf, time)
	y.resize(1, 3)
	time_raf_delay = time_length[1]*size

	if(time_raf_delay%2 == 0 or (time_raf_delay + 1)%2 == 0):
		for i in xrange(1, time.size):
			temp, time_indexes = rules(y[i - 1], initial_vals, time_indexes, time_length)
			y = np.vstack((y, temp))
			new_raf_out[i] = y[i, 2]
			if( i > time_raf_delay ):
				y[i, 2] = new_raf_out[i - time_raf_delay]
			else:
				y[i, 2] = 0
		plt.plot(time,y[:,1], time, y[:,2])
		plt.xlabel('Time')
		plt.ylabel('egfr')
		plt.axis([-0.01,time_length[0] + 0.01, -0.01, 1.01])		
		plt.show()

	else:
		print "Not a valid time_raf_delay"

if(__name__ == '__main__'):
	main()