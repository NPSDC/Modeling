import numpy as np
import skfuzzy as fuzz
import matplotlib.pyplot as plt

'''def calculate_egfr(time, egf, size):
	time_low = fuzz.zmf(time, 0, 1 - egf + 0.1)
	time_high = fuzz.smf(time, 0,1 -  egf + 0.1)
	egfr = np.linspace(0, egf, size)
	egfr_out = np.zeros(size)
	egfr_high = fuzz.gaussmf(egfr, egf, 0.01)
	egfr_low = fuzz.gaussmf(egfr, 0, 0.01)

	for i in xrange(time.size):
		a1 = time_low[i] 
		c1 = np.fmin(a1, egfr_low)
		a2 = time_high[i]
		c2 = np.fmin(a2, egfr_high)
		c_com = np.fmax(c1, c2)
		egfr_out[i] = fuzz.defuzz(egfr, c_com, 'centroid')
	#print egfr_out	
	return egfr_out

def main():
	size = 100
	time = np.linspace(0, 1, size)
	egf = np.linspace(0, 1, 11)
	egfr = np.zeros(size)

	for i in egf:
		if i == 0:
			continue
		temp_egfr = calculate_egfr(time, i, size)
		plt.plot(time, temp_egfr)
		plt.xlabel('Time')
		plt.ylabel('egfr')
		plt.axis([-0.01,1.1,-0.01,1.1])
		egfr = np.vstack((egfr, temp_egfr))
	plt.show()'''

def calculate_egfr_mfs(time, egf, egf_index):
	size = egf.size
	egfr = np.linspace(0, egf[egf_index], size)
	time_low = fuzz.zmf(time, 0, time[egf.size - 1] - time[egf_index] + time[1] )
	time_high = fuzz.smf(time, 0, time[egf.size - 1] - time[egf_index] + time[1] )
	egfr_low = fuzz.gaussmf(egfr, 0, egf[egf_index]/20)
	egfr_high = transform_gaussmf(egfr, (egf[egf_index], egf[egf_index]/20), 1)
	return ((time_low, time_high, egfr_low, egfr_high), egfr)

def calculate_raf_mfs(initial_vals):
	#Raf
	egfr_low = fuzz.zmf(initial_vals[2], initial_vals[2][0], 0.95) #0,0.95
	egfr_high = fuzz.smf(initial_vals[2], 0.15, initial_vals[2][initial_vals[2].size - 1]) #0.25,0.9
	raf_low = fuzz.gaussmf(initial_vals[3], initial_vals[3][0], initial_vals[3][initial_vals[3].size - 1]/20)
	raf_high = fuzz.gaussmf(initial_vals[3], initial_vals[3][initial_vals[3].size - 1], initial_vals[3][initial_vals[3].size - 1]/20)
	raf = (egfr_low, egfr_high, raf_low, raf_high)
	return raf

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

def calculate_egfr(time_mfs, egfr_mfs, egfr, time_index):
	a1 = time_mfs[0][time_index]
	c1 = np.fmin(a1, egfr_mfs[0])
	a2 = time_mfs[1][time_index]
	c2 = np.fmin(a2, egfr_mfs[1])
	c_com = np.fmax(c1, c2)
	egfr_val = fuzz.defuzz(egfr, c_com, 'centroid')
	return egfr_val

def calculate_raf(egfr_mfs, raf_mf, raf, egfr_index):
	a1 =  egfr_mfs[0][egfr_index]
	c1 = np.fmin(a1, raf_mf[0])
	a2 = egfr_mfs[1][egfr_index]
	c2 = np.fmin(a2, raf_mf[1])
	c_com = np.fmax(c1, c2)
	raf_val = fuzz.defuzz(raf, c_com, 'centroid')
	return raf_val

def rules(input_values, initial_values, mfs, time_index,):
	y = np.copy(input_values)
	egfr_mfs = calculate_egfr_mfs(initial_values[-1],  initial_values[0], np.searchsorted(initial_values[0], input_values[0], 'left'))
	y[1] = calculate_egfr((egfr_mfs[0][0], egfr_mfs[0][1]), (egfr_mfs[0][2], egfr_mfs[0][3]), egfr_mfs[1], time_index )
	y[2] = calculate_raf((mfs[0], mfs[1]), (mfs[2] , mfs[3]), initial_values[3], np.searchsorted(initial_values[2], input_values[1], 'left'))
	return y

def main():
	y = np.array([0.9, 0, 0])
	size = 101
	time = np.linspace(0, 2, size)
	egf = np.linspace(0, 1, size)
	egfr = np.linspace(0, y[0], size)
	raf = egfr
	egfr_out = np.zeros(size )
	'''for i in xrange(1,size):
		if i == 0:
			for t in xrange(1, size):
				egfr_out[t] = 0
			continue
		egfr_mfs = calculate_egfr_mfs(time, egf, i)
		for t in xrange(1,size):
			if t == 0:
				egfr_out[t] = 0
				continue
			egfr_out[t] = calculate_egfr((egfr_mfs[0][0], egfr_mfs[0][1]), (egfr_mfs[0][2], egfr_mfs[0][3]), egfr_mfs[1], t)
		
		if i%10 == 0:
			plt.plot(time, egfr_out)
			plt.xlabel('Time')
			plt.ylabel('egfr')
			plt.axis([-0.01,1.1,-0.01,1.1])	
	plt.show()'''
	initial_vals = (egf, egf, egfr, raf, time)
	mfs = calculate_raf_mfs(initial_vals)
	y.resize(1, 3)

	for i in xrange(1, time.size):
		temp = rules(y[i - 1], initial_vals, mfs, i)
		y = np.vstack((y, temp))
	print y[:,1:3]
	plt.plot(time,y[:,1], time, y[:,2])
	plt.xlabel('Time')
	plt.ylabel('egfr')
	plt.axis([-0.01,2.1,-0.01,1.01])	
	plt.show()

if(__name__ == '__main__'):
	main()