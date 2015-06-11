"""
	1) if egf is high then egfr is high
	2) if egf is low and hrg is high then egrf is mod
	3) if egf is low and hrg is low then egrf is low	
	
"""
import numpy as np
import skfuzzy as fuzz
import matplotlib.pyplot as plt

initial_cond = np.array([1,1,0])  ##stores the initial condition

#Universe functions
egf = np.linspace(0, 1, 30)
hrg = np.linspace(0, 1, 30)
time = np.linspace(0, 1, 30)
egfr = np.linspace(0, 2, 30)

#Membership functions for egf
egf_low = fuzz.trimf(egf, [0, 0, 1])
egf_high = fuzz.trimf(egf, [0, 1, 1])

#Membership functions for hrg
hrg_low = fuzz.trimf(hrg, [0, 0, 1])
hrg_high = fuzz.trimf(hrg, [0, 1, 1])

#Membership functions for time
time_low = fuzz.trimf(time, [0, 0, 1])
time_high = fuzz.trimf(time, [0, 1, 1])

#Membership functions for egfr
egfr_low = fuzz.trapmf(egfr, [0., 0., egfr[2], 2])
egfr_high = fuzz.trapmf(egfr, [0., egfr[-2], 2, 2])
egfr_mod = fuzz.trapmf(egfr, [0., egfr[14], egfr[16], 2])

#Combining antecedents  a stands for antecedent
R_a1 = np.fmin(egf_high, time_high) #1
R_a2 = np.fmin(egf_low, np.fmin(hrg_high, time_high)) #2
R_a3 = np.fmin(egf_low, hrg_low) #3
R_a4 = time_low #4

#Cumputing consequents
R_c1 = fuzz.relation_product(R_a1, egfr_high)
R_c2 = fuzz.relation_product(R_a2, egfr_mod)
R_c3 = fuzz.relation_product(R_a3, egfr_low) 
R_c4 = fuzz.relation_product(R_a4, egfr_low)
R_combined = np.fmax(np.fmax(R_c1, R_c2), np.fmax(R_c3, R_c4))

predicted_egfr = np.zeros_like(time)
'''plt.plot(egfr, R_combined[-1])
plt.show()
'''
for i in xrange(len(predicted_egfr)):
	predicted_egfr[i] = fuzz.defuzz(egfr, R_combined[i,:], 'centroid') 
#print fuzz.defuzz(egfr, R_combined[-1], 'centroid')	

plt.plot(time, predicted_egfr, 'go')
plt.xlabel('Time')
plt.ylabel('Egfr');
plt.show()