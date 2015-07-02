import numpy as np
import skfuzzy as fuzz	

'''The various membership functions are computed here'''

def eval_membership_function_egf(egf_values):
	'''Evaluates membership function for egf '''
	egf_high = fuzz.smf(egf_values, 0, 1)
	egf_low = fuzz.zmf(egf_values, 0, 1)
	return (egf_low, egf_high)

def eval_membership_function_hrg(hrg_values):
	'''Evaluates membership function for hrg '''
	hrg_high = fuzz.smf(hrg_values, 0 , 1)
	hrg_low = fuzz.zmf(hrg_values, 0, 1)
	return (hrg_low, hrg_high)

def eval_membership_function_egfr(egfr_values, positive_change_egfr, negative_change_egfr):
	'''Evaluates membership function for egfr '''
	egfr_high = fuzz.gaussmf(egfr_values,1, 0.1)
	egfr_low = fuzz.gaussmf(egfr_values,0, 0.3)
	egfr_high1 = fuzz.trapmf(egfr_values,(0.1,0.1,0.8,1.4))
	positive_change_egfr_high = fuzz.gaussmf(positive_change_egfr,1, 0.1)
	positive_change_egfr_low = fuzz.gaussmf(positive_change_egfr, 0, 0.1)
	negative_change_egfr_high = fuzz.gaussmf(negative_change_egfr, -1, 0.1)
	negative_change_egfr_low = fuzz.gaussmf(negative_change_egfr, 0, 0.1)
	return (egfr_low, egfr_high, positive_change_egfr_low, positive_change_egfr_high, negative_change_egfr_low, negative_change_egfr_high, egfr_high1)

def eval_membership_function_raf(raf_values, positive_change_raf, negative_change_raf):
	'''Evaluates membership function for raf '''
	raf_high = fuzz.gaussmf(raf_values, 1, 0.1)
	raf_low = fuzz.gaussmf(raf_values, 0, 0.1)
	raf_high1 = fuzz.trapmf(raf_values,(0.1,0.1,0.8,1.4))
	positive_change_raf_high = fuzz.gaussmf(positive_change_raf, 1, 0.01)
	positive_change_raf_low = fuzz.gaussmf(positive_change_raf, 0, 0.01)
	positive_change_raf_mid = fuzz.gaussmf(positive_change_raf, 0.5, 0.1)
	negative_change_raf_high = fuzz.gaussmf(negative_change_raf, -1, 0.01)
	negative_change_raf_low = fuzz.gaussmf(negative_change_raf, 0, 0.01)
	return (raf_low, raf_high, positive_change_raf_low, positive_change_raf_high, negative_change_raf_low, negative_change_raf_high, raf_high1)

def eval_membership_function_pi3k(pi3k_values, positive_change_pi3k, negative_change_pi3k):
	'''Evaluates membership function for pi3k '''
	pi3k_high = fuzz.gaussmf(pi3k_values, 1, 0.1)
	pi3k_low = fuzz.gaussmf(pi3k_values, 0, 0.1)
	pi3k_high1 = fuzz.trapmf(pi3k_values,(0.1,0.1,0.8,1.4))
	positive_change_pi3k_high = fuzz.gaussmf(positive_change_pi3k, 1, 0.1)
	positive_change_pi3k_low = fuzz.gaussmf(positive_change_pi3k, 0, 0.1)
	negative_change_pi3k_high = fuzz.gaussmf(negative_change_pi3k, -1, 0.1)
	negative_change_pi3k_low = fuzz.gaussmf(negative_change_pi3k, 0, 0.1)
	return (pi3k_low, pi3k_high, positive_change_pi3k_low, positive_change_pi3k_high, negative_change_pi3k_low, negative_change_pi3k_high, pi3k_high1)

def eval_membership_function_erk(erk_values, positive_change_erk, negative_change_erk):
	'''Evaluates membership function for erk '''
	erk_high = fuzz.gaussmf(erk_values, 1, 0.1)
	erk_low = fuzz.gaussmf(erk_values, 0, 0.1)
	erk_high1 = fuzz.trapmf(erk_values,(0.1,0.1,0.8,1.4))
	positive_change_erk_high = fuzz.gaussmf(positive_change_erk, 1, 0.1)
	positive_change_erk_low = fuzz.gaussmf(positive_change_erk, 0, 0.1)
	negative_change_erk_high = fuzz.gaussmf(negative_change_erk, -1, 0.1)
	negative_change_erk_low = fuzz.gaussmf(negative_change_erk, 0, 0.1)
	return (erk_low, erk_high, positive_change_erk_low, positive_change_erk_high, negative_change_erk_low, negative_change_erk_high, erk_high1)

def eval_membership_function_akt(akt_values, positive_change_akt, negative_change_akt):
	'''Evaluates membership function for akt '''
	akt_high = fuzz.gaussmf(akt_values, 1, 0.1)
	akt_low = fuzz.gaussmf(akt_values, 0, 0.3)
	akt_high1 = fuzz.trapmf(akt_values,(0.1,0.1,0.8,1.4))
	positive_change_akt_high = fuzz.gaussmf(positive_change_akt, 1, 0.1)
	positive_change_akt_low = fuzz.gaussmf(positive_change_akt, 0, 0.1)
	negative_change_akt_high = fuzz.gaussmf(negative_change_akt, -1, 0.1)
	negative_change_akt_low = fuzz.gaussmf(negative_change_akt, 0, 0.1)
	return (akt_low, akt_high, positive_change_akt_low, positive_change_akt_high, negative_change_akt_low, negative_change_akt_high, akt_high1)

def eval_membership_function_time(time_values):
	'''Evaluates membership function for time '''
	time_high = fuzz.smf(time_values, 0, 1)
	time_low = fuzz.zmf(time_values, 0, 1)
	return (time_low, time_high)

def eval_membership_functions(initial_values):
	'''Calls the mfs for all '''
	egf_mfs = eval_membership_function_egf(initial_values[0])
	hrg_mfs = eval_membership_function_hrg(initial_values[1])
	egfr_mfs = eval_membership_function_egfr(initial_values[2][0], initial_values[2][1], initial_values[2][2])
	raf_mfs = eval_membership_function_raf(initial_values[3][0], initial_values[3][1], initial_values[3][2])
	pi3k_mfs = eval_membership_function_pi3k(initial_values[4][0], initial_values[4][1], initial_values[4][2])
	erk_mfs = eval_membership_function_erk(initial_values[5][0], initial_values[5][1], initial_values[5][2])
	akt_mfs = eval_membership_function_akt(initial_values[6][0], initial_values[6][1], initial_values[6][2])
	time_mfs = eval_membership_function_time(initial_values[7])
	return (egf_mfs, hrg_mfs, egfr_mfs, raf_mfs, pi3k_mfs, erk_mfs, akt_mfs, time_mfs)