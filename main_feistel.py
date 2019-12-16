import function
import parameter
from model import class_model
from numpy import *
import numpy as np 
import time
import copy 

global cipher
global goal
global search_round
global bestobj
global LB
global flag
global flag_lbas
global upperbound
global order_of_round
global R_front

def initArray():
	global cipher
	global search_round
	global bestobj
	global LB
	global flag
	global flag_lbas

	bestobj = np.zeros((search_round), dtype=int)
	LB = np.zeros((search_round,2,1+search_round,1+cipher.num_of_diff_pattern_all), dtype=int32)
	flag = np.zeros((search_round,2,1+search_round,1+cipher.num_of_diff_pattern_all), dtype=int32)
	flag_lbas = np.zeros((search_round,2,1+search_round,cipher.num_of_diff_pattern_all), dtype=int32)

	# methods 1 to 3

def Method_1(r,Na,i,DP):
	# split r rounds into the first r1 rounds and the last r-r1 rounds, r1 = 1,2,..,r-1. 
	global LB

	temp = 0
	if i == 0 and DP == 0:
		for r1 in range(1,r):
			temp = max(temp, LB[r1-1,Na-1,-1,-1]+LB[r-r1-1,Na-1,-1,-1])
	elif 1 <= i < r :
		if 1 < i < r:
			temp = LB[i-1,Na-1,i-1,DP-1] + LB[r-i+1-1,Na-1,0,DP-1]
		for r1 in range(1,i): 
			temp = max(temp, LB[r1-1,Na-1,-1,-1] + LB[r-r1-1,Na-1,i-r1-1,DP-1])
		for r1 in range(i+1,r):
			temp = max(temp, LB[r1-1,Na-1,i-1,DP-1]+LB[r-r1-1,Na-1,-1,-1])
	elif i == r:
		for r1 in range(1,r-1):
			temp = max(temp, LB[r1-1,Na-1,-1,-1]+LB[r-r1-1,Na-1,r-r1-1,DP-1])
	
	LB[r-1,Na-1,i-1,DP-1] = max(LB[r-1,Na-1,i-1,DP-1],temp)	

def Method_2(r,Na,i,DP):
	global LB
	global upperbound
	global order_of_round
	global R_front

	temp = 0
	R_front = order_of_round[np.where(order_of_round==i)[0][0]-1]

	if i < order_of_round[0]:
		# split 1: split r rounds into: (1) round 1 to i; (2) round (i+1) to R_front; (3) round (R_front+1) to r
		temp = function.get_value(i,LB[i-1,Na-1,-1,-1]) + LB[R_front-i-1,Na,-1,-1] + function.get_value(r-R_front,LB[r-R_front-1,Na-1,-1,-1])

		# split 2: split r rounds into: (1) round 1 to i; (2) round i to R_front; (3) round R_front+1 to r rounds. 
		temp = max(temp, function.get_value(i,LB[i-1,Na-1,i-1,DP-1]) + LB[R_front-i+1-1,Na,0,DP-1] + function.get_value(r-R_front,LB[r-R_front-1,Na-1,-1,-1]))
		
	elif i > order_of_round[0]:
		# split 1: split r rounds into: (1) round 1 to (R_front-1); (2) round R_front to (i-1); (3) round i to r
		temp = function.get_value(R_front-1,LB[R_front-2,Na-1,-1,-1]) + LB[i-1-R_front,Na,-1,-1] + function.get_value(r-i+1,LB[r-i+1-1,Na-1,-1,-1])
		
		# split 2: split r rounds into: (1) round 1 to R_front-1; (2) round R_front to i; (3) round i to r rounds. 
		temp = max(temp, function.get_value(R_front-1,LB[R_front-1-1,Na-1,-1,-1]) + LB[i-R_front,Na,i-R_front,DP-1] + function.get_value(r-i+1,LB[r-i+1-1,Na-1,0,DP-1]))
	
	if temp >= upperbound:
		LB[r-1,Na-1,i-1,DP-1] = max(LB[r-1,Na-1,i-1,DP-1],upperbound)	

def Method_3(r,Na,i,DP,Strategy,obj_compare):
	global cipher
	global goal
	global search_round
	global LB
	global flag
	global flag_lbas

	
	if Strategy == 'Rough':
		if cipher.oriented == "bit" and cipher.branch_num_of_sbox == 2:
			return 0
		else:
			# check if the model has been solved
			if flag_lbas[r-1,Na-1,i-1,DP-1] == 0:
				if i == r:
					const_param = {"model_goal":"LBAS","model_round":r,"const_diff":["diff_pattern_front",r,cipher.diff_pattern_all[DP-1]],"const_sbox":[[1,r-1,Na]],"obj_compare":0}
				elif 1<=i<r:
					const_param = {"model_goal":"LBAS","model_round":r,"const_diff":["diff_pattern_next",i,cipher.diff_pattern_all[DP-1]],"const_sbox":[[1,i-1,Na],[i+1,r,Na]],"obj_compare":0}
				elif i == 0 and DP == 0:
					const_param = {"model_goal":"LBAS","model_round":r,"const_diff":[],"const_sbox":[[1,r,Na]],"obj_compare":0}
				model = class_model(cipher,const_param)
				
				# assign a value to the lower bound array by the optimal objective value of the model
				if goal == "AS":
					LB[r-1,Na-1,i-1,DP-1] = max(LB[r-1,Na-1,i-1,DP-1],model.model_obj)
				elif goal == "DC":
					LB[r-1,Na-1,i-1,DP-1] = max(LB[r-1,Na-1,i-1,DP-1],model.model_obj*cipher.min_weight_of_sbox)
				
				# flag the model has been solved
				flag_lbas[r-1,Na-1,i-1,DP-1] = 1
				
				for r2 in range(r+1,search_round+1):
					if i == r:
						Method_1(r2,Na,r2,DP)
					elif 1<=i<=r:
						Method_1(r2,Na,i,DP)
				with open("result/"+cipher.name+"/"+goal+"/solved_LBAS_model.txt", "a") as f:
					f.write("%s: %d\n"%(const_param,model.model_obj))
				

	elif Strategy == 'Tightest':
		# check if the model has been solved
		if flag[r-1,Na-1,i-1,DP-1] == 0:
			if i == r:
				const_param = {"model_goal":goal,"model_round":r,"const_diff":["diff_pattern_front",r,cipher.diff_pattern_all[DP-1]],"const_sbox":[[1,r-1,Na]],"obj_compare":obj_compare}
			elif 1<=i<r:
				const_param = {"model_goal":goal,"model_round":r,"const_diff":["diff_pattern_next",i,cipher.diff_pattern_all[DP-1]],"const_sbox":[[1,i-1,Na],[i+1,r,Na]],"obj_compare":obj_compare}
			elif i == 0 and DP == 0:
				const_param = {"model_goal":goal,"model_round":r,"const_diff":[],"const_sbox":[[1,r,Na]],"obj_compare":obj_compare}
			model = class_model(cipher,const_param)
			
			# assign a value to the lower bound array by the optimal objective value of the model
			LB[r-1,Na-1,i-1,DP-1] = max(LB[r-1,Na-1,i-1,DP-1],model.model_obj)
			
			# flag the model has been solved
			flag[r-1,Na-1,i-1,DP-1] = 1
			
			for r2 in range(r+1,search_round+1):
				if i == r:
					Method_1(r2,Na,r2,DP)
				elif 1<=i<=r:
					Method_1(r2,Na,i,DP)
			with open("result/"+cipher.name+"/"+goal+"/solved_model.txt", "a") as f:
				f.write("%s: %d\n"%(const_param,model.model_obj))
			

def search(r):
	global cipher
	global goal
	global bestobj
	global LB
	global upperbound
	global order_of_round		

	if r == 1:
		# when searching the simplest 1 round, we don't use the split-Method and directly solve the model 
		initLBArray(r)
		Method_3(r,0,0,0,'Tightest',25600)
		bestobj[r-1] = LB[r-1,-1,-1,-1]


	else:
		# Step 1. Generating the currently best r-round differential characteristic and an upper bound of the minimum weight UpperBound by using Technique 1;
		upperbound = genUpperBound(r)
		
 		# Step 2. initializing the lower bound array LB
 		initLBArray(r)

 		# Step 3. search Subset0
 		searchSubset0(r)

 		# Step 4. search Subset1
 		searchSubset1(r)

 		# Step 5. storing the final upper bound, namely, the mimimum weight into the array 'bestobj'
		bestobj[r-1] = upperbound
	
	with open("result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_search_result.txt", "a") as f:
		f.write("From all of above, we obtain an optimal objective value:%d.\n"%(bestobj[r-1]))
	LB[r-1] = np.maximum(LB[r-1],bestobj[r-1])

def initLBArray(r):
	global cipher
	global goal
	global search_round
	global LB
	
	for Na in [0,1]:
		# (1) assigning a value by Method 1
		Method_1(r,Na,0,0)

		# (2) assigning a value by Method 3 and 2
		max_num_of_round_to_solve_LBAS, max_num_of_round_to_solve_AS_or_DC = cipher.get_max_num_of_round_to_solve(goal)
		if r <= max_num_of_round_to_solve_LBAS:
			Method_3(r,Na,0,0,'Rough',0)
		if r <= max_num_of_round_to_solve_AS_or_DC and Na > 0:
			Method_3(r,Na,0,0,'Tightest',0)

def genUpperBound(r):
	global cipher
	global goal

	const_param_1 = {"model_goal":goal,"model_round":r,"const_diff":[],"const_sbox":"get_upperbound_1","obj_compare":25600}
	model_1 = class_model(cipher,const_param_1)
	upperbound_1 = model_1.model_obj
	const_param_2 = {"model_goal":goal,"model_round":r,"const_diff":[],"const_sbox":"get_upperbound_2","obj_compare":upperbound_1}
	model_2 = class_model(cipher,const_param_2)
	upperbound_2 = model_2.model_obj
	with open("result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_search_result.txt","a") as f:
		f.write("initialized upperbound = min(%d,%d)=%d.\n"%(upperbound_1,upperbound_2,min(upperbound_1,upperbound_2)))
 	return min(upperbound_1,upperbound_2)
		
def searchSubset0(r):
	global cipher
	global goal
	global LB
	global upperbound
	global order_of_round
	order_of_round = function.get_search_round(r)
	Na = 0
 	for i in order_of_round:
 		for DP in range(1,1+cipher.num_of_diff_pattern_all):
 			LB[r-1,Na-1,i-1,DP-1] = max(LB[r-1,Na-1,i-1,DP-1],LB[r-1,Na-1,-1,-1])
 			Method_1(r,Na,i,DP)
 			Method_2(r,Na,i,DP)
 			if LB[r-1,0-1,i-1,DP-1] < upperbound: 
				UpdateLBSubset0(r,Na,i,DP)
				if LB[r-1,0-1,i-1,DP-1] < upperbound: 
					upperbound = LB[r-1,-1,i-1,DP-1]
					with open("result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_search_result.txt", "a") as f:
						f.write("*****************************************************************\n")
						f.write("(*_*) We update upperbound = LB[%d,%d,%d,%d] = %d. \n"%(r,0,i,DP,LB[r-1,-1,i-1,DP-1]))
						f.write("*****************************************************************\n")

def UpdateLBSubset0(r,Na,i,DP):
	global cipher
	global goal
	global search_round
	global LB
	global upperbound
	global order_of_round
	global R_front


	for Strategy in ['Rough', 'Tightest']:
		r1 = 2
		while r1 <= max(i,r-i+1) and LB[r-1,Na-1,i-1,DP-1] < upperbound:
			# case1: i = order_of_round[0]
			if i == order_of_round[0]:
				if r1 <= r-i+1:
					Method_3(r1,Na,1,DP,Strategy,0)
					Method_1(r,Na,i,DP)
				if r1 <= i and LB[r-1,Na-1,i-1,DP-1] < upperbound:
					Method_3(r1,Na,r1,DP,Strategy,0)
					Method_1(r,Na,i,DP)

			
			# case2: i < order_of_round[0]
			elif i < order_of_round[0]:
				if r1 <= i and LB[r-1,Na-1,i-1,DP-1] < upperbound:
					Method_3(r1,Na,r1,DP,Strategy,0)
					Method_1(r,Na,i,DP)
					Method_2(r,Na,i,DP)
				if r1 <= R_front-i+1 and LB[r-1,Na-1,i-1,DP-1] < upperbound:
					Method_3(r1,Na+1,1,DP,Strategy,0)
					Method_2(r,Na,i,DP)
				if r1 <= r-i+1 and LB[r-1,Na-1,i-1,DP-1] < upperbound:
					Method_3(r1,Na,1,DP,Strategy,0)
					Method_1(r,Na,i,DP)
					
			
			# case3: i > order_of_round[0]
			elif i > order_of_round[0]:
				if r1 <= r-i+1:
					Method_3(r1,Na,1,DP,Strategy,0)
					Method_1(r,Na,i,DP)
					Method_2(r,Na,i,DP)
				if r1 <= i-R_front+1 and LB[r-1,Na-1,i-1,DP-1] < upperbound:
					Method_3(r1,Na+1,r1,DP,Strategy,0)
					Method_2(r,Na,i,DP)
				if r1 <= i and LB[r-1,Na-1,i-1,DP-1] < upperbound:
					Method_3(r1,Na,r1,DP,Strategy,0)
					Method_1(r,Na,i,DP)

			r1 = r1 + 1
		
		if LB[r-1,Na-1,i-1,DP-1] < upperbound:
			Method_3(r,Na,i,DP,Strategy,upperbound)
		
		if Strategy == 'Rough' and i != order_of_round[0] and LB[r-1,Na-1,i-1,DP-1] < upperbound:
			estimateAddition1(r,Na,i,DP)

		
def estimateAddition1(r,Na,i,DP):
	global cipher
	global goal
	global LB
	global upperbound
	global order_of_round
	global R_front

	if cipher.oriented == 'bit' and cipher.branch_num_of_sbox == 2:
		return 0
	else:
		if i < order_of_round[0]:
			const_param = {"model_goal":"LBAS","model_round":r,"const_diff":["diff_pattern_next",i,cipher.diff_pattern_all[DP-1]],"const_sbox":[[1,i,Na],[i+1,R_front,Na+1],[R_front+1,r,Na]],"obj_compare":0}
		elif i > order_of_round[0]:
			const_param = {"model_goal":"LBAS","model_round":r,"const_diff":["diff_pattern_next",i,cipher.diff_pattern_all[DP-1]],"const_sbox":[[1,R_front-1,Na],[R_front,i-1,Na+1],[i,r,Na]],"obj_compare":0}
		
		model = class_model(cipher,const_param)
		with open("result/"+cipher.name+"/"+goal+"/solved_LBAS_model.txt", "a") as f:
			f.write("%s: %d\n"%(const_param,model.model_obj))
		if goal == "AS":
			result = model.model_obj
		elif goal == "DC":
			result = model.model_obj*cipher.min_weight_of_sbox
		if result >= upperbound:
			LB[r-1,Na-1,i-1,DP-1] = max(LB[r-1,Na-1,i-1,DP-1],upperbound)

def searchSubset1(r):
	global search_round
	global upperbound

	Na = 1

	if LB[r-1,Na-1,-1,-1] < upperbound:
		for Strategy in ['Rough', 'Tightest']:
			r1 = 1
			while r1 <= r and LB[r-1,Na-1,-1,-1] < upperbound:
				Method_3(r1,Na,0,0,Strategy)
				for r2 in range(r1 + 1, search_round+1):
					Method_1(r2,Na,0,0)
				r1 = r1 + 1
		upperbound = min(upperbound,LB[r-1,Na-1,-1,-1])



if __name__ == "__main__": 
	global cipher
	global goal
	global search_round
	global bestobj

	cipher = parameter.cipher
	goal = parameter.goal
	search_round = cipher.get_search_round(goal)
	time_all = []

	function.gen_filefolder(cipher.name,goal)
	with open("result/"+cipher.name+"/"+goal+"/"+"solved_LBAS_model.txt", "w") as f:
		f.write("models we solved:\n")
	with open("result/"+cipher.name+"/"+goal+"/"+"solved_model.txt", "w") as f:
		f.write("models we solved:\n")


	initArray()
	# search process from 1 to search_round rounds basing on Algorithm 1
	for r in range(1,search_round+1):
		time_start = time.time()
		with open("result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_search_result.txt", "w") as f:
			f.write("search process:\n")
		with open("result/"+cipher.name+"/"+goal+"/"+"solved_LBAS_model.txt", "a") as f:
			f.write("\n*****************************************************************************\nwhen searching the %d-round cipher, models we solved:\n"%(r))
		with open("result/"+cipher.name+"/"+goal+"/"+"solved_model.txt", "a") as f:
			f.write("\n*****************************************************************************\nwhen searching the %d-round cipher, models we solved:\n"%(r))
		search(r)
		function.get_trail_feistel(cipher,goal,r)
		time_end = time.time()
		timespend = time_end - time_start
		time_all.append(timespend)
		with open("result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_search_result.txt", "a") as f:
			f.write("time is %d s.\n"%(timespend))

		with open("result/"+cipher.name+"/"+goal+"/"+"bestobj.txt","w") as f:
			f.write("bestobj = %s\n"%(str(bestobj)))
			f.write("time_all = %s\n"%(str(time_all)))
	