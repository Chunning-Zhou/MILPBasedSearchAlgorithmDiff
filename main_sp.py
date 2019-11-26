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


# initialize arrays which are used during the search
def initArray():
	global cipher
	global search_round
	global bestobj
	global LB
	global flag
	global flag_lbas

	bestobj = np.zeros((search_round), dtype=int)
	LB = np.zeros((search_round,3,2+search_round,1+cipher.num_of_diff_all[0]+cipher.num_of_diff_all[1]), dtype=int32)
	flag = np.zeros((search_round,3,2+search_round,cipher.num_of_diff_all[0]+cipher.num_of_diff_all[1]), dtype=int32)
	flag_lbas = np.zeros((search_round,3,search_round,cipher.num_of_diff_all[0]+cipher.num_of_diff_all[1]), dtype=int32)

def Method_1(r,Na,i,D):
	# split r rounds into the first r1 rounds and the last r-r1 rounds, r1 = 1,2,..,r-1. 
	global LB

	temp = 0
	if i == 0 and D == 0:
		for r1 in range(1,r):
			temp = max(temp, LB[r1-1,Na-1,-1,-1]+LB[r-r1-1,Na-1,-1,-1])
	elif 1 <= i <= r :
		if i > 1:
			temp = LB[i-2,Na-1,i-1,D-1] + LB[r-i,Na-1,0,D-1]
		for r1 in range(1,i): 
			temp = max(temp, LB[r1-1,Na-1,-1,-1] + LB[r-r1-1,Na-1,i-r1-1,D-1])
		for r1 in range(i,r):
			temp = max(temp, LB[r1-1,Na-1,i-1,D-1]+LB[r-r1-1,Na-1,-1,-1])
	elif i == r+1:
		for r1 in range(1,r):
			temp = max(temp, LB[r1-1,Na-1,-1,-1]+LB[r-r1-1,Na-1,r-r1,D-1])
	LB[r-1,Na-1,i-1,D-1] = max(LB[r-1,Na-1,i-1,D-1],temp)	

def Method_2(r,Na,i,D):
	global LB
	global upperbound
	global order_of_round
	global R_front

	temp = 0
	R_front = order_of_round[np.where(order_of_round==i)[0][0]-1]
	
	if i < order_of_round[0]:
		# split 1: split r rounds into: (1) round 1 to i; (2) round (i+1) to R_front; (3) round (R_front+1) to r
		temp = function.get_value(i,LB[i-1,Na-1,-1,-1]) + LB[R_front-i-1,Na,-1,-1] + function.get_value(r-R_front,LB[r-R_front-1,Na-1,-1,-1])
		# split 2: split r rounds into (1) 1 to i-1; (2)i to R_front; (3) (R_front+1) to r
		temp = max(temp, function.get_value(i-1,LB[i-1-1,Na-1,i-1,D-1]) + LB[R_front-i+1-1,Na,0,D-1] + function.get_value(r-R_front,LB[r-R_front-1,Na-1,-1,-1]))
		
	elif i > order_of_round[0]:
		# split 1: split r rounds into: (1) round 1 to (R_front-1); (2) round R_front to (i-1); (3) round i to r
		temp = function.get_value(R_front-1,LB[R_front-1-1,Na-1,-1,-1]) + LB[i-1-R_front+1-1,Na,-1,-1] + function.get_value(r-i+1,LB[r-i+1-1,Na-1,-1,-1])
		# split 2: split r rounds into: (1) 1 to (R_front-1); (2)R_front to i-1; (3) i to r
		temp = max(temp, function.get_value(R_front-1,LB[R_front-1-1,Na-1,-1,-1]) + LB[i-1-R_front,Na,i-1-R_front+1,D-1] + function.get_value(r-i+1,LB[r-i+1-1,Na-1,0,D-1]))
	
	if temp >= upperbound:
		LB[r-1,Na-1,i-1,D-1] = max(LB[r-1,Na-1,i-1,D-1],upperbound)	

def Method_3(r,Na,i,D,Strategy,obj_compare):
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
			if flag_lbas[r-1,Na-1,i-1,D-1] == 0:
				
				# add constraints into the model 
				if i == r + 1:
					const_param = {"model_goal":"LBAS","model_round":r,"const_diff":["output_diff",r,cipher.diff_all[D-1]],"const_sbox":[[1,r,Na]],"obj_compare":0}
				elif i == 1:
					const_param = {"model_goal":"LBAS","model_round":r,"const_diff":["input_diff",i,cipher.diff_all[D-1]],"const_sbox":[[2,r,Na]],"obj_compare":0}
				elif i == 0 and D == 0:
					const_param = {"model_goal":"LBAS","model_round":r,"const_diff":[],"const_sbox":[[1,r,Na]],"obj_compare":0}
				model = class_model(cipher,const_param)
				
				# assign a value to the lower bound array by the optimal objective value of the model
				if goal == "AS":
					LB[r-1,Na-1,i-1,D-1] = max(LB[r-1,Na-1,i-1,D-1], model.model_obj)
				elif goal == "DC":
					LB[r-1,Na-1,i-1,D-1] = max(LB[r-1,Na-1,i-1,D-1], model.model_obj * cipher.min_weight_of_sbox)
				
				# flag the model has been solved
				flag_lbas[r-1,Na-1,i-1,D-1] = 1
				
				# uptaing LB aray by using Method_1
				for r2 in range(r+1,search_round+1):
					if i == r+1:
						Method_1(r2,Na,r2+1,D)
					elif 1<=i<=r:
						Method_1(r2,Na,i,D)
				
				with open("result/"+cipher.name+"/"+goal+"/solved_LBAS_model.txt", "a") as f:
					f.write("%s: %d\n"%(const_param,model.model_obj))
				

	elif Strategy == 'Tightest':
		# check if the model has been solved
		if flag[r-1,Na-1,i-1,D-1] == 0:
			if i == r + 1:
				const_param = {"model_goal":goal,"model_round":r,"const_diff":["output_diff",r,cipher.diff_all[D-1]],"const_sbox":[[1,r,Na]],"obj_compare":obj_compare}
			elif i == 1:
				const_param = {"model_goal":goal,"model_round":r,"const_diff":["input_diff",i,cipher.diff_all[D-1]],"const_sbox":[[i+1,r,Na]],"obj_compare":obj_compare}
			elif i == 0 and D == 0 and Na == 1:
				const_param = {"model_goal":goal,"model_round":r,"const_diff":[],"const_sbox":[],"obj_compare":25600}
			elif i == 0 and D == 0 and Na > 1:
				const_param = {"model_goal":goal,"model_round":r,"const_diff":[],"const_sbox":[[1,r,Na]],"obj_compare":0}
			model = class_model(cipher,const_param)
			
			# assign a value to the lower bound array by the optimal objective value of the model
			LB[r-1,Na-1,i-1,D-1] = model.model_obj
			
			# flag the model has been solved
			flag[r-1,Na-1,i-1,D-1] = 1		
			
			# uptaing LB aray by using Method_1
			for r2 in range(r+1,search_round+1):
				if i == r+1:
					Method_1(r2,Na,r2+1,D)
				elif 1<=i<=r:
					Method_1(r2,Na,i,D)
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
		Method_3(r,1,0,0,'Tightest',25600)
		bestobj[r-1] = LB[r-1,0,-1,-1]
	else:
		# Step 1. Generating the currently best r-round differential characteristic and an upper bound of the minimum weight UpperBound by using Technique 1;
		upperbound = genUpperBound(r)
		
 		# Step 2. initializing the lower bound array LB
 		initLBArray(r)

 		# Step 3. search Subset1 and Subset2
 		searchSubset12(r)

 		# Step 4. search Subset3
 		searchSubset3(r)

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
	
	for Na in [1,2,3]:
		# (1) assigning a value by Method 1
		Method_1(r,Na,0,0)

		# (2) assigning a value by Method 3 and 2
		max_num_of_round_to_solve_LBAS, max_num_of_round_to_solve_AS_or_DC = cipher.get_max_num_of_round_to_solve(goal)
		if r <= max_num_of_round_to_solve_LBAS:
			Method_3(r,Na,0,0,'Rough',0)
		if r <= max_num_of_round_to_solve_AS_or_DC and Na > 1:
			Method_3(r,Na,0,0,'Tightest',0)
			

def genUpperBound(r):
	global cipher
	global goal

	const_param_1 = {"model_goal":goal,"model_round":r,"const_diff":[],"const_sbox":"get_upperbound_1","obj_compare":25600}
	model = class_model(cipher,const_param_1)
	upperbound_1 = model.model_obj
	const_param_2 = {"model_goal":goal,"model_round":r,"const_diff":[],"const_sbox":"get_upperbound_2","obj_compare":upperbound_1}
	model = class_model(cipher,const_param_2)
	upperbound_2 = model.model_obj
	with open("result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_search_result.txt","a") as f:
		f.write("initialized upperbound = min(%d,%d)=%d.\n"%(upperbound_1,upperbound_2,min(upperbound_1,upperbound_2)))
 		
 	return min(upperbound_1,upperbound_2)
		

def searchSubset12(r):
	global cipher
	global goal
	global LB
	global upperbound
	global order_of_round

	order_of_Na = function.get_order_of_Na(cipher,r)
 	order_of_round = function.get_search_round(r)
 	for Na in order_of_Na:
 		for i in order_of_round:
 			for D in range(1+cipher.num_of_diff_all[0]*(Na-1),1+cipher.num_of_diff_all[0]*(Na-1)+cipher.num_of_diff_all[Na-1]):
 				LB[r-1,Na-1,i-1,D-1] = max(LB[r-1,Na-1,i-1,D-1],LB[r-1,Na-1,-1,-1])
 				Method_1(r,Na,i,D)
 				Method_2(r,Na,i,D)
 				if LB[r-1,Na-1,i-1,D-1] < upperbound: 
					UpdateLBSubset12(r,Na,i,D)
					if LB[r-1,Na-1,i-1,D-1] < upperbound: 
						upperbound = LB[r-1,Na-1,i-1,D-1]
						function.get_var_from_two_submodels(cipher,goal,r,Na,i,cipher.diff_all[D-1])
						with open("result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_search_result.txt", "a") as f:
							f.write("*****************************************************************\n")
							f.write("(*_*) We update upperbound = LB[%d,%d,%d,%d] = %d. \n"%(r,Na,i,D,LB[r-1,Na-1,i-1,D-1]))
							f.write("*****************************************************************\n")

def UpdateLBSubset12(r,Na,i,D):
	global cipher
	global goal
	global search_round
	global LB
	global upperbound
	global order_of_round
	global R_front


	for Strategy in ['Rough', 'Tightest']:
		r1 = 1
		while r1 <= max(i-1,r-i+1) and LB[r-1,Na-1,i-1,D-1] < upperbound:
			# case1: i = order_of_round[0]
			if i == order_of_round[0]:
				if r1 <= r-i+1:
					Method_3(r1,Na,1,D,Strategy,25600)
					Method_1(r,Na,i,D)
				if r1 <= i-1 and LB[r-1,Na-1,i-1,D-1] < upperbound:
					Method_3(r1,Na,r1+1,D,Strategy,25600)
					Method_1(r,Na,i,D)
			
			# case2: i < order_of_round[0]
			elif i < order_of_round[0]:
				if r1 <= i-1 and LB[r-1,Na-1,i-1,D-1] < upperbound:
					Method_3(r1,Na,r1+1,D,Strategy,25600)
					Method_1(r,Na,i,D)
					Method_2(r,Na,i,D)
				if r1 <= R_front-i+1 and LB[r-1,Na-1,i-1,D-1] < upperbound:
					Method_3(r1,Na+1,1,D,Strategy,0)
					Method_2(r,Na,i,D)
				if Strategy == 'Tightest' and r1 == R_front-i+1 and r-i+1 >= 8 and Na == 1 and LB[r-1,Na-1,i-1,D-1] < upperbound:
					estimateAddition2(r,Na,i,D)
				if r1 <= r-i+1 and LB[r-1,Na-1,i-1,D-1] < upperbound:
					Method_3(r1,Na,1,D,Strategy,25600)
					Method_1(r,Na,i,D)
					
			
			# case3: i > order_of_round[0]
			elif i > order_of_round[0]:
				if r1 <= r-i+1:
					Method_3(r1,Na,1,D,Strategy,25600)
					Method_1(r,Na,i,D)
					Method_2(r,Na,i,D)
				if r1 <= i-1-R_front+1 and LB[r-1,Na-1,i-1,D-1] < upperbound:
					Method_3(r1,Na+1,r1+1,D,Strategy,0)
					Method_2(r,Na,i,D)
				if Strategy == 'Tightest' and r1 == i-1-R_front+1 and i-1 >= 8 and Na == 1 and LB[r-1,Na-1,i-1,D-1] < upperbound:
					estimateAddition2(r,Na,i,D)
				if r1 <= i-1 and LB[r-1,Na-1,i-1,D-1] < upperbound:
					Method_3(r1,Na,r1+1,D,Strategy,25600)
					Method_1(r,Na,i,D)
			r1 = r1 + 1
		
		if Strategy == 'Rough' and i != order_of_round[0] and LB[r-1,Na-1,i-1,D-1] < upperbound:
			estimateAddition1(r,Na,i,D)
		
def estimateAddition1(r,Na,i,D):
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
			const_param = {"model_goal":"LBAS","model_round":r,"const_diff":["input_diff",i,cipher.diff_all[D-1]],"const_sbox":[[1,i,Na],[i+1,R_front,Na+1],[R_front+1,r,Na]],"obj_compare":0}
		elif i > order_of_round[0]:
			const_param = {"model_goal":"LBAS","model_round":r,"const_diff":["input_diff",i,cipher.diff_all[D-1]],"const_sbox":[[1,R_front-1,Na],[R_front,i-1,Na+1],[i,r,Na]],"obj_compare":0}
		model = class_model(cipher,const_param)
		with open("result/"+cipher.name+"/"+goal+"/solved_LBAS_model.txt", "a") as f:
			f.write("%s: %d\n"%(const_param,model.model_obj))
		if goal == "AS":
			result = model.model_obj
		elif goal == "DC":
			result = model.model_obj*cipher.min_weight_of_sbox
		if result >= upperbound:
			LB[r-1,Na-1,i-1,D-1] = max(LB[r-1,Na-1,i-1,D-1],upperbound)

def estimateAddition2(r,Na,i,D):
	global LB
	global upperbound
	global order_of_round
	global R_front

	
	if i < order_of_round[0] and R_front < r:
		# consider two cases of the number of active S-boxes at round (R_front+1)
		i2 = R_front+1

		# consider case 1: round (R_front+1) has exactly Na active s-boxes
		for D2 in range(1+cipher.num_of_diff_all[0]*(Na-1),1+cipher.num_of_diff_all[0]*(Na-1)+cipher.num_of_diff_all[Na-1]):
			Method_1(r,Na,i2,D2)
			case1_obj = LB[r-1,Na-1,i2-1,D2-1]
			r1 = 1 
			while r1 <= max(i2-i-1,r-i2+1) and case1_obj < upperbound:
				if r1 <= i2-i-1:
					Method_3(r1,Na+1,r1+1,D2,'Tightest',0)
					case1_obj = max(case1_obj,get_case1_obj(r,Na,i,D,i2,D2))
				if r1 <= r-i2+1:
					Method_3(r1,Na,1,D2,'Tightest',25600)
					case1_obj = max(case1_obj,get_case1_obj(r,Na,i,D,i2,D2))
				r1 = r1 + 1
			if case1_obj < upperbound:
				break
		# consider case 2: round (R_front+1) has greater than or equal to (Na+1) active s-boxes.
		if case1_obj >= upperbound:
			Method_3(i2-i+1,Na+1,1,D,'Tightest',0)
			case2_obj = function.get_value(i-1,LB[i-1-1,Na-1,i-1,D-1]) + LB[i2-i,Na,0,D-1] + function.get_value(r-i2,LB[r-i2-1,Na-1,-1,-1])
			if case2_obj >= upperbound:
				LB[r-1,Na-1,i-1,D-1] = upperbound
				
	elif i > order_of_round[0]:
		# consider two cases of the number of active S-boxes at round (R_front-1)
		i2 = R_front-1
		
		# consider case 1: round (R_front-1) has exactly Na active s-boxes
		for D2 in range(1+cipher.num_of_diff_all[0]*(Na-1),1+cipher.num_of_diff_all[0]*(Na-1)+cipher.num_of_diff_all[Na-1]):
			Method_1(r,Na,i2,D2)
			case1_obj = LB[r-1,Na-1,i2-1,D2-1]
			r1 = 1
			while r1 <= max(i2-1,i-1-i2+1) and case1_obj < upperbound:
				if r1 <= i2-1:
					Method_3(r1,Na,r1+1,D2,'Tightest',25600)
					case1_obj = max(case1_obj,get_case1_obj(r,Na,i2,D2,i,D))
				if r1 <= i-i2:
					Method_3(r1,Na+1,1,D2,'Tightest',0)
					case1_obj = max(case1_obj,get_case1_obj(r,Na,i2,D2,i,D))
				r1 = r1 + 1
			if case1_obj < upperbound:
				break
		
		# consider case 2: round (R_front+1) has greater than or equal to (Na+1) active s-boxes.
		if case1_obj >= upperbound: 
			Method_3(i-i2,Na+1,i-i2+1,D,'Tightest',0)
			case2_obj = function.get_value(i2-1,LB[i2-2,Na-1,-1,-1]) + LB[i-i2-1,Na,i-i2,D-1] + LB[r-i,Na-1,0,D-1]
			if case2_obj >= upperbound:
				LB[r-1,Na-1,i-1,D-1] = upperbound
				
def get_case1_obj(r,Na,r1,D1,r2,D2):
	obj_list = []
	for k in range(r1,r2):
		obj_list.append(LB[k-r1,Na,0,D1-1]+function.get_value(r2-1-k,LB[r2-1-k-1,Na,r2-1-k,D2-1]))
	temp = function.get_value(r1-1,LB[r1-2,Na-1,r1-1,D1-1]) + max(obj_list) + LB[r-r2,Na-1,0,D2-1]
	return temp		


def searchSubset3(r):
	global search_round
	global upperbound

	if LB[r-1,2,-1,-1] < upperbound:
		for Strategy in ['Rough', 'Tightest']:
			r1 = 1
			while r1 <= r and LB[r-1,2,-1,-1] < upperbound:
				Method_3(r1,3,0,0,Strategy)
				for r2 in range(r1 + 1, search_round+1):
					Method_1(r2,3,0,0)
				r1 = r1 + 1
		upperbound = min(upperbound,LB[r-1,2,-1,-1])


if __name__ == "__main__": 
	global cipher
	global goal
	global search_round
	global bestobj

	cipher = parameter.cipher
	goal = parameter.goal
	search_round = cipher.get_search_round(goal)
	time_all = []

	
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
		function.get_trail_sp(cipher,goal,r)
		time_end = time.time()
		timespend = time_end - time_start
		time_all.append(timespend)
		with open("result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_search_result.txt", "a") as f:
			f.write("time is %d s.\n"%(timespend))

		with open("result/"+cipher.name+"/"+goal+"/"+"bestobj.txt","w") as f:
			f.write("bestobj = %s\n"%(str(bestobj)))
			f.write("time_all = %s\n"%(str(time_all)))
		
		