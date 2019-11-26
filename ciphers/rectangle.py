from numpy import *
import numpy as np 
import function
import copy
class class_rectangle:
	def __init__(self):
		self.name = "rectangle"
		self.structure = "sp"
		self.oriented = "bit"
		self.block_size = 64
		self.sbox_size = 4
		self.num_of_sbox = 1
		self.nibble = 16
		self.sbox = [[0x6,0x5,0xc,0xa,0x1,0xe,0x7,0x9,0xb,0x0,0x3,0xd,0x8,0xf,0x4,0x2]]
		self.pbox = [13,14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,\
		28,29,30,31,16,17,18,19,20,21,22,23,24,25,26,27,\
		33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,32,\
		48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]
		self.branch_num_of_sbox = function.hamming_weight(self.name,self.sbox[0],self.sbox_size) #2
		self.var_and_num_LBAS = {"x":[0,self.block_size],"d":[self.nibble,0],"A":[self.nibble,0]}
		self.var_and_num_AS = {"x":[self.block_size,self.block_size],"A":[self.nibble,0]}
		self.var_and_num_DC = {"x":[self.block_size,self.block_size],"p":[self.nibble*2,0]}

		# searching
		self.diff_all = self.get_diff_all()
		self.num_of_diff_all = [15,3375]
		self.min_weight_of_sbox = 2
		self.max_round_of_solve_LBAS = 0
		self.round_of_first_region2 = 28
	
	def get_search_round(self,goal):
		if goal == "AS":
			return 17
		elif goal == "DC":
			return 16

	def get_max_num_of_round_to_solve(self,goal):
		if goal == "AS":
			max_num_of_round_to_solve = 6
		elif goal == "DC":
			max_num_of_round_to_solve = 5
		return 0,max_num_of_round_to_solve
	
	def get_diff_all(self):
		diff = []
		diff_round_1 = [0 for i in range(self.nibble)]
		for k in range(1,pow(2,self.sbox_size)):
			diff_round_1[0] = k
			diff.append(copy.deepcopy(diff_round_1))
		for i in range(2,self.nibble+1):
			diff_round_2 = [0 for j in range(self.nibble)]
			for k1 in range(1,pow(2,self.sbox_size)):
				for k2 in range(1,pow(2,self.sbox_size)):
					diff_round_2[0] = k1
					diff_round_2[i-1] = k2
					diff.append(copy.deepcopy(diff_round_2))
		return diff

	def obj_fun(self,model_goal,model_filename,r1,r2): 
		if model_goal == "LBAS" or model_goal == "AS": 
			for i in range ((r1-1)*self.nibble+1,r2*self.nibble):
				with open(model_filename, "a") as f:
					f.write("A%d + "%(i))
			with open(model_filename, "a") as f:
				f.write("A%d"%(r2*self.nibble))
		elif model_goal == "DC":
			for i in range ((r1-1)*self.nibble*2+1,r2*self.nibble*2):
				if i%2 == 1:
					with open(model_filename, "a") as f:
						f.write("1 p%d + "%(i))
				elif i%2 == 0:
					with open(model_filename, "a") as f:
						f.write("2 p%d + "%(i))
			with open(model_filename, "a") as f:
				f.write("%d p%d"%(2, r2*self.nibble * 2))
	
	def gen_input_state(self): 
		state = np.arange(1,self.block_size+1)
		return state
	
	def get_state_through_sbox(self,r):
		state = np.arange(self.block_size*r+1,self.block_size*(r+1)+1)
		return state
	
	def diff_propagation_of_sbox(self,model_goal,model_filename,r,state,state_through_sbox):
		var = {}
		temp = [0 for i in range(self.sbox_size)]
		for n in range(self.nibble):
			var["x"] = [state[n],state[n+self.nibble*1],state[n+self.nibble*2],state[n+self.nibble*3]]
			var["y"] = [state_through_sbox[n],state_through_sbox[n+self.nibble*1],state_through_sbox[n+self.nibble*2],state_through_sbox[n+self.nibble*3]]
			if model_goal == "AS":
				var["A"] = self.nibble*(r-1)+n+1
				ine = np.loadtxt("txt/"+self.name + "/AS/1th_sbox_reduce_ine_greedy.txt",int)
				function.diff_propagation_of_sbox_AS(self.sbox_size,model_filename,var,ine)
			elif model_goal == "DC":
				var["p"] = [2*n+(r-1)*2*self.nibble+1, 2*n+(r-1)*2*self.nibble+2] 
				ine = np.loadtxt("txt/"+self.name + "/DC/1th_sbox_reduce_ine_greedy.txt",int)
				function.diff_propagation_of_sbox_DC(self.sbox_size,self.var_and_num_DC["p"][0]/self.nibble,model_filename,var,ine)
	
	def get_state_through_per(self,state): 
		state1 = np.arange(1,self.block_size+1)
		for i in range(self.block_size):
			state1[self.block_size-1-self.pbox[self.block_size-1-i]] = state[i]
		return state1
	
	def lowbound_of_sbox(self,model_goal,model_filename,r1,r2,lb):
		if model_goal == "LBAS" or model_goal == "AS":
			for r in range(r1,r2+1):
				for i in range(1,self.nibble):
					with open(model_filename, "a") as f:
						f.write("A%d + "%(self.nibble*(r-1)+i))
				with open(model_filename, "a") as f:
					f.write("A%d >= %d \n"%(self.nibble*r,lb))	
		elif model_goal == "DC":
			for r in range(r1,r2+1):
				for i in range(1,self.nibble):
					with open(model_filename, "a") as f:
						f.write("p%d + "%(self.nibble*2*(r-1)+2*i))
				with open(model_filename, "a") as f:
					f.write("p%d >= %d\n"%(self.nibble*r*2,lb))	

	def fix_diff(self,model_goal,model_filename,diff_type,r,diff_value):
		if diff_type == "input_diff":
			state = self.get_state_through_sbox(r-1)
			if r > 1:
				state = self.get_state_through_per(state)
		elif diff_type == "output_diff":
			state = self.get_state_through_per(self.get_state_through_sbox(r))
		for i in range(self.nibble):
			value = diff_value[i]
			value_bin = bin(value).replace("0b","").zfill(self.sbox_size)
			for j in range(self.sbox_size):
				with open(model_filename, "a") as f:
					f.write("x%d = %d \n"%(state[self.nibble*j+i],int(value_bin[j])))