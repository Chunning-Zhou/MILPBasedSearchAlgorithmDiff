from numpy import *
import numpy as np 
import function
import copy
from itertools import combinations
class class_present:
	def __init__(self):
		self.name = "present"
		self.structure = "sp"
		self.oriented = "bit"
		self.block_size = 64
		self.sbox_size = 4
		self.num_of_sbox = 1
		self.nibble = 16
		self.sbox = [[0xc, 0x5, 0x6, 0xb, 0x9, 0x0, 0xa, 0xd, 0x3, 0xe, 0xf, 0x8, 0x4, 0x7, 0x1, 0x2]] 
		self.pbox = [0, 16, 32, 48, 1, 17, 33, 49, 2, 18, 34,50, 3, 19, 35, 51,\
		4, 20, 36, 52, 5, 21, 37, 53, 6, 22, 38, 54, 7, 23, 39, 55,\
		8, 24, 40, 56, 9, 25, 41, 57, 10, 26, 42, 58, 11, 27, 43, 59,\
		12, 28, 44, 60, 13, 29, 45, 61, 14, 30, 46, 62, 15, 31, 47, 63]
		self.branch_num_of_sbox = 3
		self.var_and_num_LBAS = {"x":[0,self.block_size],"d":[self.nibble,0],"A":[self.nibble,0]}
		self.var_and_num_AS = {"x":[self.block_size,self.block_size],"d":[self.nibble,0], "A":[self.nibble,0]}
		self.var_and_num_DC = {"x":[self.block_size,self.block_size],"d":[self.nibble,0], "p":[self.nibble*2,0]} #,"d":[self.nibble,0]

		# searching
		self.diff_all = self.get_diff_all()
		self.num_of_diff_all = [240,27000]
		self.min_weight_of_sbox = 2
		self.round_of_first_region2 = 6

	def get_search_round(self,goal):
		if goal == "AS":
			return 31
		elif goal == "DC":
			return 18

	def get_max_num_of_round_to_solve(self,goal):
		return 15,0


	def get_diff_all(self):
		diff = []
		position_1 = list(combinations(np.arange(1,self.nibble+1),1))
		for i in position_1:
			diff_round_1 = [0 for l in range(self.nibble)]
			j = list(i)[0]
			for k in range(1,pow(2,self.sbox_size)):
				diff_round_1[j-1] = k
				diff.append(copy.deepcopy(diff_round_1))
		position_2 = list(combinations(np.arange(1,self.nibble+1),2))
		for i in position_2:
			diff_round_2 = [0 for l in range(self.nibble)]
			j1 = list(i)[0]
			j2 = list(i)[1]
			for k1 in range(1,pow(2,self.sbox_size)):
				for k2 in range(1,pow(2,self.sbox_size)):
					diff_round_2[j1-1] = k1
					diff_round_2[j2-1] = k2
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
		for n in range(self.nibble):
			var["x"] = state[self.sbox_size*n:self.sbox_size*(n+1)]
			var["y"] = state_through_sbox[self.sbox_size*n:self.sbox_size*(n+1)]
			var["d"] = self.nibble*(r-1)+n+1
			if model_goal == "LBAS":
				var["A"] = self.nibble*(r-1)+n+1
				function.diff_propagation_of_sbox_LBAS_bit(self.sbox_size,self.branch_num_of_sbox,model_filename,var)
			else:
				for i in range(self.sbox_size):
					with open(model_filename, "a") as f:
						f.write("x%d + "%(var["x"][i]))
				for i in range(self.sbox_size-1):
					with open(model_filename, "a") as f:			
						f.write("x%d + "%(var["y"][i]))
				with open(model_filename, "a") as f:	
					f.write("x%d - %d d%d >= 0\n"%(var["y"][self.sbox_size-1],self.branch_num_of_sbox,var["d"]))
				for i in range(self.sbox_size):
					with open(model_filename, "a") as f:	
						f.write("d%d - x%d >= 0\n"%(var["d"],var["x"][i]))
				for i in range(self.sbox_size):
					with open(model_filename, "a") as f:	
						f.write("d%d - x%d >= 0\n"%(var["d"],var["y"][i]))
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
					f.write("x%d = %d \n"%(state[self.sbox_size*i+j],int(value_bin[j])))