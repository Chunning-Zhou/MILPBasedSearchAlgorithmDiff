from numpy import *
import numpy as np 
import function
import copy
from itertools import combinations

class class_lblock:
	def __init__(self):
		self.name = "lblock"
		self.structure = "feistel"
		self.oriented = "byte"
		self.block_size = 64
		self.sbox_size = 4
		self.num_of_sbox = 8
		self.nibble = 8
		self.sbox = [[13, 10, 15, 0, 14, 4, 9, 11, 2, 1, 8, 3, 7, 5, 12, 6],\
					[11, 9, 4, 14, 0, 15, 10, 13, 6, 12, 5, 7, 3, 8, 1, 2],\
					[2, 13, 11, 12, 15, 14, 0, 9, 7, 10, 6, 3, 1, 8, 4, 5],\
					[14, 5, 15, 0, 7, 2, 12, 13, 1, 8, 4, 9, 11, 10, 6, 3],\
					[7, 6, 8, 11, 0, 15, 3, 14, 9, 10, 12, 13, 5, 2, 4, 1],\
					[1, 14, 7, 12, 15, 13, 0, 6, 11, 5, 9, 3, 2, 4, 8, 10],\
					[4, 11, 14, 9, 15, 13, 0, 10, 7, 12, 5, 6, 2, 8, 1, 3],\
					[14, 9, 15, 0, 13, 4, 10, 11, 1, 2, 8, 3, 7, 6, 12, 5]]
		self.pbox_left_DC = [8,9,10,11,0,1,2,3,12,13,14,15,4,5,6,7,\
		24,25,26,27,16,17,18,19,28,29,30,31,20,21,22,23]
		self.pbox_right_DC = [8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,\
		25,26,27,28,29,30,31,0,1,2,3,4,5,6,7]
		self.pbox_left_LBAS = [2,0,3,1,6,4,7,5]
		self.pbox_right_LBAS = [2,3,4,5,6,7,0,1]
		self.var_and_num_LBAS = {"A":[self.block_size/8,self.block_size/8],"d":[self.block_size/8,-self.block_size/8]}
		self.var_and_num_AS = {"x":[self.block_size,self.block_size/2],"A":[self.nibble,0],"d":[self.block_size/2,-self.block_size/2]}
		self.var_and_num_DC = {"x":[self.block_size,self.block_size/2],"p":[self.nibble*2,0],"d":[self.block_size/2,-self.block_size/2]}

		# searching
		self.diff_pattern_all = self.get_diff_pattern_all()
		self.num_of_diff_pattern_all = 255
		self.min_weight_of_sbox = 2

	def get_search_round(self,goal):
		if goal == "AS":
			return 20
		elif goal == "DC":
			return 16

	def get_max_num_of_round_to_solve(self,goal):
		return 16,0

 	def get_diff_pattern_all(self):
		diff_pattern = []
		for i in range(1,self.nibble+1):
			position = list(combinations(np.arange(1,self.nibble+1),i))
			for j in position:
				diff_pattern_round = [0 for l in range(self.nibble)]
				for k in j:
					diff_pattern_round[k-1] = 1
				diff_pattern.append(copy.deepcopy(diff_pattern_round))
		return diff_pattern
		
	def obj_fun(self,model_goal,model_filename,r1,r2): 
		if model_goal == "LBAS":
			r1 = r1+1
			r2 = r2+1
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
	def gen_input_state(self,model_goal): 
		if model_goal == "LBAS":
			state = np.arange(1,self.nibble*2+1)
		elif model_goal == "AS" or model_goal == "DC":
			state = np.arange(1,self.block_size+1)
		return state
	def get_state_through_sbox(self,model_goal,r):
		if model_goal == "LBAS":
			state = np.arange(self.block_size*r+1,self.block_size*(r+1)+1)
		elif model_goal == "AS" or model_goal == "DC":
			state = np.arange(self.block_size*r+1,self.block_size*r+1+self.block_size/2)
		return state
	def diff_propagation_of_sbox(self,model_goal,model_filename,r,state,state_through_sbox):
		var = {}
		for n in range(self.nibble):
			var["x"] = state[self.sbox_size*n:self.sbox_size*(n+1)]
			var["y"] = state_through_sbox[self.sbox_size*n:self.sbox_size*(n+1)]
			if model_goal == "AS":
				var["A"] = self.nibble*(r-1)+n+1
				ine = np.loadtxt("txt/"+self.name + "/AS/"+str(n+1)+"th_sbox_reduce_ine_greedy.txt",int)
				function.diff_propagation_of_sbox_AS(self.sbox_size,model_filename,var,ine)
			elif model_goal == "DC":
				var["p"] = [2*n+(r-1)*2*self.nibble+1, 2*n+(r-1)*2*self.nibble+2] 
				ine = np.loadtxt("txt/"+self.name + "/DC/"+str(n+1)+"th_sbox_reduce_ine_greedy.txt",int)
				function.diff_propagation_of_sbox_DC(self.sbox_size,self.var_and_num_DC["p"][0]/self.nibble,model_filename,var,ine)
	def get_state_through_per_right(self,model_goal,state): 
		if model_goal == "LBAS":
			state1 = zeros(self.nibble, dtype = int16)
			for i in range(self.nibble):
				state1[self.nibble-1-self.pbox_right_LBAS[self.nibble-1-i]] = state[i]
		elif model_goal == "AS" or model_goal == "DC":
			state1 = zeros(self.block_size/2, dtype = int16)
			for i in range(self.block_size/2):
				state1[self.block_size/2-1-self.pbox_right_DC[self.block_size/2-1-i]] = state[i]
		return state1
	def get_state_through_per_left(self,model_goal,state):
		if model_goal == "LBAS":
			state1 = np.arange(1,self.nibble+1)
			for i in range(self.nibble):
				state1[self.nibble-1-self.pbox_left_LBAS[self.nibble-1-i]] = state[i]
		elif model_goal == "AS" or model_goal == "DC":
			state1 = np.arange(1,self.block_size/2+1)
			for i in range(self.block_size/2):
				state1[self.block_size/2-1-self.pbox_left_DC[self.block_size/2-1-i]] = state[i]
		return state1
		########################### xor ###############################	
	def get_state_through_xor(self,model_goal,r):
		if model_goal == "LBAS":
			state1 = np.arange((r+1)*self.block_size/8+1,(r+2)*self.block_size/8+1)
		elif model_goal == "AS" or model_goal == "DC":
			state1 = np.arange(r*self.block_size+self.block_size/2+1,r*self.block_size+self.block_size/2+1+self.block_size/2)
		return state1
	def diff_propagation_of_xor(self,model_goal,model_filename,r,state1,state2,state3):
		var = {}
		if model_goal == "LBAS":
			for i in range(self.nibble):
				var["x"] = state1[i]
				var["y"] = state2[i]
				var["z"] = state3[i]
				var["d"] = self.nibble*(r-1)+i+1
				function.diff_propagation_of_xor_word(model_filename,var)				
		elif model_goal == "AS" or model_goal == "DC":
			for i in range(self.block_size/2):
				var["x"] = state1[i]
				var["y"] = state2[i]
				var["z"] = state3[i]
				var["d"] = self.block_size/2*(r-1)+i+1
				function.diff_propagation_of_xor_bit(model_filename,var)	

	def lowbound_of_sbox(self,model_goal,model_filename,r1,r2,lb):
		if model_goal == "LBAS" or model_goal == "AS":
			if model_goal == "LBAS":
				r1 = r1+1
				r2 = r2+1
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

	def fix_diff(self,model_goal,model_filename,diff_type,r,diff_pattern):
		diff_pattern_copy = copy.deepcopy(diff_pattern)
		if diff_type == "diff_pattern_front":
			r1 = r-1
			diff_pattern_copy[0:2] = copy.deepcopy(diff_pattern[6:8])
			diff_pattern_copy[2:8] = copy.deepcopy(diff_pattern[0:6])
		elif diff_type == "diff_pattern_next":
			r1 = r+1
		if model_goal == "LBAS" or model_goal == "AS":
			if model_goal == "LBAS":
				r = r+1
				r1 = r1+1
			for i in range(1,self.nibble+1):
				with open(model_filename, "a") as f:
					f.write("A%d = 0\n"%(self.nibble*(r-1)+i))
			for i in range(1,self.nibble+1):
				if diff_pattern_copy[i-1] == 1:
					with open(model_filename, "a") as f:
						f.write("A%d = 1\n"%(self.nibble*(r1-1)+i))
				else:
					with open(model_filename, "a") as f:
						f.write("A%d = 0\n"%(self.nibble*(r1-1)+i))
		elif model_goal == "DC":
			for i in range(1,self.nibble+1):
				with open(model_filename, "a") as f:
					f.write("p%d = 0\n"%(self.nibble*(r-1)*2+i*2))
			for i in range(1,self.nibble+1):
				if diff_pattern_copy[i-1] == 1:
					with open(model_filename, "a") as f:
						f.write("p%d = 1\n"%(self.nibble*(r1-1)*2+i*2))
				else:
					with open(model_filename, "a") as f:
						f.write("p%d = 0\n"%(self.nibble*(r1-1)*2+i*2))