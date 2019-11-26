from gurobipy import *
import time
import copy
import function
import math
import re

class class_model:
	def __init__(self,cipher,model_param):
		self.cipher = cipher
		self.model_goal = model_param["model_goal"]
		self.model_round = model_param["model_round"]
		self.const_sbox = model_param["const_sbox"]
		self.const_diff = model_param["const_diff"]
		self.obj_compare = model_param["obj_compare"]
		self.model_filename = "model/"+self.cipher.name+"/"+self.model_goal+"/"+str(self.model_round)+"_round_"+str(self.const_sbox)+str(self.const_diff)+"_model.lp"
		
		self.build_model = self.build_model()
		self.model_obj = self.solve_model()

	def build_model(self):
		self.obj_fun = self.obj_fun()
		self.constraint = self.constraint()
		self.binary = self.binary()

	def obj_fun(self): 
		with open(self.model_filename, "w") as f:
			f.write("Minimize\n")
		self.cipher.obj_fun(self.model_goal,self.model_filename,1,self.model_round)
		with open(self.model_filename, "a") as f:
			f.write("\n")
	
	def constraint(self):
		with open(self.model_filename, "a") as f:
			f.write("Subject To\n")
		if self.const_sbox != []:
			if self.const_sbox == "get_upperbound_1":
				self.constraint_upperbound(1)
				self.const_sbox = []
			elif self.const_sbox == "get_upperbound_2":
				self.constraint_upperbound(2)
				self.const_sbox = []
			else:
				for i in range(len(self.const_sbox)):
					r1 = self.const_sbox[i][0]
					r2 = self.const_sbox[i][1]
					Na = self.const_sbox[i][2]
					if self.cipher.structure == "sp" and Na >= 2:
						self.cipher.lowbound_of_sbox(self.model_goal,self.model_filename,r1,r2,Na)
					elif self.cipher.structure == "feistel" and Na >= 1:
						self.cipher.lowbound_of_sbox(self.model_goal,self.model_filename,r1,r2,Na)
		if self.const_diff != []:
			diff_type = self.const_diff[0]
			diff_round = self.const_diff[1]
			diff_value = self.const_diff[2]
			self.cipher.fix_diff(self.model_goal,self.model_filename,diff_type,diff_round,diff_value)
			
		if self.cipher.structure == "sp":
			self.constraints_state_sp()
		elif self.cipher.structure == "feistel":
			self.constraints_state_feistel()
		self.input_non_zero() #at least one sbox must be active
		
	def constraints_state_sp(self):
		state = self.cipher.gen_input_state()
		for r in range (1,self.model_round+1):
			state_through_sbox = self.cipher.get_state_through_sbox(r)
			self.cipher.diff_propagation_of_sbox(self.model_goal,self.model_filename,r,state,state_through_sbox)
			state_through_per = self.cipher.get_state_through_per(state_through_sbox)
			state = copy.deepcopy(state_through_per)
	
	def constraints_state_feistel(self):
		state = self.cipher.gen_input_state(self.model_goal)
		if self.model_goal == "LBAS" and self.cipher.oriented == "byte":
			state_r = state[0:len(state)/2]
			state_l = state[len(state)/2:len(state)]
		else:
			state_l = state[0:len(state)/2]
			state_r = state[len(state)/2:len(state)]
		for r in range (1,self.model_round+1):
			if self.model_goal == "LBAS" and self.cipher.oriented == "byte":
				sc_state_l = copy.deepcopy(state_l) 
			else:
				sc_state_l = self.cipher.get_state_through_sbox(self.model_goal,r)
			if self.cipher.name == "lblock":
				ls_state_r = self.cipher.get_state_through_per_right(self.model_goal,state_r)
				p_state_l = self.cipher.get_state_through_per_left(self.model_goal,sc_state_l)
				self.cipher.diff_propagation_of_sbox(self.model_goal,self.model_filename,r,state_l,sc_state_l)
				if r < self.model_round:
					state_r = state_l
					xor_state = self.cipher.get_state_through_xor(self.model_goal,r)
					self.cipher.diff_propagation_of_xor(self.model_goal,self.model_filename,r,p_state_l,ls_state_r,xor_state)
					state_l = xor_state
			elif self.cipher.name == "twine":
				self.cipher.diff_propagation_of_sbox(self.model_goal,self.model_filename,r,state_l,sc_state_l)
				if r < self.model_round:
					xor_state = self.cipher.get_state_through_xor(self.model_goal,r)
					self.cipher.diff_propagation_of_xor(self.model_goal,self.model_filename,r,sc_state_l,state_r,xor_state)
					state_r =  self.cipher.get_state_through_per_right(self.model_goal,state_l)
					state_l =  self.cipher.get_state_through_per_left(self.model_goal,xor_state)	
					
	
	def constraint_upperbound(self,order):
		filename = "result/"+self.cipher.name+"/"+self.model_goal+"/" + str(self.model_round-1) +"_round_[][]_optimal_solution.txt"
		if order == 1:
			add_num = 0
		elif order == 2:
			add_num = self.cipher.nibble
		if self.model_goal == "AS":
			fr = open (filename,"r")
			for v in fr:
				if v[0] == "A":
					var_index,var_value = int(re.findall(r'(-?[\d]+)',v)[0]),int(re.findall(r'(-?[\d]+)',v)[1])
					with open(self.model_filename, "a") as f:
						f.write("A%d = %d\n"%(var_index+add_num,var_value))
		elif self.model_goal == "DC":
			if self.cipher.name == "present" or self.cipher.name == "rectangle" or self.cipher.name == "lblock" or self.cipher.name == "twine":
				fr = open (filename,"r")
				for v in fr:
					if v[0] == "p": #<= (self.model_round-1)*self.cipher.nibble*self.cipher.num_of_p_var:
						var_index,var_value = int(re.findall(r'(-?[\d]+)',v)[0]),int(re.findall(r'(-?[\d]+)',v)[1])
						if var_index%2 == 0:
							with open(self.model_filename, "a") as f:
								f.write("p%d = %d\n"%(var_index+add_num*2,var_value))
			elif self.cipher.name == "gift":
				var = [0 for i in range(self.cipher.nibble*3*(self.model_round-1))]
				fr = open (filename,"r")
				for v in fr:
					if v[0] == "p": #<= (self.model_round-1)*self.cipher.nibble*self.cipher.num_of_p_var:
						var_index,var_value = int(re.findall(r'(-?[\d]+)',v)[0]),int(re.findall(r'(-?[\d]+)',v)[1])
						var[var_index-1] = var_value
				for i in range(self.cipher.nibble*3*(self.model_round-1)):
					if i%3 == 0:
						with open(self.model_filename, "a") as f:
							f.write("p%d + p%d + p%d = %d\n"%(i+1+add_num*3,i+2+add_num*3,i+3+add_num*3,max(var[i:i+3])))

	def input_non_zero(self):
		if self.model_goal == "LBAS" and self.cipher.oriented == "byte":
			for i in range (1,self.cipher.nibble*2):
				with open(self.model_filename, "a") as f:
					f.write("A%d + "%(i))
			with open(self.model_filename, "a") as f:
				f.write("A%d >= 1\n"%(self.cipher.nibble*2))
		else:
			for i in range (1,self.cipher.block_size):
				with open(self.model_filename, "a") as f:
					f.write("x%d + "%(i))
			with open(self.model_filename, "a") as f:
				f.write("x%d >= 1\n"%(self.cipher.block_size))
	def binary(self): 
		with open(self.model_filename, "a") as f:
			f.write("Binary\n")
		if self.model_goal == "LBAS":
			var_dict = self.cipher.var_and_num_LBAS.copy()
			var_key = self.cipher.var_and_num_LBAS.keys()
		elif self.model_goal == "AS":	
			var_dict = self.cipher.var_and_num_AS.copy()
			var_key = self.cipher.var_and_num_AS.keys()
		elif self.model_goal == "DC":	
			var_dict = self.cipher.var_and_num_DC.copy()
			var_key = self.cipher.var_and_num_DC.keys()
		for i in var_key:
			for j in range(1,self.model_round*var_dict[i][0]+var_dict[i][1]+1):
				with open(self.model_filename, "a") as f:
					f.write(i+"%d\n"%(j))
		with open(self.model_filename, "a") as f:
			f.write("End")
	
	def solve_model(self):
		time_start = time.time()
		m = read(self.model_filename)
		m.Params.MIPFocus = 2 # MIPFocus
		m.optimize()
		time_end = time.time()
		timespend = time_end - time_start
		if m.Status == 2:
			temp = int(round(m.objVal))
		elif m.Status == 3:
			temp = 25600
		if temp < self.obj_compare:
			if self.cipher.structure == "sp":
				optimal_solution_file = "txt/"+self.cipher.name+"/"+self.model_goal+"/optimal_solution_of_submodel/"+str(self.model_round)+"_round_"+str(self.const_sbox)+str(self.const_diff)+"_optimal_solution.txt"									
			elif self.cipher.structure == "feistel":
				optimal_solution_file = "txt/"+self.cipher.name+"/"+self.model_goal+"/optimal_solution_of_submodel/"+str(self.model_round)+"_round_[][]_optimal_solution.txt"									
			f = open (optimal_solution_file,"w")
			f.write("solving the model " + str(self.model_filename) +"\n")
			f.write('obj is: %d\n'%(m.objVal))
			f.write('time is %d s.\n'%(timespend))
			for v in m.getVars():
				f.write("%s = %d\n"%(v.varName,int(round(v.x))))
			f.close()
		function.remove_file(self.model_filename)
		return temp