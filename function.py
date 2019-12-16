import os
import re
from numpy import *
import numpy as np 
import math
import copy
import shutil

def gen_filefolder(cipher_name,goal):
	filefolder = "model/"+cipher_name+"/LBAS/"
	if not os.path.isdir(filefolder):
		os.makedirs(filefolder)
	filefolder = "model/"+cipher_name+"/"+str(goal)+"/"
	if not os.path.isdir(filefolder):
		os.makedirs(filefolder)
	filefolder = "result/"+cipher_name+"/"+str(goal)+"/"
	if not os.path.isdir(filefolder):
		os.makedirs(filefolder)
	filefolder = "txt/"+cipher_name+"/"+goal+"/optimal_solution_of_submodel/"
	if not os.path.isdir(filefolder):
		os.makedirs(filefolder)

def remove_file(filename):
	if os.path.exists(filename):
		os.remove(filename)

def diff_propagation_of_sbox_LBAS_bit(sbox_size,branch_num_of_sbox,model_filename,var):
	for i in range(sbox_size-1):
		with open(model_filename, "a") as f:
			f.write("x%d + "%(var["x"][i]))
	# constraint 1
	with open(model_filename, "a") as f:
		f.write("x%d - A%d >= 0\n"%(var["x"][sbox_size-1],var["A"]))
	for i in range(sbox_size):
		with open(model_filename, "a") as f:
			f.write("A%d - x%d >= 0\n"%(var["A"],var["x"][i]))
	# constraint 2
	for i in range(sbox_size):
		with open(model_filename, "a") as f:
			f.write("x%d + "%(var["x"][i]))
	for i in range(sbox_size-1):
		with open(model_filename, "a") as f:			
			f.write("x%d + "%(var["y"][i]))
	with open(model_filename, "a") as f:	
		f.write("x%d - %d d%d >= 0\n"%(var["y"][sbox_size-1],branch_num_of_sbox,var["d"]))
	for i in range(sbox_size):
		with open(model_filename, "a") as f:	
			f.write("d%d - x%d >= 0\n"%(var["d"],var["x"][i]))
	for i in range(sbox_size):
		with open(model_filename, "a") as f:	
			f.write("d%d - x%d >= 0\n"%(var["d"],var["y"][i]))
	# constraint 3
	with open(model_filename, "a") as f:	
		f.write("%d x%d + %d x%d + %d x%d + %d x%d - x%d - x%d - x%d - x%d >= 0\n"%(sbox_size,var["y"][0],sbox_size,var["y"][1],sbox_size,var["y"][2],sbox_size,var["y"][3],var["x"][0],var["x"][1],var["x"][2],var["x"][3]))
		f.write("%d x%d + %d x%d + %d x%d + %d x%d - x%d - x%d - x%d - x%d >= 0\n"%(sbox_size,var["x"][0],sbox_size,var["x"][1],sbox_size,var["x"][2],sbox_size,var["x"][3],var["y"][0],var["y"][1],var["y"][2],var["y"][3]))

def diff_propagation_of_sbox_AS(sbox_size,model_filename,var,ine):
	vector_size = 2 * sbox_size 
	symbal = [" " for i in range(vector_size+1)]
	for i in range(sbox_size-1):
		with open(model_filename, "a") as f:
			f.write("x%d + "%(var["x"][i]))
	with open(model_filename, "a") as f:
		f.write("x%d - A%d >= 0\n"%(var["x"][sbox_size-1],var["A"]))
	for i in range(sbox_size):
		with open(model_filename, "a") as f:
			f.write("A%d - x%d >= 0\n"%(var["A"],var["x"][i]))
	for i in range(len(ine)):
		for k in range (0,vector_size):
			if ine[i,k] < 0:
				symbal[k] = '-'
			elif ine[i,k] >= 0:
				symbal[k] = '+'
		for k1 in range(sbox_size):
			with open(model_filename, "a") as f:
				f.write("%s %d x%d "%(symbal[k1], abs(ine[i,k1]),var["x"][k1]))
		for k2 in range (sbox_size,sbox_size*2):
			with open(model_filename, "a") as f:
				f.write("%s %d x%d "%(symbal[k2], abs(ine[i,k2]),var["y"][k2-sbox_size]))
		if ine[i,vector_size] > 0:
			symbal[vector_size] = '-'
		elif ine[i,vector_size] <= 0:
			symbal[vector_size] = '+'
		with open(model_filename, "a") as f:
			f.write(" >= %s %d \n"%(symbal[vector_size], abs(ine[i,vector_size])))

def diff_propagation_of_sbox_DC(sbox_size,num_of_p_var,model_filename,var,ine):
	vector_size = 2 * sbox_size + num_of_p_var 
	symbal = [" " for i in range(vector_size+1)]
	for i in range(len(ine)):
		for k in range (0,vector_size):
			if ine[i,k] < 0:
				symbal[k] = '-'
			elif ine[i,k] >= 0:
				symbal[k] = '+'
		for k1 in range(sbox_size):
			with open(model_filename, "a") as f:
				f.write("%s %d x%d "%(symbal[k1], abs(ine[i,k1]),var["x"][k1]))
		for k2 in range (sbox_size,sbox_size*2):
			with open(model_filename, "a") as f:
				f.write("%s %d x%d "%(symbal[k2], abs(ine[i,k2]),var["y"][k2-sbox_size]))
		for k3 in range(sbox_size*2,vector_size):
			with open(model_filename, "a") as f:
				f.write("%s %d p%d "%(symbal[k3], abs(ine[i,k3]),var["p"][k3-sbox_size*2]))
		if ine[i,vector_size] > 0:
			symbal[vector_size] = '-'
		elif ine[i,vector_size] <= 0:
			symbal[vector_size] = '+'
		with open(model_filename, "a") as f:
			f.write(" >= %s %d \n"%(symbal[vector_size], abs(ine[i,vector_size])))

def diff_propagation_of_xor_word(model_filename,var):
	with open(model_filename, "a") as f:
		f.write("A%d + A%d + A%d - 2 d%d >= 0\n"%(var["x"], var["y"], var["z"], var["d"]))
		f.write("d%d - A%d >= 0\n"%(var["d"],var["x"]))
		f.write("d%d - A%d >= 0\n"%(var["d"],var["y"]))
		f.write("d%d - A%d >= 0\n"%(var["d"],var["z"]))
	
def diff_propagation_of_xor_bit(model_filename,var):
	with open(model_filename, "a") as f:
		f.write("x%d + x%d + x%d - 2 d%d >= 0\n"%(var["x"], var["y"], var["z"], var["d"]))
		f.write("d%d - x%d >= 0\n"%(var["d"],var["x"]))
		f.write("d%d - x%d >= 0\n"%(var["d"],var["y"]))
		f.write("d%d - x%d >= 0\n"%(var["d"],var["z"]))
		f.write("x%d + x%d + x%d <= 2\n"%(var["x"], var["y"], var["z"]))


#################################################################################
# searching function
def get_order_of_Na(cipher,r):
	if r < cipher.round_of_first_region2:
		temp = [1,2]
	elif r >= cipher.round_of_first_region2:
		temp = [2,1]
	return temp
	
def get_value(r,value):
	if r == 0:
		return 0
	elif r >= 1:
		return value

def get_search_round(r):
	middle = int((r+2)/2)
	array_1 = middle - np.arange(1,middle)
	array_2 = np.arange(middle+1,r+1)
	temp = [middle]
	if r%2 == 0:
		for i in range(max(len(array_1),len(array_2))):
			if i < len(array_1):
				temp.append(array_1[i])
			if i < len(array_2):
				temp.append(array_2[i])
	elif r%2 == 1:
		for i in range(max(len(array_1),len(array_2))):
			if i < len(array_2):
				temp.append(array_2[i])
			if i < len(array_1):
				temp.append(array_1[i])
	return array(temp)
	
def trans_diff_pattern(cipher,Na):
	diff_pattern = []
	position = list(combinations(np.arange(1,cipher.nibble+1),Na))
	for i in range(len(position)):
		pattern = [0 for i in range(cipher.nibble)]
		for j in positoin[i]:
			pattern[j-1] = 1 
		diff_pattern.append(copy.deepcopy(pattern))
	return diff_pattern

def read_txt(filename):
	line = []
	with open(filename) as f:
		for j in f.readlines():
			line.append(j)
	return line		

def get_var_from_two_submodels(cipher,goal,r,Na,i,diff):
	filename_1 = "txt/"+cipher.name+"/"+goal+"/optimal_solution_of_submodel/"+str(i-1)+"_round_"+str([[1,i-1,Na]])+str(["output_diff",i-1,diff])+"_optimal_solution.txt"									
	filename_2 = "txt/"+cipher.name+"/"+goal+"/optimal_solution_of_submodel/"+str(r-i+1)+"_round_"+str([[2,r-i+1,Na]])+str(["input_diff",1,diff])+"_optimal_solution.txt"									
	filename = "txt/"+cipher.name+"/"+goal+"/optimal_solution_of_submodel/"+"/"+str(r)+"_round_[][]_optimal_solution.txt"
	with open(filename,"w") as f:
		f.write("get variables from "+str(i-1)+"_round_"+str([[1,i-1,Na]])+str(diff)+"_and_"+str(r-i+1)+"_round_"+str([[2,r-i+1,Na]])+str(diff)+"_optimal_solutions\n")
	if goal == "AS":
		x_var = [0 for j in range (r*cipher.var_and_num_AS["x"][0]+cipher.var_and_num_AS["x"][1])]
		x_var_2 = [0 for j in range ((r-i+1)*cipher.var_and_num_AS["x"][0]+cipher.var_and_num_AS["x"][1])]
		A_var = [0 for j in range (r*cipher.var_and_num_AS["A"][0]+cipher.var_and_num_AS["A"][1])]
	elif goal == "DC":
		x_var = [0 for j in range (r*cipher.var_and_num_DC["x"][0]+cipher.var_and_num_DC["x"][1])]
		x_var_2 = [0 for j in range ((r-i+1)*cipher.var_and_num_DC["x"][0]+cipher.var_and_num_DC["x"][1])]
		p_var = [0 for j in range (r*cipher.var_and_num_DC["p"][0]+cipher.var_and_num_DC["p"][1])]
	if i > 1:
		var = read_txt(filename_1)
		for j in var:
			if j[0] == "x":
				x_index = int(re.findall(r'(-?[\d]+)',j)[0])
				x_value = int(re.findall(r'(-?[\d]+)',j)[1])
				x_var[x_index-1] = x_value
			elif j[0] == "A" and goal == "AS":
				A_index = int(re.findall(r'(-?[\d]+)',j)[0])
				A_value = int(re.findall(r'(-?[\d]+)',j)[1])
				A_var[A_index-1] = A_value
			elif j[0] == "p"and goal == "DC":
				p_index = int(re.findall(r'(-?[\d]+)',j)[0])
				p_value = int(re.findall(r'(-?[\d]+)',j)[1])
				p_var[p_index-1] = p_value
	var = read_txt(filename_2)
	for j in var:
		if j[0] == "x":
			x_index = int(re.findall(r'(-?[\d]+)',j)[0])
			x_value = int(re.findall(r'(-?[\d]+)',j)[1])
			x_var_2[x_index-1] = x_value
		if j[0] == "A" and goal == "AS":
			A_index = int(re.findall(r'(-?[\d]+)',j)[0])
			A_value = int(re.findall(r'(-?[\d]+)',j)[1])
			A_var[A_index+cipher.nibble*(i-1)-1] = A_value
		if j[0] == "p" and goal == "DC":
			p_index = int(re.findall(r'(-?[\d]+)',j)[0])
			p_value = int(re.findall(r'(-?[\d]+)',j)[1])
			p_var[p_index+(i-1)*cipher.var_and_num_DC["p"][0]-1] = p_value
	state1 = cipher.gen_input_state()
	state2 = cipher.gen_input_state()	
	for j in range (1,r+2):
		if j >= i: 
			for k in range(len(state2)):
				x_var[state1[k]-1] = x_var_2[state2[k]-1]
			Sstate2 = cipher.get_state_through_sbox(j-i+1)
			Pstate2 = cipher.get_state_through_per(Sstate2)
			state2 = Pstate2		
		Sstate1 = cipher.get_state_through_sbox(j)
		Pstate1 = cipher.get_state_through_per(Sstate1)
		state1 = Pstate1
	with open(filename,"a") as f:
		for j in range(len(x_var)):
			f.write("x%d = %d\n"%(j+1,x_var[j]))
		if goal =="AS":
			for j in range(len(A_var)):
				f.write("A%d = %d\n"%(j+1,A_var[j]))	
		elif goal =="DC":
			for j in range(len(p_var)):
				f.write("p%d = %d\n"%(j+1,p_var[j]))		

def get_trail_sp(cipher,goal,r):
	if goal == "AS":
		x_var = [0 for j in range (r*cipher.var_and_num_AS["x"][0]+cipher.var_and_num_AS["x"][1])]
	elif goal == "DC":
		x_var = [0 for j in range (r*cipher.var_and_num_DC["x"][0]+cipher.var_and_num_DC["x"][1])]	
	shutil.copy("txt/"+cipher.name+"/"+goal+"/optimal_solution_of_submodel/"+str(r)+"_round_[][]_optimal_solution.txt", "result/"+cipher.name+"/"+goal+"/")
	
	var = read_txt("txt/"+cipher.name+"/"+goal+"/optimal_solution_of_submodel/"+str(r)+"_round_[][]_optimal_solution.txt")
	for i in var:
		if i[0] == "x":
			x_index = int(re.findall(r'(-?[\d]+)',i)[0])
			x_value = int(re.findall(r'(-?[\d]+)',i)[1])
			x_var[x_index-1] = x_value
	file_path = "result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_best_path.txt"
	with open(file_path,"w") as f:
		f.write("x_var list is " + str(x_var)+"\n")
	state = cipher.gen_input_state()
	for r in range (1,r+1):
		with open(file_path,"a") as f:
			f.write("\n%dth round input diff of Sbox: \n"%(r))
		for i in range(len(state)):
			if i%4 == 0:
				with open(file_path,"a") as f:
					f.write(" ")
			with open(file_path,"a") as f:
				f.write(str(x_var[state[i]-1]))
		state_through_sbox = cipher.get_state_through_sbox(r)		
		# print "state_through_sbox",state_through_sbox
		with open(file_path,"a") as f:
			f.write("\n%dth round output diff of Sbox: \n"%(r))
		for i in range(len(state)):
			if i%4 == 0:
				with open(file_path,"a") as f:
					f.write(" ")
			with open(file_path,"a") as f:
				f.write(str(x_var[state_through_sbox[i]-1]))
		state_through_per = cipher.get_state_through_per(state_through_sbox)
		state = copy.deepcopy(state_through_per)
		# print "state_through_per",state
		

def get_trail_feistel(cipher,goal,r):
	if goal == "AS":
		x_var = [0 for j in range (r*cipher.var_and_num_AS["x"][0]+cipher.var_and_num_AS["x"][1])]
	elif goal == "DC":
		x_var = [0 for j in range (r*cipher.var_and_num_DC["x"][0]+cipher.var_and_num_DC["x"][1])]	
	shutil.copy("txt/"+cipher.name+"/"+goal+"/optimal_solution_of_submodel/"+"/"+str(r)+"_round_[][]_optimal_solution.txt", "result/"+cipher.name+"/"+goal+"/")
	var = read_txt("txt/"+cipher.name+"/"+goal+"/optimal_solution_of_submodel/"+"/"+str(r)+"_round_[][]_optimal_solution.txt")
	
	for i in var:
		if i[0] == "x":
			x_index = int(re.findall(r'(-?[\d]+)',i)[0])
			x_value = int(re.findall(r'(-?[\d]+)',i)[1])
			x_var[x_index-1] = x_value
	file_path = "result/"+cipher.name+"/"+goal+"/"+str(r)+"_round_best_path.txt"
	with open(file_path,"w") as f:
		f.write("x_var list is " + str(x_var)+"\n")
	state = cipher.gen_input_state(goal)
	state_l = state[0:len(state)/2]
	state_r = state[len(state)/2:len(state)]
	for r in range(1,r+1):
		with open(file_path,"a") as f:
			f.write("\n%dth round left: input diff of Sbox: \n"%(r))
		for i in range(len(state_l)):
			if i%4 == 0:
				with open(file_path,"a") as f:
					f.write(" ")
			with open(file_path,"a") as f:
				f.write(str(x_var[state_l[i]-1]))
		sc_state_l = cipher.get_state_through_sbox(goal,r)
		with open(file_path,"a") as f:
			f.write("\n%dth round left: output diff of Sbox: \n"%(r))
		for i in range(len(sc_state_l)):
			if i%4 == 0:
				with open(file_path,"a") as f:
					f.write(" ")
			with open(file_path,"a") as f:
				f.write(str(x_var[sc_state_l[i]-1]))
		with open(file_path,"a") as f:
			f.write("\n%dth round right: \n"%(r))
		for i in range(len(state_r)):
			if i%4 == 0:
				with open(file_path,"a") as f:
					f.write(" ")
			with open(file_path,"a") as f:
				f.write(str(x_var[state_r[i]-1]))
		if cipher.name == "lblock":
			state_r = state_l
			xor_state = cipher.get_state_through_xor(goal,r)
			state_l = xor_state
		elif cipher.name == "twine":
			xor_state = cipher.get_state_through_xor(goal,r)
			state_r =  cipher.get_state_through_per_right(goal,state_l)
			state_l = cipher.get_state_through_per_left(goal,xor_state)	