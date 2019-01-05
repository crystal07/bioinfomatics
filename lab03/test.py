import sys
import collections
import math
import re
import datetime

if len(sys.argv) < 2 :
  print("gene file name needed")
  exit(-1)

gene_file_name = sys.argv[1]
if len(sys.argv) >= 3 :
	model_file_name = sys.argv[2]
else :
	model_file_name = "model.txt"

with open(gene_file_name) as f :
	sequence = f.read().split()[1:]
	sequence = "".join(sequence)


gene_character_list = ["A", "C", "G", "T"]
states = ['start_1', 'start_2', 'start_3', 'end_1', 'end_2', 'end_3', 'genic', 'intergenic']
# start_p = {'start_1': 0.5, 'start_2': 0.0, 'start_3': 0.0, 'end_1' : 0.0, 'end_2' : 0.0, 'end_3' : 0.0, 'intergenic' : 0.5, 'genic' : 0.0}
# start_p = {'start_1': 0.5, 'start_2': sys.float_info.min * sys.float_info.epsilon, 'start_3': sys.float_info.min * sys.float_info.epsilon, 'end_1' : sys.float_info.min * sys.float_info.epsilon, 'end_2' : sys.float_info.min * sys.float_info.epsilon, 'end_3' : sys.float_info.min * sys.float_info.epsilon, 'intergenic' : 0.5, 'genic' : sys.float_info.min * sys.float_info.epsilon}
# start_p = {'start_1' : 0.5, 'intergenic' : 0.5}
start_p = {0 : 0.5, 6 : 0.5}
states = [0, 1, 2, 3, 4, 5, 6, 7]

state_mapping = {
	'start_1' : 0,
	'start_2' : 1,
	'start_3' : 2,
	'end_1' : 3,
	'end_2' : 4,
	'end_3' : 5,
	'intergenic' : 6,
	'genic' : 7
}

trans_p = {
	0 : {1 : 1},
	1 : {2 : 1},
	2 : {7 : 1},
	3 : {4 : 1},
	4 : {5 : 1},
	5 : {6 : 1},
	6 : {0 : 0.001, 6 : 0.999},
	7 : {7 : 0.999, 3 : 0.001},
}

emit_p = {}
with open(model_file_name) as f :
	f.readline()
	probabilities = f.read().split('\n')
	for i in range(len(probabilities))[::2] :
		a_prob, c_prob, g_prob, t_prob = map(float, probabilities[i+1].split())
		emit_p[state_mapping[probabilities[i]]] = {'A' : a_prob, 'C' : c_prob, 'G' : g_prob, 'T' : t_prob}

def viterbi(gene, states, start_p, trans_p, emit_p) :
	print(datetime.datetime.now())
	prob_list = [[]]
	prev_list = [[]]
	for i in states : 
		if i in start_p :
			prob_list[0].append(math.log(start_p[i]) + math.log(emit_p[i][gene[0]]))

		else :
			prob_list[0].append(float('-inf'))
		prev_list[0].append(None)

	for i in range(len(gene)) :
		prob_list.append([])
		prev_list.append([])
		for state in states :
			
			max_list = [prob_list[i][prev_state] + math.log(trans_p[prev_state][state]) + math.log(emit_p[state][gene[i]]) if state in trans_p[prev_state] else float('-inf') for prev_state in states]
			prob_list[i + 1].append(max(max_list))
			prev_list[i + 1].append(max_list.index(prob_list[i + 1][-1]))
			# if state >= 6 :
			# 	print(state, max_list.index(prob_list[i + 1][-1]))

			# print(True if prob_list[i] == prob_list[i + 1] else False)
			# print(gene[i], state, prob_list[i + 1][-1] if prev_list[i] == prev_list[i + 1] else False)

	print(datetime.datetime.now())

	with open("output_table.txt", "w") as f :
		for i in range(len(prob_list)) :
			f.write("({})\n".format(" ".join(map(str, zip(prob_list[i], prev_list[i])))))

	max_prev = prev_list[-1][prob_list[-1].index(max(prob_list[-1]))]
	x = [str(-1) for i in range(len(prob_list))]
	x[-1] = str(max_prev)
	for i in reversed(range(len(prob_list))) :
		max_prev = prev_list[i][max_prev]
		# print("case", i)
		# print(prob_list[i], max(prob_list[i]))
		# print(prev_list[i], max_prev)
		# x.insert(0, max_prev)
		x[i] = str(max_prev)

	print(datetime.datetime.now())

	with open("output_test.txt", "w") as f :
		f.write(" ".join(x[1:]))
	
	return x

result = viterbi(sequence, states, start_p, trans_p, emit_p)
result_pair = re.finditer('01?2?[7]*3?4?5?', "".join(result[1:]))
with open("result.txt", "w") as f :
	index = 0
	for m in result_pair :
		f.write("gene {} : {} {}\n".format(index, m.start(), m.end()))
		index = index + 1
