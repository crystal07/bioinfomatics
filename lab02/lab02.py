import sys

sequence1, sequence2, score_file = sys.argv[1 : 4]
if len(sys.argv) >= 5 :
	output_file = sys.argv[4]
else :
	output_file = "output.txt"

def parse_sequence(filename) :
	with open(filename) as f:
		sequence_name = f.readline().strip('\n')[1:]
		sequence = f.readline().strip('\n').upper()
	return {sequence_name : sequence}

sequence1 = parse_sequence(sequence1)
sequence2 = parse_sequence(sequence2)

score = {}
with open(score_file) as f :
	for i in range(3) :
		tmp = f.readline().strip('\n').split('=')
		score[tmp[0]] = int(tmp[1])
print(score)

def get_lcs_table(local, seq1, seq2, score) :
	# table initialize
	lcs_table = [[-10000000 for _ in range(len(seq2) + 1)] for __ in range(len(seq1) + 1)]
	backtracking_matrix = [["" for _ in range(len(seq2) + 1)] for __ in range(len(seq1) + 1)]
	for i in range(0, len(seq2) + 1) :
		lcs_table[0][i] = score['gap'] * i
		if local == True :
			if lcs_table[0][i] < 0 :
				lcs_table[0][i] = 0
		backtracking_matrix[0][i] = '-'
	for i in range(1, len(seq1) + 1) :
		lcs_table[i][0] = score['gap'] * i
		if local == True :
			if lcs_table[i][0] < 0 :
				lcs_table[i][0] = 0
		backtracking_matrix[i][0] = '|'

	if local == True :
		best_score = 0

	# calculate score
	for i in range(1, len(seq1) + 1) :
		for j in range(1, len(seq2) + 1) :
			# for gap
			if (lcs_table[i][j - 1] + score['gap']) > lcs_table[i][j] :
				lcs_table[i][j] = lcs_table[i][j - 1] + score['gap']
				backtracking_matrix[i][j] = '-'

			# for gap
			if (lcs_table[i - 1][j] + score['gap']) > lcs_table[i][j] :
				lcs_table[i][j] = lcs_table[i - 1][j] + score['gap']
				backtracking_matrix[i][j] = '|'

			# for match
			if seq1[i - 1] == seq2[j - 1] and (lcs_table[i - 1][j - 1] + score['match']) > lcs_table[i][j] :
				lcs_table[i][j] = lcs_table[i - 1][j - 1] + score['match']
				backtracking_matrix[i][j] = '`'

			# for mismatch
			elif seq1[i - 1] != seq2[j - 1] and (lcs_table[i - 1][j - 1] + score['mismatch']) > lcs_table[i][j] :
				lcs_table[i][j] = lcs_table[i - 1][j - 1] + score['mismatch']
				backtracking_matrix[i][j] = '`'

			# if want to get local alignment, compare with 0
			if local == True :
				if lcs_table[i][j] < 0 :
					lcs_table[i][j] = 0

				if lcs_table[i][j] > best_score :
					best_score = lcs_table[i][j]
					last_info = [i, j]

	if local == True :
		return lcs_table, backtracking_matrix, last_info
	else :
		return lcs_table, backtracking_matrix, None

def get_backtrack(local, backtracking_matrix, score, seq1, seq2, last=[0,0]) :

	if local == True :
		i = last[0]
		j = last[1]
	else :
		i = len(seq1)
		j = len(seq2)

	new_sequence = [[], []]
	match_list = [0, 0, 0] # match, mismatch, gap
	while True :
		if backtracking_matrix[i][j] == '`' :
			new_sequence[0].append(seq1[i - 1])
			new_sequence[1].append(seq2[j - 1])
			if seq1[i - 1] == seq2[j - 1] :
				match_list[0] = match_list[0] + 1
			else :
				match_list[1] = match_list[1] + 1
			i = i - 1
			j = j - 1
		elif backtracking_matrix[i][j] == '|' :
			new_sequence[0].append(seq1[i - 1])
			new_sequence[1].append('-')
			match_list[2] = match_list[2] + 1
			i = i - 1
		else :
			new_sequence[0].append('-')
			new_sequence[1].append(seq2[j - 1])
			match_list[2] = match_list[2] + 1
			j = j - 1

		if (i == 0 and j == 0) :
			break

		if local == True :
			if score[i][j] == 0 :
				break

	new_sequence[0] = "".join(reversed(new_sequence[0]))
	new_sequence[1] = "".join(reversed(new_sequence[1]))

	return new_sequence, match_list

def print_table(table) :
	for i in range(len(table)) :
		for j in range(len(table[i])) :
			print(table[i][j], end=" ")
		print()
print(sequence1)
print(sequence2)
lcs_table, backtracking_matrix, last = get_lcs_table(False, list(sequence1.values())[0], list(sequence2.values())[0], score)
# print_table(lcs_table)
new_sequence, match_list = get_backtrack(False, backtracking_matrix, lcs_table, list(sequence1.values())[0], list(sequence2.values())[0])
print("".join(new_sequence[0]), "".join(new_sequence[1]))
# print(match_list)

with open("global_" + output_file, "w") as f :
	f.write("matches: {}\n".format(match_list[0]))
	f.write("mismatchs: {}\n".format(match_list[1]))
	f.write("gabs: {}\n".format(match_list[2]))
	f.write("score: {}\n\n".format(match_list[0] * score['match'] + match_list[1] * score['mismatch'] + match_list[2] * score['gap']))
	f.write(">{}\n{}\n\n".format(list(sequence1.keys())[0], new_sequence[0]))
	f.write(">{}\n{}\n\n".format(list(sequence2.keys())[0], new_sequence[1]))

print(sequence1, sequence2)
lcs_table, backtracking_matrix, last = get_lcs_table(True, list(sequence1.values())[0], list(sequence2.values())[0], score)
# print_table(lcs_table)
new_sequence, match_list = get_backtrack(True, backtracking_matrix, lcs_table, list(sequence1.values())[0], list(sequence2.values())[0], last)
print("".join(new_sequence[0]), "".join(new_sequence[1]))
# print(match_list)

with open("local_" + output_file, "w") as f :
	f.write("matches: {}\n".format(match_list[0]))
	f.write("mismatchs: {}\n".format(match_list[1]))
	f.write("gabs: {}\n".format(match_list[2]))
	f.write("score: {}\n\n".format(match_list[0] * score['match'] + match_list[1] * score['mismatch'] + match_list[2] * score['gap']))
	f.write(">{}\n{}\n\n".format(list(sequence1.keys())[0], new_sequence[0]))
	f.write(">{}\n{}\n\n".format(list(sequence2.keys())[0], new_sequence[1]))

