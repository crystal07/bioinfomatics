import random
import operator
import sys
import datetime

if len(sys.argv) < 3 :
    print("program, input file name, output file name needed")
    exit(-1)

if len(sys.argv) > 3 :
    iteration = int(sys.argv[3])
else :
    iteration = 1

# get data from test file
test_file = sys.argv[1]
with open(test_file) as f :
    k = int(f.readline())
    t = int(f.readline())
    m = int(f.readline())
    sequences = f.read().upper().split()

# ramdom k-mers
def get_random_motif(sequence) :
    i = random.randint(0, len(sequence) - (k + 1))
    return sequence[i : i + k]

# profile calculation
def get_profile(motifs, index) :
    profile = [{"A" : 0, "G" : 0, "C" : 0, "T" : 0} for i in range(k)]
    for i in range(k) :
        index_character_list = [motifs[j][i] for j in range(len(motifs)) if j != index]
        profile[i]["A"] = (index_character_list.count("A") + 1) / (t + 4)
        profile[i]["G"] = (index_character_list.count("G") + 1) / (t + 4)
        profile[i]["C"] = (index_character_list.count("C") + 1) / (t + 4)
        profile[i]["T"] = (index_character_list.count("T") + 1) / (t + 4)
        # print(profile[i])
    return profile

# new motif from profile
def get_new_motif(profile, motif) :
    index = 0
    max_probability = 0
    for i in range(len(motif) - k) :
        new_probabitily = get_probability(profile, motif[i : i + k])
        if new_probabitily > max_probability :
            max_probability = new_probabitily
            index = i
    return motif[index : index + k]

def get_profile_motif(profile) :
    return "".join([max(i.items(), key=operator.itemgetter(1))[0] for i in profile])

def get_probability(profile, motif) :
    probability = 1
    for i in range(k) :
        probability = probability * profile[i][motif[i]]
    return probability

# calculate score
def get_score(motifs, consensus) :
    score = 0
    for i in motifs :
        for j in range(k) :
            if consensus[j] != i[j] :
                score = score + 1
    return score

# gibbs sampling
def gibbs_sapmling() :
    motifs = [get_random_motif(i) for i in sequences]
    best_motifs = motifs.copy()
    best_score = t * k + 100
    tmp_motifs = motifs.copy()
    for j in range(3000) :
        i = random.randint(0, t - 1)
        profile = get_profile(tmp_motifs, i)
        tmp_motifs[i] = get_new_motif(profile, sequences[i])
        profile_motif = get_profile_motif(profile)
        new_score = get_score(tmp_motifs, profile_motif)
        if new_score < best_score :
            # print("best motif changed in j iteration", j)
            # if best_motifs == tmp_motifs :
            #     print(best_motifs)
            #     break
            best_motifs = tmp_motifs.copy()
            best_profile_motif = profile_motif
            best_score = new_score
    return best_motifs, best_profile_motif

for i in range(iteration) :
    print(i, datetime.datetime.now())
    motifs, profile_motif = gibbs_sapmling()
    # for i in best_motifs :
    #     print(i)
    if i != 0 :
        if get_score(best_motifs, best_profile_motif) > get_score(motifs, profile_motif) :
            best_motifs = motifs
            best_profile_motif = profile_motif
    else :
        best_motifs = motifs
        best_profile_motif = profile_motif


profile = get_profile(best_motifs, -1)

profile_a = []
proflie_c = []
profile_g = []
profile_t = []
for i in profile :
    profile_a.append(i["A"])
    proflie_c.append(i["C"])
    profile_g.append(i["G"])
    profile_t.append(i["T"])

output_file = sys.argv[2]
with open(output_file, 'w') as f :
    # write motifs
    for i in best_motifs :
        f.write('{}\n'.format(i))
    # write profiles
    f.write('{}\n'.format(' '.join(map(str, profile_a))))
    f.write('{}\n'.format(' '.join(map(str, proflie_c))))
    f.write('{}\n'.format(' '.join(map(str, profile_g))))
    f.write('{}\n'.format(' '.join(map(str, profile_t))))

# print motifs
# for i in best_motifs :
#     print(i)
# for i in profile :
#     print(i)
