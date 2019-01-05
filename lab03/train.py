import sys
import collections
import re

gene_character_list = ['A', 'C', 'G', 'T']
transition_probability = {'self' : 0.999, 'other' : 0.001, 'progress' : 1}

# get file name
if (len(sys.argv) < 3 ) :
  print("gene file name, intergene file name needed")
  exit(-1)

gene_file_name = sys.argv[1]
inter_gene_file_name = sys.argv[2]

# get gene data from training file
with open(gene_file_name) as f :
  gene_data_list = [i for i in f.read().split()][1::2]

with open(inter_gene_file_name) as f :
  inter_gene_data_list = [i for i in f.read().split()][1::2]

# substituation not gene sequence character to gene sequence character
def substitute_non_gene_character(data_list) :
  for i in range(len(data_list)) :
    for j in range(len(data_list[i])) :
      if data_list[i][j] not in gene_character_list :
        data_list[i] = data_list[i][:j] +'C' + data_list[i][j + 1:]

substitute_non_gene_character(gene_data_list)
substitute_non_gene_character(inter_gene_data_list)

def get_gene_counter(sequence) :
  emission_probability = {"A" : 1, "C" : 1,  "G" : 1, "T" : 1}
  start_emission_probability = [{"A" : 1, "C" : 1,  "G" : 1, "T" : 1}, {"A" : 1, "C" : 1,  "G" : 1, "T" : 1}, {"A" : 1, "C" : 1,  "G" : 1, "T" : 1}]
  end_emission_probability = [{"A" : 1, "C" : 1,  "G" : 1, "T" : 1}, {"A" : 1, "C" : 1,  "G" : 1, "T" : 1}, {"A" : 1, "C" : 1,  "G" : 1, "T" : 1}]

  for i in sequence :
    for j in range(3) :
      start_emission_probability[j][i[j]] = start_emission_probability[j][i[j]] + 1

    emission_probability['A'] = i[3:-3].count('A') + emission_probability['A']
    emission_probability['C'] = i[3:-3].count('C') + emission_probability['C']
    emission_probability['G'] = i[3:-3].count('G') + emission_probability['G']
    emission_probability['T'] = i[3:-3].count('T') + emission_probability['T']

    for j in range(3) :
      end_emission_probability[j][i[-3 + j]] = end_emission_probability[j][i[-3 + j]] + 1
  return emission_probability, start_emission_probability, end_emission_probability

def get_counter(sequence) :
  emission_probability = {"A" : 1, "C" : 1,  "G" : 1, "T" : 1}
  for i in sequence :
    emission_probability['A'] = i.count('A') + emission_probability['A']
    emission_probability['C'] = i.count('C') + emission_probability['C']
    emission_probability['G'] = i.count('G') + emission_probability['G']
    emission_probability['T'] = i.count('T') + emission_probability['T']
  return emission_probability

gene_region, start_region, end_region = get_gene_counter(gene_data_list)

total_gene_counter = sum(list(map(int, list(gene_region.values()))))
gene_region_probability = [str(i / total_gene_counter) for i in gene_region.values()]
start_region_probability = [[str(j / (len(gene_data_list) + 4)) for j in start_region[i].values()] for i in range(3)]
end_region_probability = [[str(j / (len(gene_data_list) + 4)) for j in end_region[i].values()] for i in range(3)]

inter_gene_region = get_counter(inter_gene_data_list)
total_gene_counter = sum(list(map(int, list(inter_gene_region.values()))))
inter_gene_region_probability = [str(i / total_gene_counter) for i in inter_gene_region.values()]

with open("model.txt", "w") as f :
  f.write("emission probability\n")
  for i in range(3) :
    f.write("start_{}\n{}\n".format(str(i + 1), " ".join(start_region_probability[i])))
  for i in range(3) :
    f.write("end_{}\n{}\n".format(str(i + 1), " ".join(end_region_probability[i])))
  f.write("genic\n{}\n".format(" ".join(gene_region_probability)))
  f.write("intergenic\n{}".format(" ".join(inter_gene_region_probability)))
