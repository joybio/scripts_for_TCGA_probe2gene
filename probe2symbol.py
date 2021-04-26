#!/root/miniconda/bin/python

import re
import optparse
from optparse import OptionParser
import csv
import collections
from collections import defaultdict

parser = OptionParser("Usage: %prog -a annotation_file -p probe_matrix -o gene_Matrix")
parser.add_option("-a","--annotation",dest="ann",
		help = "platform annotation file: e.g. GPL570-55999")
parser.add_option("-p","--probe",dest="probe",
		help = "probe matrix: e.g. *series_matrix.txt")
parser.add_option("-o","--output",dest="out",
		help = "output file: e.g. gene_matrix")
(options,args) = parser.parse_args()

probe_file = open(options.probe,"r")
annotation_file = open(options.ann,"r")
out_file = open(options.out,"w")

sample_list = []
sample_probe_dict = defaultdict(list)
for i in probe_file:
	if i.startswith("#") or i.startswith("!"):
		pass
	elif i.startswith("ID"):
		i = i.strip().replace("\"","").split("\t")
		for j in range(len(i)):
			if j == 0:
				sample_list.append(i[j])
			else:
				sample_name = re.split(r"[._]",i[j])
				sample_name = sample_name[0]
				sample_list.append(sample_name)
			j += 1
	else:
		i = i.strip().split("\t")
		for j in range(len(i)):
			i[j] = i[j].replace("\"","")
			sample_probe_dict[sample_list[j]].append(i[j])
			j += 1
'''
			if sample_name[j] not in sample_probe_dict.keys():
				sample_probe_dict[sample_name[j]] = i[j].split("\n") 
				# 以不存在的\n split，将字符串转为list
				#sample_probe_dict[sample_name[j]] = i[j]
			else:
				sample_probe_dict[sample_name[j]].append(i[j])
				#sample_probe_dict[sample_name[j]] += "," +i[j]
				#后面再用， split成list
'''
probe_file.close()
#print(sample_probe_dict)
probe_symbol = defaultdict(str)
for i in annotation_file:
	if i.startswith("#") or i.startswith("!"):
		pass
	else:
		if i.startswith("ID"):
			key = re.findall("\t([g,G]ene [s,S]ymbol)\t",i)
			#print(key[0])
			i = i.strip().split("\t")
			#print(i)
			pos = i.index(key[0])
#			print(pos)
		else:
			i = i.strip().split("\t")
			#if re.search(".*\s*.*",i[pos])
			#	pass
			if len(i) > pos:
				if re.search(".+\s+.+",i[pos]) and i[pos] == "":
					pass
#				print(i[pos])
				else:
					gene_symbol = i[pos].split(" /// ")
					probe_symbol[i[0]] = gene_symbol[0]
annotation_file.close()
#print(probe_symbol)


probe_name = sample_probe_dict["ID_REF"] #list of probe name
del sample_probe_dict["ID_REF"] # remove key = "ID_REF"
gene_list = defaultdict(int) #list of gene name
sample_gene_exp = defaultdict(str) #keys = sample values = dict of gene expression 
for i in sample_probe_dict.keys():
	gene_exp_ave = defaultdict(list) #dict of expression average
	gene_count = defaultdict(int) #dict of gene_name number
	gene_exp_sum = defaultdict(int) #dict of expression summary
	exp = sample_probe_dict[i] #list of probe expressiom
	#print(exp)
	#exp = sample_probe_dict[i].split(",")
	for j in range(len(probe_name)):
		if probe_symbol[probe_name[j]] != None and probe_symbol[probe_name[j]] != "":
			gene_name = probe_symbol[probe_name[j]] #gene name <=> probe name
			gene_list[gene_name] += 1
			gene_count[gene_name] += 1
			gene_exp_sum[gene_name] += float(exp[j])
			#print(gene_count)
			j += 1
	for k in gene_count.keys():
		gene_exp_ave[k] = str(round(gene_exp_sum[k]/gene_count[k],1))
	sample_gene_exp[i] = gene_exp_ave
	#sample_gene_exp[i] = gene_exp_ave
sample_list.remove("ID_REF")
out_file.write("geneNames" + "\t" + "\t".join(sample_list) + "\n")
#for i in sample_list:
#	sample_gene_expr = sample_gene_exp[i]
for i in gene_list.keys():
	out_file.write(i + "\t")
	for j in range(len(sample_list) - 1):
		out_file.write(sample_gene_exp[sample_list[j]][i] + "\t")
	out_file.write(sample_gene_exp[sample_list[-1]][i] + "\n")
out_file.close()


