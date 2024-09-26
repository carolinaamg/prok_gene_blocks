import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO
import shutil

rearrangements_file = sys.argv[1]
core_genome = sys.argv[2]
genomes = sys.argv[3]
output_folder = sys.argv[4]

genome_list = list()
with open(rearrangements_file) as input_file:
	for line in input_file:
		line = line.rstrip()
		if line.startswith("Query"):
			pass
		else:
			tabs = line.split("\t")
			query = tabs [0]
			if query not in genome_list:
				genome_list.append(query)
input_file.close()

protein2genome = dict()
with open(core_genome) as core_genome:
	for line in core_genome:
		line = line.rstrip()
		tabs = line.split("\t")
		fam = tabs [0]
		info = tabs [1:len(tabs)]

		for i in info:
			item = i.split("&")
			genome = item[0].replace(".faa","")
			protein = item[1]
			protein2genome[protein] = genome
core_genome.close()

genome2seqs = defaultdict(list)
for file in os.listdir(genomes):
	if file.endswith(".aln"):
		file_path = os.path.join(genomes, file)
			
		fam_entries = dict()
		for record in SeqIO.parse(file_path, "fasta"):
			genome_name = protein2genome[str(record.id)]
			record.id = genome_name
			fam_entries[genome_name] = record
			seq_len = len(record.seq)

		for genome in genome_list:
			if genome in fam_entries:
				entry = str(fam_entries[genome].seq)
				genome2seqs[genome].append(entry)
			else:
				entry = "-"*seq_len
				genome2seqs[genome].append(entry)

name = rearrangements_file.split("/")
output_name = name[-1] + ".aln"
final_output_name = os.path.join(output_folder, output_name)
final_output_name = final_output_name.replace(".tsv","")
output_alignment = open(final_output_name, "w")
print(final_output_name)
for key,value in genome2seqs.items():
	final_seq = "".join(value)
	# for i in value:
	output_alignment.write(">"+key + "\n")
	output_alignment.write(final_seq + "\n")
output_alignment.close()

