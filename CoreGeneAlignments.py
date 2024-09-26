import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO
import shutil

rearrangements_file = sys.argv[1]
core_genome = sys.argv[2]
genomes = sys.argv[3]

fam2genomes = defaultdict(list)
fam2proteins = defaultdict(list)
protein2fam = dict()
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
			fam2genomes[fam].append(genome)
			fam2proteins[fam].append(protein)
			protein2fam[protein] = fam
core_genome.close()

genome_list = list()
with open(rearrangements_file) as rearrangements_file:
	for line in rearrangements_file:
		line = line.rstrip()
		if line.startswith("Query"):
			pass
		else:
			tabs = line.split("\t")
			query = tabs [0]
			ref = tabs[1]
			if query not in genome_list:
				genome_list.append(query)
			if ref not in genome_list:
				genome_list.append(rearrangements_file)			
rearrangements_file.close()

genome2proteins = defaultdict(list)
for genome in genome_list:
	for k,v in fam2genomes.items():
		proteins = fam2proteins[k]
		if genome in v:
			genome_ind = v.index(genome)
			genome_protein = proteins[genome_ind]
			genome2proteins[genome].append(genome_protein)
		else:
			pass

fam2sequences = defaultdict(list)
for key, value in genome2proteins.items():
	genome_path = os.path.join(genomes, key + ".genes.fna")
	protein_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))

	for i in value:
		record = protein_dict[i]
		family = protein2fam[i]
		fam2sequences[family].append(record)

for k,v in fam2sequences.items():
	output = os.path.join(genomes, k + ".fna")
	alignment = os.path.join(genomes, k + ".aln")
	if len(v)>=5:
		SeqIO.write(v, output, "fasta")
		
		cmd = "mafft --auto "+ output
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open(alignment, "a"), stderr=open("log_file.txt", "a"))
		print(cmd)		
