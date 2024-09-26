import os
import sys
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import more_itertools as mit

file = sys.argv[1] #families_core.txt"
folder_genomes = sys.argv[2] #Proteus_cibarius"
output = open(sys.argv[3], "w")

###I first find out what families/protein each genome has
all_genomes = list()
all_families = list()
core2genes = defaultdict(list)
core2fams = defaultdict(list)
with open(file) as file:
    for line in file:
        line = line.rstrip()
        tabs = line.split("\t")
        fam = tabs [0]
        genomes = tabs [1:len(tabs)]
        all_families.append(fam)

        for x in genomes:
            info = x.split("&")
            genome_name = info[0].replace(".faa","")
            protein = info[1]
            core2genes[genome_name].append(protein)
            core2fams[genome_name].append(fam)

            if genome_name not in all_genomes:
                all_genomes.append(genome_name)
file.close()

all_pairs = list()
all_pairs_proteins = defaultdict(list)
protein2fam = dict()
for x in all_genomes:
    for y in all_genomes:
        if x != y:
            pair_name = x+"~"+y
            all_pairs.append(pair_name)

            x_fams = core2fams[x]
            y_fams = core2fams[y]
            x_prots = core2genes[x]
            y_prots = core2genes[y]

            for fam in all_families:
                if fam in x_fams and fam in y_fams:
                    x_index = x_fams.index(fam)
                    y_index = y_fams.index(fam)

                    x_protein = x_prots[x_index]
                    y_protein = y_prots[y_index]

                    all_pairs_proteins[pair_name].append(x_protein)
                    all_pairs_proteins[pair_name].append(y_protein)
                    protein2fam[x_protein] = fam
                    protein2fam[y_index] = fam

species = folder_genomes.split("/")
species = species[-1]
###Analizing pairs
for pair in all_pairs:
    pair_name = pair.split("~")
    query = pair_name[0]
    ref = pair_name[1]
    pair_proteins = all_pairs_proteins[pair]

    query_genome = os.path.join(folder_genomes, query + ".genes.fna")
    query_genome_file = os.path.join(folder_genomes, query+ ".fna")
    ref_genome = os.path.join(folder_genomes, ref + ".genes.fna")
    
    local2actual = dict()
    global_count = 0
    for record in SeqIO.parse(query_genome_file, "fasta"):
        local_count = 0
        sequence = str(record.seq)
        for i in sequence:
            global_count = global_count+1
            local_count = local_count+1

            position_id = str(record.id) + "_" + str(local_count)
            local2actual[position_id] = global_count

    query_fams = list()
    ref_fams = list()
    fams2prots = dict()
    gene2realposition = dict()
    for record in SeqIO.parse(query_genome, "fasta"):
        if str(record.id) in all_pairs_proteins[pair]:
            info = record.description.split(" # ")
            start = info[1]
            end = info[2]
            contig_name = info[0].split("_")
            contig_name.pop()
            contig_name = "_".join(contig_name)

            actual_start = local2actual[contig_name + "_" + start]
            actual_end = local2actual[contig_name + "_" +end]

            # #Recording the actual position
            position = str(actual_start) + "_" + str(actual_end)
            gene2realposition[str(record.id)] = position
            
            fam_name = protein2fam[str(record.id)]
            query_fams.append(fam_name)
            fams2prots[fam_name] = str(record.id)

    for record in SeqIO.parse(ref_genome, "fasta"):
        if str(record.id) in all_pairs_proteins[pair]:
            fam_name = protein2fam[str(record.id)]
            ref_fams.append(fam_name)

    ###Checking whether the first fam in query is not the last gene in ref
    if query_fams[0] == ref_fams[-1]: 
        item = ref_fams[-1]
        ref_fams.insert(0,ref_fams.pop(ref_fams.index(item))) ##Changing the position of the last gene in ref

    order = list()
    for family in query_fams:
        ref_index = ref_fams.index(family)
        order.append(ref_index)

    all_blocks_translocated = list()
    all_blocks_inverted = list()
    reversed_families = list()
    for group in mit.consecutive_groups(order):
        block = list(group)
        if len(block) > 1:
            block = list(set(block))
            all_blocks_translocated.append(block)
        else:
            reversed_families.append(block[0])
                     
    for rearrangement in all_blocks_translocated:
        genes = list()
        families = list()
        for i in rearrangement:
            fam = ref_fams[i]
            protein = fams2prots[fam]
            genes.append(protein)
            families.append(fam)


        first_prot = gene2realposition[genes[0]]
        last_prot = gene2realposition[genes[-1]]

        all_genes = ";".join(genes)
        all_fams = ";".join(families)
        output.write(query + "\t" +ref+ "\t" + species + "\t" + "Translocation"+ "\t" + all_fams+ "\t" + all_genes + "\t" + first_prot+ "\t" + last_prot+ "\n")
        # print(query + "\t" +ref+ "\t" + species + "\t" + "Translocation"+ "\t" + all_fams+ "\t" + all_genes + "\t" + first_prot+ "\t" + last_prot)

    reversed_families.reverse()
    for group in mit.consecutive_groups(reversed_families):
        block = list(group)
        if len(block) > 1:
            all_blocks_inverted.append(block)
        else:
            pass

    for rearrangement in all_blocks_inverted:
        genes = list()
        families = list()
        for i in rearrangement:
            fam = ref_fams[i]
            protein = fams2prots[fam]
            genes.append(protein)
            families.append(fam)


        first_prot = gene2realposition[genes[0]]
        last_prot = gene2realposition[genes[-1]]

        all_genes = ";".join(genes)
        all_fams = ";".join(families)
        output.write(query + "\t" +ref+ "\t" + species + "\t" + "Inversion"+ "\t" + all_fams+ "\t" + all_genes + "\t" + last_prot+ "\t" + first_prot+ "\n")
        # print(query + "\t" +ref+ "\t" + species + "\t" + "Inversion"+ "\t" + all_fams+ "\t" + all_genes + "\t" + last_prot+ "\t" + first_prot)        

