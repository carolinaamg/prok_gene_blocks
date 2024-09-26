import os, sys, re, subprocess, shlex
from collections import defaultdict

input_folder = sys.argv[1]

outfolder = input_folder
for i in os.listdir(input_folder):
	if i.endswith(".fna"): 
		fasta = os.path.join(input_folder, i)
		core = re.sub(".fna", "", i)

		prot = os.path.join(outfolder, core+".faa")
		nucl = os.path.join(outfolder, core+".genes.fna")

		cmd = "prodigal -i "+ fasta +" -a "+ prot +" -d "+ nucl
		print (i)
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("prodigal.out", "w"), stderr=open("prodigal.err", "w"))

print("DONE.")

