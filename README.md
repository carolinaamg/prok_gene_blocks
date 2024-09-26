# prok_gene_blocks

# Hi, there! 👋  
This repository contains the code used to identify blocks of core genes for the publication **"Prevalence and Dynamics of Genome Rearrangements in Bacteria
and Archaea"** by Carolina A. Martinez-Gutierrez and Louis-Marie Bobay. Note that the manuscript was just submitted! Before identifying gene blocks, you first have to find the core genome of your microbial species. To do so, you need to install and run *Corecruncher*. The program will take yout fasta files containing predicted gene amino acid sequences (e.g., output files of Prodigal). 

The authors are working on getting the code developed neat and tidy to be more user fiendly than what it currently is 😀   

➡️ Example of how to run *Prodigal*:  
⁉️ *Prodigal* information:   
Manuscript for citation: https://link.springer.com/article/10.1186/1471-2105-11-119  
Installation and how to run: https://github.com/hyattpd/Prodigal

I wrote a little script to run *Prodigal* on your input files:
```
python RunProdigal.py my_species_folder 
```
**my_species_folder** is the folder in which your assembled genomes are **(with the extension .fna)**. Output files will go there too.    
You will get two output files for each assembly: .genes.fna will have predicted genes with nucleotide sequences and .faa will have predicted genes with amino acid sequences.

➡️ Example of how to run *Corecruncher*:  

⁉️ *Corecruncher* information:    
Manuscript for citation: https://academic.oup.com/mbe/article-abstract/38/2/727/5901545  
Installation and how to run: https://github.com/lbobay/CoreCruncher  

```
python corecruncher_master.py -ext .faa -freq 90 -score 80 -in  my_species_folder -out my_species_core_genome
```
**-freq**: is the minimum frequency of a gene in your species to be considered a core gene.   
**-score**: is the minimum score of a given core gene to be considered.   
**my_species_core_genome** will have your output core genome!  

Now, let's get your blocks of core genes identified...
Note that this code will do a pairwise analysis, meaning that it will compare everything vs everything to identify your blocks of core genes. 

➡️ Run the code FindRearrangementsProk.py as follows:   

```
python FindBlocksProk.py families_core.txt my_species_folder output_file_blocks.tsv
```
**families_core.txt** is the ouuput file of *Corecruncher* and it relates each protein to a family.   
**my_species_folder** is the folder in which you have your genome files. Make sure you also have the assembly files in it!    
**output_file_blocks.tsv** will contain all the information of the blocks identified.     

If you would like to build a core genome alignment to make a phylogenomic tree and ask questions regarding the evolution of your blocks, *we got you covered*!    

➡️ Run the code CoreGeneAlignments.py as follows:   

```
python CoreGeneAlignments.py output_file_blocks.tsv families_core.txt my_species_folder
```
**output_file_blocks.tsv** is the output file of FindRearrangementsProk.py. Note that ConcatenatedAlignment.py will only include in the alignments the genomes listed in the first column.    
**families_core.txt** is the ouuput file of *Corecruncher* and it relates each protein to a family.   
**my_species_folder** is the folder in which you have your genome files. This is were your core gene alignments will be stored.    

⁉️ CoreGeneAlignments will use MAFFT to align sequences (make sure you have it on your path!).   
Publication here: https://academic.oup.com/nar/article/33/2/511/2549118    

➡️ Run the code ConcatenatedAlignment.py as follows:   

```
python ConcatenatedAlignment.py output_file_blocks.tsv families_core.txt my_species_folder individual_alignments_output_folder
```
**output_file_blocks.tsv** is the output file of FindRearrangementsProk.py.    
**families_core.txt** is the ouuput file of *Corecruncher* and it relates each protein to a family.    
**my_species_folder** is the folder in which you have your genome files. This is were your core gene alignments will be stored during the previous step.    
**individual_alignments_output_folder** is the folder in which you want your output concatenated alignment.     

### Note that the code provided requires multiple python libraries: subprocess, shlex, collections, Biopython SeqIO, and shutil. 
