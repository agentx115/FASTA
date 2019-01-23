# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 15:25:31 2019

@author: Rhiannon
"""
#%%
from Bio import SeqIO
#from Bio.Seq import Seq
#%%

#read in sequence as fasta
fname = "sequence.txt"
dna = str(SeqIO.read(fname, "fasta").seq)
#%%
#create codons list

codons = {}
fcodons = "codons.txt"
open_codons = open(fcodons, "r")
codon_content = open_codons.read()
open_codons.close()

codon_list = codon_content.split("\n")

for i in range(0, len(codon_list)-1, 2):
    codons[codon_list[i]] = codon_list[i+1]
#%%
#get reverse strand
#own way:
reverse_dna = dna.replace("a", "x")
reverse_dna = reverse_dna.replace("t", "a")
reverse_dna = reverse_dna.replace("x", "t")
reverse_dna = reverse_dna.replace("c", "d")
reverse_dna = reverse_dna.replace("g", "c")
reverse_dna = reverse_dna.replace("d", "g")
reverse_dna = reverse_dna[::-1]

#BioPy ver:
#much shorter!!
#dna1 = Seq(dna)
#reverse_dna = dna1.reverse_complement()

    
 #%%   
#translate sequence

def translate_seq(dna_seq):
    trans =""
    for triplet in range(0, len(dna_seq), 3):
        code = dna_seq[triplet:(triplet+3)]
        if code in codons:
            trans = trans + codons[code]
        else:
            trans = trans + "-"
    return trans

#%%
#get reading frames

rf_list = []

for i in range(0,3):
    rf_list.append(dna[i:])
    rf_list.append(reverse_dna[i:])

for sequence in rf_list:
    translated = translate_seq(sequence)
    if rf_list.index(sequence)%2 == 0:
        print("Forwards Strand")
    else:
        print("Reverse Strand")
    print("DNA seq: ",sequence)
    print("Translated seq: ", translated)
    print("\n")
