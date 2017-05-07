# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:32:26 2017

@author: Owner
"""
import Chapter1_3 
#http://rosalind.info/problems/

#1A 
text = 'GCGCG'
pattern = 'GCG'
r = Chapter1_3.patternCount(text,pattern)

#1B
text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
k = 4
Chapter1_3.frequentWords(text,k)
Chapter1_3.fasterFrequentWords(text,k)
#1C
text = 'AAAACCCGGT'
Chapter1_3.reverseComplement(text)

#1D
text = 'GATATATGCATATACTT'
p = 'ATAT'
Chapter1_3.naive(p,text)

#1E
text = 'CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC'
Chapter1_3.clumpFinding(text,5,4,75)

#1F
text = 'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG'
Chapter1_3.skew(text)

#1G
s1 = 'GGGCCGTTGGT'
s2 = 'GGACCGTTGAC'
Chapter1_3.hammingDistance(s1,s2)

#1H
text = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC'
p = 'ATTCTGGA'
Chapter1_3.approximatePatternMatch(text,p,3)

#1I
text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
Chapter1_3.frequentWordsWithMismatches(text,4,1)

#1J
text = ''

#1K

#1L

#1M

#1N

#2A
s1 = 'ATTTGGC'
s2 = 'TGCCTTA'
s3 = 'CGGTATC'
s4 = 'GAAAATT'
l = []
l.append(s1)
l.append(s2)
l.append(s3)
l.append(s4)
Chapter1_3.motifEnumeration(l,3,1)

#2B
s1 = 'AAATTGACGCAT'
s2 = 'GACGACCACGTT'
s3 = 'CGTCAGCGCCTG'
s4 = 'GCTGAAGCCCGG'
s5 = 'AGTACGGGACAG'
l = []
l.append(s1)
l.append(s2)
l.append(s3)
l.append(s4)
l.append(s5)
Chapter1_3.medianString(l,3)


#2H
s1 = 'TTACCTTAAC'
s2 = 'GATATCTGTC'
s3 = 'ACGGCGTTCG'
s4 = 'CCCTAAAGAG'
s5 = 'CGTCAGAGGT'
l = []
l.append(s1)
l.append(s2)
l.append(s3)
l.append(s4)
l.append(s5)
Chapter1_3.medianString(l,3)
Chapter1_3.distanceBetweenPatternAndStrings('AAA',l)

#3A
g = Chapter1_3.readGenome('rosalind_ba3a.txt')
kmers = Chapter1_3.kmerCompositionLexigographic(g,50)
kmers

#3B
reads = Chapter1_3.getSequencingReads('rosalind_ba3b.txt')
t = Chapter1_3.genomePath(reads)
print(t)

#for k in r:
#	    print(k)
#     
#for k, v in r.items():
#	    print(k,v)

import networkx as nx
G=nx.Graph()
G.add_node("spam")
G.add_edge(1,2)
G.add_edge(2,3)
print(list(G.nodes()))
print(list(G.edges()))

dna = 'TAATGCCATGGGATGTT'
d = Chapter1_3.kmersFromDNA(dna,3)

import networkx as nx
G=nx.Graph()
for i in d:
    G.add_node(i)
    for j in d:
        if Chapter1_3.suffix(i) == Chapter1_3.prefix(j):
            G.add_edge(i,j)
            

            



