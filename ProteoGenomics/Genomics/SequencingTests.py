import networkx as nx
import Chapter1_3 
import chpt4

a = 'ATGCG'
b = 'GCATG'
c = 'CATGC'
d = 'AGGCA'
e = 'GGCAT'
list = []
list.append(a)
list.append(b)
list.append(c)
list.append(d)
list.append(e)


#Rosalind test
#s = 'AAGATTCTCTAC'

#book test
s = 'TAATGCCATGGGATGTT'
kmers = Chapter1_3.kmersFromDNA(s,3)

#long DNA string test
#d = Chapter1_3.readFasta("test.txt")
#s = d.get("test")

g2 = chpt4.debruijn(12,s)

#save results to file
#chpt4.savedebruijn(g2)

#s = Chapter1_3.readSequenceArray("test.txt")
#g = chpt4.overlap(list)
chpt4.assembleGenome(g2)
chpt4.assembleGenome(
#chpt4.saveOverlapGraph(g)




 
















