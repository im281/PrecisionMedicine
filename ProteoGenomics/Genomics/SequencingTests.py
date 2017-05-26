import networkx as nx
import Chapter1_3 
import chpt4
import debruijn

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




#book test
s = 'TAATGCCATGGGATGTT'
kmers = Chapter1_3.kmersFromDNA(s,3)
g3 = debruijn.construct_graph(kmers,2)

g = debruijn.output_contigs(g3)





 
















