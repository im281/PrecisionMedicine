#import networkx as nx
#import Chapter1_3 
#G=nx.Graph()
#G.add_node("spam")
#G.add_edge(1,2)
#G.add_edge(2,3)
#G.add_edge(3,4)
#print(list(G.nodes()))
#print(list(G.edges()))
l = []
a = 'ATGCG'
b = 'GCATG'
c = 'CATGC'
d = 'AGGCA'
e = 'GGCAT'
l.append(a)
l.append(b)
l.append(c)
l.append(d)
l.append(e)

import networkx as nx
import Chapter1_3 
import chpt4

#dna = Chapter1_3.readSequenceArray("test.txt")
#dna = Chapter1_3.readFasta("test.txt")
#g = chpt4.overlap(l)
#list = []
#for i in g.edges():
#    print(i[0] + " -> " + i[1])
#    list.append(i[0] + " -> " + i[1])
#Chapter1_3.WriteResults(list)





#s = 'AAGATTCTCTAC'
#s = 'TAATGCCATGGGATGTT'
#s = 'AAGATTCTCTAC'
#dna = Chapter1_3.kmersFromDNA(s,3)
#g = deBruijn.construct_graph(dna,2)
#v = deBruijn.output_contigs(g)




d = Chapter1_3.readFasta("test.txt")
s = d.get("test")

g2 = chpt4.debruijn(12,s)
n = g2.nodes()
e = g2.edges()


list = []
for i in g2.edge.items():
    temp = []
    for j in i[1]:
        temp.append(j)    
    myString = ",".join(temp)
    if myString != "":
        list.append(i[0] + " -> " + myString)
        print(i[0] + " -> " + myString)
        
Chapter1_3.WriteResults(list)













#G = Chapter1_3.kmerGraph(dna,3)
#print(list(G.nodes()))
#print(list(G.edges()))

#dna = Chapter1_3.readFasta("test.txt")
#s = dna.get("test")
#G = Chapter1_3.kmerdeBruijnGraph(dna,3)

#edgeList = []
#for i in G.edges():
#    edgeList.append(i[0] + " -> " + i[1])

#Chapter1_3.WriteResults(edgeList)
