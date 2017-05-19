import networkx as nx
import Chapter1_3 
import chpt4







#s = 'AAGATTCTCTAC'
#s = 'TAATGCCATGGGATGTT'
#s = 'AAGATTCTCTAC'
#dna = Chapter1_3.kmersFromDNA(s,3)
#g = deBruijn.construct_graph(dna,2)
#v = deBruijn.output_contigs(g)



d = Chapter1_3.readFasta("test.txt")
s = d.get("test")

g2 = chpt4.debruijn(12,s)
chpt4.savedebruijn(g2)

n = g2.nodes()
e = g2.edges()

def debruijn(g):
    list = []
    for i in g.edge.items():
        temp = []
        for j in i[1]:
            temp.append(j)    
            myString = ",".join(temp)
            if myString != "":
                list.append(i[0] + " -> " + myString)
                print(i[0] + " -> " + myString)     
    Chapter1_3.WriteResults(list)












