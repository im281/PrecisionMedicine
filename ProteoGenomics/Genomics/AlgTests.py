from AlgsSedgewickWayne.Bag import *
from AlgsSedgewickWayne.Digraph import *

b = Bag()

b.add("test")
b.add("Another test")
test  = ""

g = Digraph(5)
g.addEdge("test","test2")
g.addEdge(2,3)
t = g.get_edges()
