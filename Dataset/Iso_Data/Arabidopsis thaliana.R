library("igraph")
library("bnlearn")
library("qgraph")

# import Expression dataset
dataset <- read.csv(file = 'ExpressionData.csv', row.names = 1)


# use a structure learning algorithm

pc <- pc.stable(dataset)


# visualize the graph
qgraph(pc, legend.cex = 0.5,
       asize=1,edge.color="black", vsize= 3)

