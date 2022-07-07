#Transform a given pedigree into a graph
## pedigree contains first 3 columns as "ID", "Dam", and "Sire"
## Additional columns in pedigree denote attributes of graph nodes
## `attributes` argument takes character vector of pedigree column names 
pedToGraph <- function(pedigree, attributes = NULL)
{
  if(!is.null(attributes))
  {
    attrCol <- match(attributes, names(pedigree))
    if(length(attrCol) > 0) nodeAtt <- pedigree[, c(1, attrCol)]
  }
  
  dnmiss <- which(!is.na(pedigree[, 2]))
  snmiss <- which(!is.na(pedigree[, 3]))
  
  relations <- data.frame(parent = c(pedigree[dnmiss, 2],
                                     pedigree[snmiss, 3]),
                          id = pedigree[c(dnmiss, snmiss), 1])
  graph <- graph_from_data_frame(relations, directed = TRUE,
                                 vertices = nodeAtt)
  return(graph)
}
# example
## library(nadiv)
## library(igraph)
## pedToGraph(FG90, attributes = "sex")


#return the average connectivity of a graph
avgEdgeCon <- function(graph)
{
  vcnt <- vcount(graph)
  for (i in 1:(vcnt-1))
  { 
    sm <- sum(sapply(seq.int(i + 1, vcnt, 1), FUN = edge_connectivity,
      graph = graph, source = i, checks = TRUE))
  }
  
  avgCon <- sm / choose(vcnt, 2)
  return(avgCon)
}
# example
## library(nadiv)
## library(igraph)
## g <- pedToGraph(FG90, attributes = "sex")
## avgEdgeCon(g)




#return the variance of the connectivity of a graph
varEdgeCon <- function(graph, mean = avgEdgeCon(graph), sum = 0)
{
  for (i in 1:(vcount(graph)-1))
  {
    for (j in (i+1):vcount(graph))
    {
      sum <- sum + (edge_connectivity(graph, source = i, target = j)
                   - mean)**2
    }
  }
  
  varCon <- sum / (choose(vcount(graph), 2) - 1)
  return(varCon)
}



#returns the average reproductive success of a pedigree graph
avgRepSucs <- function(graph, nodes = NULL)
{
  if(is.null(nodes)) nodes <- seq.int(1, vcount(graph), 1) 
  
  sm <- sapply(nodes, FUN = degree,
    graph = graph, mode = "out", loops = FALSE, normalized = FALSE)
  
  return(sum(sm) / length(nodes))
}
# example
## library(nadiv)
## library(igraph)
## g <- pedToGraph(FG90, attributes = "sex")
## ## calculate across all nodes in g
## avgRepSucs(g)
## ## only calculate over first 3 generations of pedigree
## gen <- genAssign(FG90)
## inclNodes <- match(FG90[gen < 3, 1], vertex_attr(g, name = "name"))
## avgRepSucs(g, nodes = inclNodes)



