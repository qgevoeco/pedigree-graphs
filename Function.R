#Transform a given pedigree into a graph
pedToGraph <- function(pedigree)
{
  nodeAtt <- pedigree[, c("id", "sex")]
  dnmiss <- which(!is.na(pedigree[, 2]))
  snmiss <- which(!is.na(pedigree[, 3]))
  
  relations <- data.frame(parent = c(pedigree[dnmiss,2],
                                     pedigree[snmiss, 3]),
                          id = pedigree[c(dnmiss, snmiss), 1])
  graph <- graph_from_data_frame(relations, directed = TRUE,
                                 vertices = nodeAtt)
  return(graph)
}

#return the average connectivity of a graph
avgEdgeCon <- function(graph, sum = 0)
{
  for (i in 1:(vcount(graph)-1))
  {
    for (j in (i+1):vcount(graph))
    {
      sum = sum + edge_connectivity(graph, source = i, target = j)
    }
  }
  
  avgCon = sum / choose(vcount(graph), 2)
  return(avgCon)
}

#return the variance of the connectivity of a graph
varEdgeCon <- function(graph, mean = avgEdgeCon(graph), sum = 0)
{
  for (i in 1:(vcount(graph)-1))
  {
    for (j in (i+1):vcount(graph))
    {
      sum = sum + (edge_connectivity(graph, source = i, target = j)
                   - mean)**2
    }
  }
  
  varCon = sum / (choose(vcount(graph), 2) - 1)
  return(varCon)
}

#returns the average reproductive success of a pedigree graph
avgRepSucs <- function(graph, sum = 0, count = 0)
{
  for (i in 1:vcount(graph))
  {
    degV <- degree(graph, v = i, mode = "out", loops = FALSE)
    if (degV == 0)
    {
      next
    }
    
    sum = sum + degV
    count = count + 1
  }
  
  return(sum / count)
}
