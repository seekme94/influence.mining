#' This function inputs a graph object, percentage and a seed method and returns k% nodes as seed using given method
#' @name select_seed
#' @param graph is the igraph object
#' @param budget number of influential nodes to be fetched. Default value is 1
#' @param seed_method is the selection method for seed (initial nodes). Value can be "random", "degree", "closeness", "betweenness", "coreness", "eigenvector", "a-degree", "a-closeness", "a-betweenness", "a-coreness", "a-eigenvector"
#' @return set of nodes
select_seed <- function (graph, budget, seed_method=c("random", "degree", "closeness", "betweenness", "coreness", "eigenvector", "a-degree", "a-closeness", "a-betweenness", "a-coreness", "a-eigenvector")) {
  nodes <- V(graph)
  seed <- NULL
  # Calculate the actual number of nodes to select as seed
  if (seed_method == "random") {
    # Select random nodes
    seed <- as.vector(sample(x=nodes, budget))
  }
  else if (seed_method == "degree") {
    # Calculate in/out degrees of all nodes
    degrees <- degree(graph, V(graph), mode="all", loops=FALSE, normalized=FALSE)
    data <- as.data.frame(degrees)
    data$nodes <- rownames(data)
    seed <- tail(data[order(data$degrees),], budget)$nodes
  }
  else if (seed_method == "closeness") {
    # Calculate in/out closeness of all nodes, normalized between 0 and 1
    closenesses <- closeness(graph, V(graph), mode="all", normalized=TRUE)
    data <- as.data.frame(closenesses)
    data$nodes <- rownames(data)
    seed <- tail(data[order(data$closenesses),], budget)$nodes
  }
  else if(seed_method == "betweenness") {
    # Calculate betweenness centrality (for huge data sets, use betweenness.estimate() and give some max value of path length as cutoff)
    betweennesses <- betweenness(graph, V(graph), directed=FALSE)
    data <- as.data.frame(betweennesses)
    data$nodes <- rownames(data)
    seed <- tail(data[order(data$betweennesses),], budget)$nodes
  }
  else if (seed_method == "coreness") {
    # Calculate in/out closeness of all nodes
    coreness <- graph.coreness(graph, mode="all")
    data <- as.data.frame(coreness)
    data$nodes <- rownames(data)
    seed <- tail(data[order(data$coreness),], budget)$nodes
  }
  else if (seed_method == "eigenvector") {
    # Calculate eigenvectors of the graph
    eigen <- evcent(graph, directed=FALSE)
    eigenvectors <- eigen$vector
    eigenvalue <- eigen$value
    data <- as.data.frame(eigenvectors)
    data$nodes <- rownames(data)
    # Select nodes with highest eigenvector centralities
    seed <- tail(data[order(data$eigenvectors),], budget)$nodes
  }
  else if (seed_method %in% c("a-degree", "a-closeness", "a-betweenness", "a-coreness", "a-eigenvector")) {
    seed <- select_adaptive_seed(graph, budget, seed_method)
  }
  as.numeric(seed)
}

#' This function inputs an igraph object and returns budget% nodes as seed using the adaptive method provided
#' @name select_adaptive_seed
#' @param graph is the igraph object
#' @param budget defines what percentage of most influential nodes out of all nodes is required as output. Default value is 1
#' @param seed_method is the selection method for seed (initial nodes). Value can be a-degree", "a-closeness", "a-betweenness", "a-coreness", "a-eigenvector"
#' @return seed object containing set of vertices
select_adaptive_seed <- function (graph, budget, seed_method=c("a-degree", "a-closeness", "a-betweenness", "a-coreness", "a-eigenvector")) {
  # Add node labels (or the graph resets)
  V(graph)$name <- V(graph)
  seed <- NULL
  # Calculate the actual number of nodes to select as seed
  graph <- largest_component(graph)
  for (i in 1:budget) {
    param <- NULL
    if (seed_method == "a-degree") {
      param <- degree(graph, V(graph), mode="all")
    }
    else if (seed_method == "a-closeness") {
      param <- closeness(graph, V(graph), mode="all")
    }
    else if(seed_method == "a-betweenness") {
      param <- betweenness(graph, V(graph))
    }
    else if (seed_method == "a-coreness") {
      param <- graph.coreness(graph, mode="all")
    }
    else if (seed_method == "a-eigenvector") {
      # Calculate eigenvectors of the graph
      eigen <- evcent(graph, directed=FALSE)
      param <- eigen$vector
    }
    max_node <- which.max(param)
    seed <- c(seed, V(graph)[max_node]$name)
    graph <- delete.vertices(graph, max_node)
    graph <- largest_component(graph)
  }
  seed
}

#' @title Returns CI value of given graph and node
#' @name collective_influence
#' @param graph the igraph object
#' @param neighborhood_distance is the distance to which the neighborhood nodes are searched for
#' @param node_id is the ID of the target node
#' @param method is the metric to calculate sum of influence. Default is "degree"
#' @return influence as product of degree of target node and total sum of degrees of neighborhood
#' TODO: extend the function and include adaptive methods as well as other centrality methods
#' @import igraph
#' @export
collective_influence <- function(graph, neighborhood_distance, node_id, method=c("degree")) {
  neighbors_at_distance <- neighborhood(graph, neighborhood_distance, nodes=node_id, mode="all")[[1]]
  neighbors_at_distance_discount <- neighborhood(graph, neighborhood_distance - 1, nodes=node_id, mode="all")[[1]]
  # find all the nodes lying at given distance
  neighbors_only_at_distance <- setdiff(neighbors_at_distance, neighbors_at_distance_discount)
  # calculate the degree of all the nodes lying at distance
  degrees <- degree(graph, neighbors_only_at_distance)
  # convert the list result into vector
  degree_sum <- as.vector(degrees)
  # subtract one from each degree, sum the result and return
  total_sum <- sum(degree_sum - 1)
  node_degree <- (degree(graph,node_id)[[1]]) - 1
  ans <- node_degree * total_sum
  ans
}

#' Returns ranks from 1 to highest rank before the graph is discontinued, for nodes in given graph using adaptive variation of given method
#' @name get_adaptive_ranking
#' @param g the igraph object
#' @param method the adaptive method to use. Allowed values are degree, closeness, betweenness, coreness, eigenvector, eccentricity, pagerank, collective_influence
#' @return vector of ranks
get_adaptive_ranking <- function(g, method=c("degree", "closeness", "betweenness", "coreness", "eigenvector", "eccentricity", "pagerank", "collective_influence")) {
  V(g)$name <- V(g)
  V(g)$rank <- -1
  current_rank <- 1
  graph <- g
  while (TRUE) {
    graph <- largest_component(graph)
    param <- 0
    if (method == "a-degree") {
      param <- degree(graph, V(graph), mode="all")
    } else if (method == "a-closeness") {
      param <- closeness(graph, V(graph), mode="all")
    } else if (method == "a-betweenness") {
      param <- betweenness(graph, V(graph))
    } else if (method == "a-coreness") {
      param <- graph.coreness(graph, mode="all")
    } else if (method == "a-eigenvector") {
      param <- evcent(graph, directed=FALSE)$vector
    } else if (method == "a-eccentricity") {
      param <- eccentricity(graph)
    } else if (method == "pagerank") {
      param <- page_rank(graph, directed=TRUE)$vector
    } else if (method == "collective_influence") {
      param <- sapply(V(graph), function(x) { collective_influence(graph, neighborhood_distance=2, x) })
    }
    max_nodes <- which.max(param)
    V(g)[V(g)$name == V(graph)[max_nodes]$name]$rank <- current_rank
    graph <- delete.vertices(graph, max_nodes)
    current_rank <- current_rank + 1
    if (vcount(graph) <= 1) {
      break
    }
  }
  V(g)$rank[V(g)$rank == -1] <- current_rank
  V(g)$rank
}
