#' @title Setup logging
#' @name setup
#' @description This function should be called before any other
#' @param logging flag to enable loggging. Default is TRUE
#' @import logging
setup <- function(logging=TRUE) {
  if (logging) {
    basicConfig()
    addHandler(writeToFile, logger="influence_maximization", file="output.log")
  }
}

#' @title Mine the most influential nodes in given graph
#' @name influence
#' @param graph is the igraph object
#' @param budget number of influential nodes to be fetched. Default value is 1
#' @param prob probability at which a node influences its neighbours
#' @param steps is the time steps for which, the diffusion process should run. Provide NULL for exhaustive run. Default value is 1
#' @param test_method specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"
#' @param heuristic specifies the heuristic method used for influence calculation. Required only when optimal_solution is FALSE
#' @param centrality_method is the centrality algorithm to use when heuristic is "CENTRALITY" or "ADAPTIVE_CENTRALITY". Value must be "DEGREE", "BETWEENNESS", "CLOSENESS" or "EIGENVECTOR"
#' @param parallel when true, executes the funtion using multiple CPU cores. Default value is TRUE
#' @param optimal_solution should be TRUE if influential nodes are to be derived using optimal algorithm. Caution! This is the slowest apporach
#' @param logging when true, a complete log is stored in output.log file
#' @return object containing: 1. Vector of influential nodes. 2. Measure of influence. 3. Elapsed time in seconds.
#' @import logging igraph
#' @export
influence <- function(graph, budget=1, prob=0.5, steps=1, optimal_solution=FALSE,
                      test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"),
                      heuristic=c("GREEDY", "PAGERANK", "COLLECTIVE_INFLUENCE", "CORENESS", "CENTRALITY", "ADAPTIVE_CENTRALITY"),
                      centrality_method=c("DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR"),
                      parallel=TRUE, logging=TRUE) {
  output <- NULL
  if (logging) {
    loginfo(paste("influence function parameters: graph_size=", vcount(graph), ", budget=", budget, ", prob=", prob,
                  ", steps=", steps, ", test_method=", test_method, ", parallel=", parallel, ", optimal_solution=", optimal_solution, sep=''))
  }
  if (optimal_solution) {
    output <- optimal_influential(graph=graph, budget=budget, prob=prob, test_method=test_method, parallel=parallel)
  } else {
    if (heuristic == "GREEDY") {
      output <- greedy_influential(graph=graph, budget=budget, steps=steps, test_method=test_method, prob=prob)
    } else if (heuristic == "PAGERANK") {
      output <- pagerank_influential(graph=graph, budget=budget, test_method=test_method)
    } else if (heuristic == "COLLECTIVE_INFLUENCE") {
      output <- collective_influence_influential(graph=graph, budget=budget, test_method=test_method)
    } else if (heuristic == "CORENESS") {
      output <- coreness_influential(graph=graph, budget=budget, test_method=test_method)
    } else if (heuristic == "CENTRALITY") {
      output <- centrality_influential(graph=graph, budget=budget, test_method=test_method, centrality_method=centrality_method)
    } else if (heuristic == "ADAPTIVE_CENTRALITY") {
      output <- adaptive_centrality_influential(graph=graph, budget=budget, test_method=test_method, centrality_method=centrality_method)
    }
  }
  if (logging) {
    loginfo(paste("Elapsed time (milliseconds)", output$time))
  }
  output
}

#' @title Returns the most influential nodes in a graph using adaptive centrality-based heuristics
#' @name adaptive_centrality_influential
#' @param graph is the igraph object
#' @param budget number of influential nodes to be fetched. Default value is 1
#' @param test_method specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"
#' @param centrality_method defines the centrality method to be used. Value must be:
#' @return object containing: 1. Vector of influential nodes. 2. Measure of influence. 3. Elapsed time in seconds.
#' @import igraph
#' @export
adaptive_centrality_influential <- function(graph, budget=1, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"),
                                   centrality_method=c("DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR")) {
  start <- as.numeric(Sys.time())
  # Preserve original graph as this object will be overwritten
  g <- graph
  V(graph)$name <- V(graph)
  influential_nodes <- NULL
  # Calculate the actual number of nodes to select
  for (i in 1:budget) {
    # Get the node with highest score
    max_node <- which.max(get_centrality_scores(g, centrality_method=centrality_method))
    influential_nodes <- c(influential_nodes, V(g)[max_node]$name)
    g <- delete.vertices(g, max_node)
    g <- largest_component(g)
    # Break if the graph is already disconnected
    if (vcount(g) == 1) {
      break
    }
  }
  end <- as.numeric(Sys.time())
  output <- NULL
  output$influential_nodes <- V(graph)[as.numeric(influential_nodes)]
  output$influence <- get_influence(graph, output$influential_nodes, test_method=test_method)
  output$time <- (end - start)
  output
}

#' @title Returns the most influential nodes in a graph using centrality-based heuristics
#' @name centrality_influential
#' @param graph is the igraph object
#' @param budget number of influential nodes to be fetched. Default value is 1
#' @param test_method specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"
#' @param centrality_method defines the centrality method to be used. Value must be:
#' @return object containing: 1. Vector of influential nodes. 2. Measure of influence. 3. Elapsed time in seconds.
#' @import igraph utils
#' @export
centrality_influential <- function(graph, budget=1, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"),
                                   centrality_method=c("DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR")) {
  start <- as.numeric(Sys.time())
  # Get Collective influence score of all nodes
  centrality <- get_centrality_scores(graph, centrality_method=centrality_method)
  x <- data.frame(centrality=centrality)
  x$node_id <- rownames(x)
  # Get budget nodes
  influential <- tail(x[order(x$centrality),], budget)$node_id
  end <- as.numeric(Sys.time())
  output <- NULL
  output$influential_nodes <- V(graph)[as.numeric(influential)]
  output$influence <- get_influence(graph, output$influential_nodes, test_method=test_method)
  output$time <- (end - start)
  output
}

#' @title Returns the most influential nodes in a graph using Collective influence heuristic
#' @name collective_influence_influential
#' @param graph is the igraph object
#' @param budget number of influential nodes to be fetched. Default value is 1
#' @param test_method specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"
#' @return object containing: 1. Vector of influential nodes. 2. Measure of influence. 3. Elapsed time in seconds.
#' @import igraph utils
#' @export
collective_influence_influential <- function(graph, budget=1, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC")) {
  start <- as.numeric(Sys.time())
  # Get Collective influence score of all nodes
  ci <- sapply(V(graph), function(x) { collective_influence(graph, neighborhood_distance=2, x) })
  x <- data.frame(ci=ci)
  x$node_id <- rownames(x)
  # Get budget nodes
  influential <- tail(x[order(x$ci),], budget)$node_id
  end <- as.numeric(Sys.time())
  output <- NULL
  output$influential_nodes <- V(graph)[as.numeric(influential)]
  output$influence <- get_influence(graph, output$influential_nodes, test_method=test_method)
  output$time <- (end - start)
  output
}

#' @title Returns the most influential nodes in a graph using Pagerank heuristic
#' @name pagerank_influential
#' @param graph is the igraph object
#' @param budget number of influential nodes to be fetched. Default value is 1
#' @param test_method specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"
#' @return object containing: 1. Vector of influential nodes. 2. Measure of influence. 3. Elapsed time in seconds.
#' @import igraph utils
#' @export
pagerank_influential <- function(graph, budget=1, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC")) {
  start <- as.numeric(Sys.time())
  # Get Pagerank of all nodes
  pagerank <- page_rank(graph)$vector
  x <- data.frame(pagerank=pagerank)
  x$node_id <- rownames(x)
  # Get budget nodes
  influential <- tail(x[order(x$pagerank),], budget)$node_id
  end <- as.numeric(Sys.time())
  output <- NULL
  output$influential_nodes <- V(graph)[as.numeric(influential)]
  output$influence <- get_influence(graph, output$influential_nodes, test_method=test_method)
  output$time <- (end - start)
  output
}

#' @title Returns the most influential nodes in a graph using Coreness heuristic
#' @name coreness_influential
#' @param graph is the igraph object
#' @param budget number of influential nodes to be fetched. Default value is 1
#' @param test_method specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"
#' @return object containing: 1. Vector of influential nodes. 2. Measure of influence. 3. Elapsed time in seconds.
#' @import igraph utils
#' @export
coreness_influential <- function(graph, budget=1, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC")) {
  start <- as.numeric(Sys.time())
  # Get coreness of all nodes
  coreness <- graph.coreness(graph, mode="all")
  # Get most core nodes
  influential <- V(graph)[which(coreness == max(coreness))]
  # If the number exceeds the given budget, then pick top degree nodes within influential
  if (length(influential) > budget) {
    x <- data.frame(degree=degree(graph, influential))
    x$node_id <- rownames(x)
    # Get budget nodes
    influential <- tail(x[order(x$degree),], budget)$node_id
  }
  end <- as.numeric(Sys.time())
  output <- NULL
  output$influential_nodes <- V(graph)[as.numeric(influential)]
  output$influence <- get_influence(graph, output$influential_nodes, test_method=test_method)
  output$time <- (end - start)
  output
}

#' @title Implements Greedy algorithm for Influence Maximization
#' @name greedy_influential
#' @param graph is the igraph object
#' @param budget number of influential nodes to be fetched. Default value is 1
#' @param steps is the time steps for which, the diffusion process should run. If exhaustive run is required, provide a high value (like 100). Default value is 1
#' @param prob is the probability of activation of a neighbour node. Default is 0.5
#' @param test_method specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"
#' @return output containing summary
#' @examples {greedy_influential(erdos.renyi.game(500, 0.005), budget=5, steps=99, prob=0.5, "RESILIENCE")}
#' @import igraph
#' @export
greedy_influential <- function(graph, budget, steps, prob=0.5, test_method) {
  start <- as.numeric(Sys.time())
  # Save list of nodes
  nodes <- V(graph)
  influence <- 0
  seed <- NULL
  while (length(seed) < budget) {
    max_influence <- 0
    most_influential <- NULL
    current <- NULL
    # For all nodes except seed
    for (node in setdiff(nodes, seed)) {
      # Find infuence of node with existing nodes in seed
      current <- get_influence(graph, c(seed, node), test_method)
      # If current node causes more influence than maximum so far, then swap
      if (current > max_influence) {
        most_influential <- node
        max_influence <- current
      }
    }
    # At the end, we should have node with maximum influence to add to influential_nodes
    seed <- c(seed, most_influential)
  }
  end <- as.numeric(Sys.time())
  output <- NULL
  output$influential_nodes <- V(graph)[seed]
  output$time <- (end - start)
  output$influence <- get_influence(graph, output$influential_nodes, test_method=test_method)
  output
}

#' @title Implements optimal algorithm for Influence Maximization
#' @name optimal_influential
#' @param graph is the igraph object
#' @param budget number of influential nodes to be fetched. Default value is 1
#' @param prob probability at which a node influences its neighbours
#' @param test_method specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"
#' @param parallel when true, executes the funtion using multiple CPU cores. Default value is TRUE
#' @return object containing: 1. Vector of influential nodes. 2. Measure of influence. 3. Elapsed time in seconds.
#' @import igraph iterpc foreach parallel
#' @export
optimal_influential <- function(graph, budget, prob=0.5, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"), parallel=TRUE) {
  start <- as.numeric(Sys.time())
  combinations <- getall(iterpc(vcount(graph), budget))
  # Add another column to store total spread
  combinations <- cbind(combinations, 0)
  if (parallel) {
    if (requireNamespace("parallel", quietly=TRUE) && requireNamespace("snow", quietly=TRUE) && requireNamespace("doSNOW", quietly=TRUE)) {
      cores <- detectCores() - 1
      cl <- snow::makeCluster(cores)
      doSNOW::registerDoSNOW(cl)
      loginfo(paste("Calculating spread under", test_method))
      # foreach requires us to define each packages and function name used within it
      spreads <- foreach (i = 1:nrow(combinations), .packages=c("igraph"), .export=c("get_influence", "simulate_ic", "simulate_lt")) %dopar% {
        nodes <- combinations[i, 1:budget]
        get_influence(graph, nodes, test_method)
      }
      combinations[,(budget + 1)] <- unlist(spreads)
      # Unregister cluster
      snow::stopCluster(cl)
    }
  }
  else {
    for (i in 1:nrow(combinations)) {
      nodes <- combinations[i,1:budget]
      # Save spread to last column
      combinations[i,(budget + 1)] <- get_influence(graph, nodes, test_method)
    }
  }
  end <- as.numeric(Sys.time())
  output <- NULL
  output$influence <- max(combinations[, (budget + 1)])
  influentials <- combinations[combinations[, (budget + 1)] == output$influence,][-(budget + 1)]
  output$influential_nodes <- V(graph)[influentials]
  output$time <- (end - start)
  output
}

#' @title Calculates influence of k nodes under Independent Cascade model
#' @name influence_ic
#' @param graph is the igraph object
#' @param seed is the initial seed nodes passed
#' @param steps is the time steps for which, the diffusion process should run. If exhaustive run is required, provide a high value (like 100). Default value is 1
#' @param prob is the probability of activation of a neighbour node. This is applicable only to IC model currently
#' @return output containing summary, including no. of nodes activated and time taken
influence_ic <- function(graph, seed, steps, prob) {
  # Algorithm: Independent Cascade model takes a network (graph) as input and some budget (k).
  # From G, k fraction of nodes are initially activated by some method. Next, we attempt to activate more nodes in the neighbourhood of these nodes.
  # Each active node attempts to activate each of its neighbour nodes with a global probability p (this is 0.5 for coin toss method)
  # Whether an attempt succeeds or fails, a node cannot be attempted for activation twice by any of the active neighbours.

  # Save the start time
  start <- as.numeric(Sys.time())
  # Read graph from file
  G <- graph
  # Save list of nodes
  nodes <- V(graph)
  # Save list of edges
  edges <- E(graph)
  influence <- 0
  output <- NULL
  output$initial_seed <- c(seed)
  attempted <- seed
  for (t in 1:steps) {
    # If all nodes have been attempted, then break
    if (length(attempted) >= length(nodes) - length(seed)) {
      break
    }
    active <- NULL
    for (v in seed) {
      # Select all neighbours of v, exempting nodes that have already been attempted
      neighbours <- setdiff(neighbors(G, v), attempted)
      if(length(neighbours) == 0) {
        next
      }
      # Store all nodes in active that had successful trial
      activated <- unlist(lapply(neighbours, function(neighbours)
        if (runif(1) >= (1 - prob)) {neighbours}))
      attempted <- unique(c(attempted, neighbours))
      active <- c(active, activated)
    }
    seed <- c(seed, active)
    #print(c("Active in step", t, "=", length(active)))
    influence <- influence + length(active)
  }
  end <- as.numeric (Sys.time())
  # Summary
  output$influence <- influence
  output$time <- (end - start)
  output$activated <- seed
  output
}

#' @title Calculates influence of k nodes under Linear Threshold model
#' @name influence_lt
#' @param graph is the igraph object
#' @param seed is the initial seed nodes passed
#' @param steps is the time steps for which, the diffusion process should run. If exhaustive run is required, provide a high value (like 100). Default value is 1
#' @param threshold is minimum threshold required to activate a node under observation
#' @return output containing summary, including no. of nodes activated and time taken
influence_lt <- function(graph, seed, steps, threshold) {
  # Algorithm: Linear Threshold model takes a network (graph) as input and some budget (k).
  # From G, k fraction of nodes are initially activated randomly. Then we attempt to activate more nodes in the neighbourhood of these nodes.
  # A node v actiates only if sum of weights of its active neighbour nodes equals or exceeds its threshold (assigned randomly here).
  # In the given function, if the fraction of active nodes in neighbourhood equals or exceeds the threshold, the inactive node becomes active
  # The process continues for t steps, in each step, the nodes activated in step t-1 also take part in diffusion process

  # Save the start time
  start <- as.numeric(Sys.time())
  # Read graph from file
  G <- graph
  # Save list of nodes
  nodes <- V(graph)
  # Save list of edges
  edges <- E(graph)
  influence <- 0
  output <- NULL
  output$initial_seed <- c(seed)
  attempted <- seed
  activated <- NULL
  for (t in 1:steps) {
    # If all nodes have been attempted, then break
    if (length(attempted) >= length(nodes) - length(seed)) {
      break
    }
    active <- NULL
    # Select all nodes having at least one neighbour in seed nodes
    inactive <- unlist(lapply(seed, function(seed) {neighbors(G, seed)}))
    # Remove nodes that have already been attempted
    inactive <- setdiff(inactive, attempted)
    # Filter duplicates
    inactive <- unique(inactive)
    for (u in inactive) {
      # Every seed node in the neighbourhood will attempt to activate u with probability p
      neighbours <- neighbors(G, u)
      active_neighbours <- intersect(neighbours, seed)
      if (length(neighbours) == 0) {
        next
      }
      # If ratio of active nodes in neighbourhood of u is greater than or equal to threshold, then activate u
      ratio <- (length(active_neighbours) / length(neighbours))
      if (ratio >= threshold) {
        active <- c(active, u)
      }
      # Active or not, this node has been attempted
      attempted <- c(attempted, u)
    }
    #print (paste("Attempted on in this step:", length(inactive), "Activated:", length(active)))
    activated <- c(activated, active)
    seed <- active
  }
  end <- as.numeric (Sys.time())
  # Summary
  output$influence <- length(activated)
  output$time <- (end - start)
  output
}

#' @title Calculates spread under IC model
#' @name ic_spread
#' @param graph is the weighted igraph object
#' @param seed is a set of seed (initial nodes)
#' @param runs is the number of times the loop should run
#' @return output average spread
#' @examples {
#' graph <- erdos.renyi.game(500, 0.005)
#' ic_spread(graph, seed=c(2,5,9,23), runs=10)
#' }
#' @import igraph
#' @export
ic_spread <- function (graph, seed, runs=100) {
  total <- 0
  for (i in 1:runs) {
    active <- NULL
    count <- 0
    # Activate seed nodes
    for (node in seed) {
      count <- count + 1
      active <- c(active, node)
    }
    count <- count + simulate_ic(graph, active);
    total <- total + count;
  }
  round(total / runs, 5)
}

#' @title Calculates spread under IC model
#' @name ic_spread_plus
#' @param graph is the igraph object
#' @param seed is a set of seed (initial nodes)
#' @param runs is the number of times the loop should run
#' @param best_node is the best known node
#' @return output average spread
#' @examples {
#' graph <- erdos.renyi.game(500, 0.005)
#' ic_spread_plus(graph, seed=c(2,5,9,23), runs=10, best_node=2)
#' }
#' @import igraph
#' @export
ic_spread_plus <- function (graph, seed, runs=100, best_node=0) {
  total <- 0
  for (i in 1:runs) {
    active <- NULL
    count <- 0
    # Activate seed nodes
    for (node in seed) {
      count <- count + 1
      active <- c(active, node)
    }
    count <- count + simulate_ic(graph, active);
    total <- total + count
    #print(paste('Spread for run #', i, count))
    # Compute next step based on previous best
    if (best_node > 0 & (!best_node %in% active)) {
      active <- c(active, best_node)
      count <- 1
      count <- count + simulate_ic(graph, active);
      total <- total + count
    }
  }
  round(total / runs, 5)
}

#' @title Simulates influence spread under Independent Cascade model
#' @name simulate_ic
#' @param graph is the weighted igraph object
#' @param active represents number of active nodes in the graph
#' @return number of nodes activated during simulation
#' @import igraph
#' @export
simulate_ic <- function(graph, active) {
  # Algorithm: given a weighted graph G and a set of active nodes V,
  # each node u in V attempts to activate its neighbours with probability equal to the weight on its edge.
  # If a coin toss with this probability is successful, then the inactive neighbour gets activated.
  # Once active, a node does not deactivate
  count <- 0
  # If the graph is unweighted, then default the weights to 1
  if (!is_weighted(graph)) {
    E(graph)$weight <- 0.5
  } else {
    E(graph)$weight <- normalize_trait(E(graph)$weight)
  }
  tried <- NULL
  for (node in active) {
    # Fetch neighbours of node
    neighbour_nodes <- neighbors(graph, node)
    # Remove already activated nodes from neighbours
    neighbour_nodes <- neighbour_nodes[!neighbour_nodes %in% active]
    # Remove already tried to be influenced from neighbours
    neighbour_nodes <- neighbour_nodes[!neighbour_nodes %in% tried]
    if (length(neighbour_nodes) == 0) {
      next
    }
    # Try to activate inactive neighbours according to the weight on edge
    for (j in 1:length(neighbour_nodes)) {
      weight <- E(graph, P=c(node, neighbour_nodes[j]))$weight
      if (runif(1) <= weight) {
        count <- count + 1
      }
      tried <- c(tried, neighbour_nodes[j])
    }
  }
  count
}

#' @title Calculates spread under LT model
#' @name lt_spread
#' @param graph is the weighted igraph object
#' @param seed is a set of seed (initial nodes)
#' @param runs is the number of times the loop should run
#' @return output average spread
#' @import igraph
#' @export
lt_spread <- function (graph, seed, runs=100) {
  total <- 0
  for (i in 1:runs) {
    active <- NULL
    count <- 0
    # Activate seed nodes
    for (node in seed) {
      count <- count + 1
      active <- c(active, node)
    }
    count <- count + simulate_lt(graph, active);
    total <- total + count;
  }
  round(total / runs, 5)
}

#' @title Simulates influence spread under Linear Threshold model
#' @name simulate_lt
#' @param graph is the weighted igraph object
#' @param active represents number of active nodes in the graph
#' @param threshold is the linear threshold between 0 and 1 with which a node is influenced
#' @return number of nodes activated during simulation
#' @import igraph
#' @export
simulate_lt <- function(graph, active, threshold=0.5) {
  # Algorithm: given a weighted graph G and a set of active nodes V,
  # each inactive node in the graph gets a chance to be activated with the probability being the collective weights on its edges with active nodes.
  # If a coin toss with probability as the sum of weights of active neighbours is greater than given threshold, then the inactive node gets activated.
  # Once active, a node does not deactivate
  count <- 0
  # If the graph is unweighted, then default the weights to 1
  if (!is_weighted(graph)) {
    E(graph)$weight <- 0.5
  }
  inactive <- unlist(lapply(active, function(active) {neighbors(graph, active)}))
  inactive <- setdiff(inactive, active)
  for (u in inactive) {
    neighbours <- neighbors(graph, u)
    active_neighbours <- intersect(neighbours, active)
    if (length(neighbours) == 0) {
      next
    }
    # If ratio of active nodes in neighbourhood of u is greater than or equal to threshold, then activate u
    ratio <- (length(active_neighbours) / length(neighbours))
    if (ratio >= threshold) {
      count <- count + 1
    }
    active <- c(active, u)
  }
  count
}

#' @title Returns the centrality values of nodes in a graph using given method
#' @name get_centrality_scores
#' @param graph is the igraph object
#' @param centrality_method defines the centrality method to be used. Value must be:
#' @return vector of centrality scores
#' @import igraph
#' @examples {
#' graph <- erdos.renyi.game(500, 0.05)
#' get_centrality_scores(graph, centrality_method="DEGREE")
#' get_centrality_scores(graph, centrality_method="BETWEENNESS")
#' get_centrality_scores(graph, centrality_method="CLOSENESS")
#' get_centrality_scores(graph, centrality_method="EIGENVECTOR")
#' }
#' @export
get_centrality_scores <- function(graph, centrality_method=c("DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR")) {
  if (centrality_method == "DEGREE") {
    # Calculate in/out degrees of all nodes
    degree(graph, V(graph), mode="all", loops=FALSE, normalized=FALSE)
  }
  else if(centrality_method == "BETWEENNESS") {
    # Calculate betweenness centrality (for huge data sets, use betweenness.estimate() and give some max value of path length as cutoff)
    betweenness(graph, V(graph), directed=FALSE)
  }
  else if (centrality_method == "CLOSENESS") {
    # Calculate in/out closeness of all nodes, normalized between 0 and 1
    closeness(graph, V(graph), mode="all", normalized=TRUE)
  }
  else if (centrality_method == "EIGENVECTOR") {
    # Calculate eigenvectors of the graph
    eigen <- evcent(graph, directed=FALSE)
    eigen$vector
  }
}

#' @title Calculates resilience of network
#' @name resilience
#' @param graph is the weighted igraph object
#' @param nodes is a set of nodes to check resilience of
#' @return number of remaining nodes in largest connected component after removing nodes
#' @examples {
#' graph <- erdos.renyi.game(500, 0.005)
#' resilience(graph, nodes=V(graph)[c(2,5,9,23)])
#' }
#' @import igraph
#' @export
resilience <- function (graph, nodes) {
  graph <- delete.vertices(graph, nodes)
  graph <- largest_component(graph)
  vcount(graph)
}

#' @title Quantifies the influence of a set of nodes in a graph
#' @name get_influence
#' @param graph is the igraph object
#' @param nodes the set of nodes to calculate influence for
#' @param test_method specifies the method to measure influence. Value "RESILIENCE" (number of total nodes REMOVED (NOT THE REMAINING ones as in original resilience function) from the graph); "INFLUENCE_IC" (see simulate_ic method); "INFLUENCE_LT" (see simulate_lt method). Default is "RESILIENCE"
#' @return vector of resiliences of provided combinations
#' @examples {
#' graph <- erdos.renyi.game(100, 0.2)
#' get_influence(graph=graph, nodes=V(graph)[1:10], test_method="RESILIENCE")
#' }
#' @import igraph
#' @export
get_influence <- function(graph, nodes, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC")) {
  if (test_method == "RESILIENCE") {
    vcount(graph) - resilience(graph, nodes)
  } else if (test_method == "INFLUENCE_IC") {
    simulate_ic(graph, nodes)
  } else if (test_method == "INFLUENCE_LT") {
    simulate_lt(graph, nodes)
  }
}

#' @title Returns CI value of given graph and node
#' @name collective_influence
#' @param graph the igraph object
#' @param neighborhood_distance is the distance to which the neighborhood nodes are searched for
#' @param node_id is the ID of the target node
#' @param method is the metric to calculate sum of influence. Default is "degree"
#' @return influence as product of degree of target node and total sum of degrees of neighborhood
#' TODO: extend the function and include adaptive methods as well as other centrality methods
#' @examples {collective_influence(graph=erdos.renyi.game(100, 0.2), neighborhood_distance=2, 1)}
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
