#' @title Gives a summary of common metrics of given graph
#' @name graph_summary
#' @param graph is the igraph object
#' @param plot uses tkplot to plot the \code{graph}. Default is FALSE
#' @return object containing summary
#' @examples {
#' graph_summary(erdos.renyi.game(20, 0.2))
#' }
#' @import igraph
#' @export
graph_summary <- function(graph, plot=FALSE) {
  o <- NULL
  degrees <- degree(graph)
  o$edges <- ecount(graph)
  o$vertices <- vcount(graph)
  o$vertex_edge_ratio <- o$vertices / o$edges
  o$connected <- is.connected(graph)
  o$average_degree <- mean(degrees)
  o$average_path_length <- average.path.length(graph)
  o$highest_degree <- max(degrees)
  o$density <- graph.density(graph)
  o$diameter <- diameter(graph)
  o$transitivity <- transitivity(graph)
  o$assortativity <- assortativity.degree(graph)
  o$average_distance <- mean_distance(graph)
  o$graph_triads <- length(triangles(graph))
  o$girth <- girth(graph)$girth
  o$power_law <- fit_power_law(graph)$alpha
  if (plot) {
    hist(degree(graph))
    tkplot(graph)
  }
  o
}

#' Calculates several traits from given graph and returns as data frame
#'
#' @name get_graph_traits
#' @param graph is the igraph object
#' @param normalize uses pnorm function to normalize the traits. Default is FALSE
#' @param graph_traits is the vector of several graph metrices to calculate
#' @param node_traits is the vector of several node metrices to calculate
#' @param verbose prints the progress log on console when TRUE. Default is FALSE
#' @return data frame containing graph and its traits
#' @import igraph
#' @export
get_graph_traits <- function(graph, normalize=FALSE,
    graph_traits=c("SIZE", "EDGES", "AVERAGE_DEGREE", "MAX_DEGREE", "AVERAGE_PATH_LENGTH", "CLUSTERING_COEFFICIENT", "DIAMETER", "DENSITY", "ASSORTATIVITY", "AVERAGE_DISTANCE", "TRIADS", "GIRTH"),
    node_traits=c("DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR", "ECCENTRICITY", "CORENESS", "PAGERANK", "COLLECTIVE_INFLUENCE", "ADAPTIVE_DEGREE", "ADAPTIVE_BETWEENNESS", "ADAPTIVE_CLOSENESS", "ADAPTIVE_EIGENVALUE", "ADAPTIVE_ECCENTRICITY", "ADAPTIVE_CORENESS", "ADAPTIVE_PAGERANK", "ADAPTIVE_COLLECTIVE_INFLUENCE"), verbose=FALSE) {
  # First, fetch all the node traits
  if (verbose) {
    print("Computing centrality traits...")
  }
  data <- NULL
  data$name <- 1:vcount(graph) - 1
  data$degree <- get_centrality_scores(graph, "DEGREE", normalize=normalize)
  if ("BETWEENNESS" %in% node_traits)
    data$betweenness <- get_centrality_scores(graph, "BETWEENNESS", normalize=normalize)
  if ("CLOSENESS" %in% node_traits)
    data$closeness <- get_centrality_scores(graph, "CLOSENESS", normalize=normalize)
  if ("EIGENVECTOR" %in% node_traits)
    data$eigenvalue <- get_centrality_scores(graph, "EIGENVECTOR", normalize=normalize)
  if ("ECCENTRICITY" %in% node_traits)
    data$eccentricity <- get_centrality_scores(graph, "ECCENTRICITY", normalize=normalize)
  if (verbose) {
    print("Computing node heuristic traits...")
  }
  if ("CORENESS" %in% node_traits) {
    data$coreness <- coreness(graph)
    if (normalize) {
      data$coreness <- normalize_trait(data$coreness)
    }
  }
  if ("PAGERANK" %in% node_traits) {
    data$pagerank <- page_rank(graph)$vector
    if (normalize) {
      data$pagerank <- normalize_trait(data$pagerank)
    }
  }
  if ("COLLECTIVE_INFLUENCE" %in% node_traits) {
    data$ci <- sapply(V(graph), function(x) { collective_influence(graph, neighborhood_distance=2, x) })
  }
  if (verbose) {
    print("Computing graph traits...")
  }
  if ("SIZE" %in% graph_traits)
    data$graph_size <- vcount(graph)
  if ("EDGES" %in% graph_traits)
    data$graph_edges <- ecount(graph)
  if ("AVERAGE_DEGREE" %in% graph_traits)
    data$graph_avg_degree <- mean(data$degree)
  if ("MAX_DEGREE" %in% graph_traits)
    data$graph_max_degree <- max(data$degree)
  if ("AVERAGE_PATH_LENGTH" %in% graph_traits)
    data$graph_apl <- average.path.length(graph)
  if ("CLUSTERING_COEFFICIENT" %in% graph_traits)
    data$graph_clust_coef <- transitivity(graph)
  if ("DIAMETER" %in% graph_traits)
    data$graph_diameter <- diameter(graph)
  if ("DENSITY" %in% graph_traits)
    data$graph_density <- graph.density(graph)
  if ("ASSORTATIVITY" %in% graph_traits)
    data$graph_assortativity <- assortativity.degree(graph)
  if ("AVERAGE_DISTANCE" %in% graph_traits)
    data$graph_avg_distance <- mean_distance(graph)
  if ("TRIADS" %in% graph_traits)
    data$graph_triads <- length(triangles(graph))
  if ("GIRTH" %in% graph_traits)
    data$graph_girth <- girth(graph)$girth
  data
}

#' @title Normalizes the numeric values passed in \code{x} between 0 and 1
#' @name normalize_trait
#' @param x is data to be normalized
#' @return normalized \code{x}
#' @import igraph
#' @export
normalize_trait <- function(x) {
  stats::pnorm(x, mean(x), sd(x))
}

#' @title  Normalize a set of numeric variables in a dataset between 0 and 1. Non-numeric data will be skipped
#' @name normalize_data
#' @param data is data frame to be normalized
#' @param columns is list of columns to be normalized
#' @return normalized data frame
#' @import igraph
#' @export
normalize_data <- function(data, columns) {
  for (column in columns) {
    # Skip non-numeric data
    if (mode(data[, column]) != "numeric") {
      next
    }
    x <- data[, column]
    x <- stats::pnorm(x, mean(x), sd(x))
    data[, column] <- x
  }
  data
}

#' @title Plots degree distribution of given graph
#' @name plot_degree_distribution
#' @param graph is the igraph object
#' @import igraph
#' @export
plot_degree_distribution <- function(graph) {
  degree = degree(graph, mode="all")
  distribution = degree.distribution(graph, mode="all", cumulative=FALSE)
  degree = 1:max(degree)
  probability = distribution[-1]
  # Remove blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # plot
  plot(probability ~ degree, log="xy", xlab="Degree (log)", ylab="Probability (log)", col=1, main="Degree Distribution")
}

#' @title  Plots degree distribution and returns power-law exponent of given graph
#' @name fit_power_law
#' @param graph is the igraph object
#' @import igraph
#' @export
fit_power_law = function(graph) {
  distribution = degree.distribution(graph, mode="all", cumulative=FALSE)
  degree = 1:max(degree(graph, mode="all"))
  probability = distribution[-1]
  # Remove blank values
  nonzero = which(probability != 0)
  probability = probability[nonzero]
  degree = degree[nonzero]
  # plot
  plot(probability ~ degree, log="xy", xlab="Degree (log)", ylab="Probability (log)", col=1, main="Degree Distribution")
  # Return alpha, the exponent of fitted power-law
  igraph::fit_power_law(degree(graph))[2]
}

#' @title Generates a tree-structured graph
#' @name generate_tree
#' @param size is the number of nodes
#' @param children is the number of children each node has (in addition to a parent node)
#' @param direction defines whether the edges are directed inwards, outwards or undirected. Possible values can be 'in', 'out' and 'undirected' (default)
#' @return igraph object
#' @import igraph
#' @export
generate_tree <- function(size, children=2, direction='undirected') {
  graph.tree(size, children, mode=direction)
}

#' @title Generates a ring-structured graph, in which nodes are connected with neighbours within given distance
#' @name generate_ring
#' @param size is the number of nodes
#' @param distance defines maximum distance each node has to its farthest direct neighbour
#' @return igraph object
#' @import igraph
#' @export
generate_ring <- function(size, distance) {
  connect.neighborhood(graph.ring(size), distance)
}

#' @title Generates a fully connected undirected graph
#' @name generate_clique
#' @param size is the number of nodes
#' @return igraph object
#' @import igraph
#' @export
generate_clique <- function(size) {
  graph.full(size)
}

#' @title Generates a Erdos Renyi random graph
#' @name generate_random
#' @param size is the number of nodes
#' @param probability is the probability of edge formation between nodes
#' @param directed generates directed graph when TRUE. Default value is FALSE
#' @param allow_cycles produces loops in the graph when TRUE. Default value is FALSE
#' @return igraph object
#' @import igraph
#' @export
generate_random <- function(size, probability=0.2, directed=FALSE, allow_cycles=FALSE) {
  erdos.renyi.game(size, probability, directed=directed, loops=allow_cycles)
}

#' @title Generates a Watts & Strogatz small-world graph by rewiring a random graph, while keeping the degree distribution consistent
#' @name generate_small_world
#' @param size is the number of nodes
#' @param probability is the probability of edge formation between nodes
#' @param directed generates directed graph when TRUE. Default value is FALSE
#' @param allow_cycles produces loops in the graph when TRUE. Default value is FALSE
#' @return igraph object
#' @export
generate_small_world <- function(size, probability=0.1, directed=FALSE, allow_cycles=FALSE) {
  graph <- generate_random(size, probability, directed, allow_cycles)
  iterations <- size * 10
  rewire(graph, with=keeping_degseq(allow_cycles, niter=iterations))
}

#' @title Generates a Barabasi scale-free graph
#' @name generate_scale_free
#' @param size is the number of nodes
#' @param preference is the power of preferencial attachment. Default is linear, i.e. 1
#' @param directed generates directed graph when TRUE. Default value is FALSE
#' @param allow_cycles produces loops in the graph when TRUE. Default value is FALSE
#' @return igraph object
#' @export
generate_scale_free <- function(size, preference=1, directed=FALSE, allow_cycles=FALSE) {
  barabasi.game(size, power=preference, directed=directed)
}

#' @title Generates a Holme-Kim Network
#' @name generate_holme_kim
#' @description Simulate a scale-free network with relatively high clustering, comparing to B-A networks (Holme and Kim, 1999).
#' @param size is the number of nodes of the network
#' @param m is the number of nodes to which a new node connects at each iteration
#' @param triad_prob is Triad formation probability after each preferential attachment mechanism
#' @param directed whether the graph is directed or not. Default is FALSE
#' @details The Holme-Kim network model is a simple extension of B-A model. It adds an additional step, called "Triad formation", with the probability \emph{pt} that compensates the low clustering in B-A networks.
#' @return A list containing the nodes of the network and their respective neighbors.
#' @author Xu Dong, Nazrul Shaikh
#' @examples {generate_holme_kim (1000, 20, 0.1)}
#' @references Holme, Petter, and Beom Jun Kim. "Growing scale-free networks with tunable clustering."Physical review E65, no. 2 (2002): 026107.
#' @import igraph
#' @export
generate_holme_kim <- function(size, m, triad_prob=0.1, directed=FALSE) {
  if (size < 0 | size %% 1 != 0) stop("Parameter 'n' must be positive integer", call. = FALSE)
  if (m < 1 | m %% 1 != 0) stop("Parameter 'm' must be integer  greater than 1", call. = FALSE)
  if (triad_prob < 0 | triad_prob > 1) stop("Parameter 'pt' must be in (0,1)", call. = FALSE)
  graph <- list()
  graph[size] <- list(NULL)
  ## Create the m0 graph + (m + 1) node
  graph[[m + 1]] <- seq(m)
  for ( k in seq(m)) {
    graph[[k]] <- m + 1
  }
  df <- c(rep(1, m), m, rep(0, size - m - 1))
  for (i in (m + 2):size){
    pa.neighbor <- sample(seq(size), 1, prob=df)
    graph[[i]] <- pa.neighbor
    graph[[pa.neighbor]] <- c(graph[[pa.neighbor]], i)
    df[pa.neighbor] <- df[pa.neighbor] + 1
    for (j in seq(2, m)) {
      pool <- setdiff(graph[[pa.neighbor]], c(i, graph[[i]]))
      if (stats::runif(1) <= triad_prob && length(pool) != 0) {
        tf.neighbor <- sample(pool, 1)
        graph[[i]] <- c(graph[[i]], tf.neighbor)
        graph[[tf.neighbor]] = c(graph[[tf.neighbor]], i)
        df[tf.neighbor] <- df[tf.neighbor] + 1
      } else {
        pa.neighbor <- sample(seq(size)[-graph[[i]]], 1, prob=df[-graph[[i]]])
        graph[[i]] <- c(graph[[i]], pa.neighbor)
        graph[[pa.neighbor]] <- c(graph[[pa.neighbor]], i)
        df[pa.neighbor] <- df[pa.neighbor] + 1
      }
    }
    df[i] <- m
  }
  mode <- ifelse(directed, "in", "total")
  graph.adjlist(graph, mode="total")
}

#' @title Returns largest connected component in a network
#' @name largest_component
#' @param graph is the igraph object
#' @return largest component igraph object
#' @import igraph
#' @export
largest_component <- function(graph) {
  gclust = igraph::clusters(graph)
  lcc = induced.subgraph(graph, V(graph)[which(gclust$membership == which.max(gclust$csize))])
  lcc
}

#' @title Returns the centrality values of nodes in a graph using given method
#' @name get_centrality_scores
#' @param graph is the igraph object
#' @param centrality_method defines the centrality method to be used. Value must be: "DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR", "ECCENTRICITY"
#' @param normalize scales the values in the output vector between 0 and 1
#' @return vector of centrality scores
#' @import igraph
#' @examples {
#' graph <- erdos.renyi.game(500, 0.05)
#' get_centrality_scores(graph, centrality_method="DEGREE")
#' get_centrality_scores(graph, centrality_method="BETWEENNESS")
#' get_centrality_scores(graph, centrality_method="CLOSENESS")
#' get_centrality_scores(graph, centrality_method="EIGENVECTOR")
#' get_centrality_scores(graph, centrality_method="ECCENTRICITY")
#' }
#' @export
get_centrality_scores <- function(graph, centrality_method=c("DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR", "ECCENTRICITY"), normalize=FALSE) {
  x <- NULL
  if (centrality_method == "DEGREE") {
    # Calculate in/out degrees of all nodes
    x <- degree(graph, V(graph), mode="all", loops=FALSE, normalized=FALSE)
  }
  else if(centrality_method == "BETWEENNESS") {
    # Calculate betweenness centrality (TODO: for huge data sets, use betweenness.estimate() and give some max value of path length as cutoff)
    x <- betweenness(graph, V(graph), directed=FALSE)
  }
  else if (centrality_method == "CLOSENESS") {
    # Calculate in/out closeness of all nodes
    x <- closeness(graph, V(graph), mode="all", normalized=FALSE)
  }
  else if (centrality_method == "EIGENVECTOR") {
    # Calculate eigenvectors of the graph
    eigen <- evcent(graph, directed=FALSE)
    x <- eigen$vector
  }
  else if (centrality_method == "ECCENTRICITY") {
    x <- eccentricity(graph, mode="all")
  }
  if (normalize) {
    normalize_trait(x)
  }
  x
}

#' Returns ranks from 1 to highest rank before the graph is discontinued, using scores of nodes (e.g. Degree, Pagerank, Coreness, etc.)
#' @name get_adaptive_rank
#' @param graph the igraph object
#' @param ranking_method the adaptive method to use. Value must be "DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR", "ECCENTRICITY", "CORENESS", "PAGERANK", "COLLECTIVE_INFLUENCE"
#' @return vector of ranks
#' @import igraph
#' @export
get_adaptive_rank <- function(graph, ranking_method=c("DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR", "ECCENTRICITY", "CORENESS", "PAGERANK", "COLLECTIVE_INFLUENCE")) {
  V(g)$name <- V(g)
  V(g)$rank <- -1
  current_rank <- 1
  graph <- g
  while (TRUE) {
    graph <- largest_component(graph)
    param <- 0
    if (ranking_method %in% c("DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR", "ECCENTRICITY")) {
      param <- get_centrality_scores(graph, centrality_method=ranking_method)
    } else if (ranking_method == "CORENESS") {
      param <- graph.coreness(graph, mode="all")
    } else if (ranking_method == "PAGERANK") {
      param <- page_rank(graph, directed=TRUE)$vector
    } else if (ranking_method == "COLLECTIVE_INFLUENCE") {
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
