#' @title Finds communities in the given graph and returns the graph after adding a vector "group" to its vertices
#' @name find_communities
#' @param graph is the igraph object
#' @param plot whether to plot the \code{graph} using \link{tkplot}. Default is TRUE
#' @param method is the method to find communities. Value can be "multilevel", "edgebetweenness", "fastgreedy", "eigenvector", "spinglass", "walktrap", "labelpropagation", "clique", "largescale"
#' @return \code{graph} object with membership attribute attached representing communities
#' @import igraph graphics
#' @importFrom grDevices rainbow
#' @export
find_communities <- function(graph, plot=TRUE, method=c("multilevel", "edgebetweenness", "fastgreedy", "eigenvector", "spinglass", "walktrap", "labelpropagation", "clique", "largescale")) {
  # multilevel.community: based on Louvaine's algorithm, it is better at scaling and avoids formation of super communities. In this method, instead of merging communities, nodes are moved between communities such that each node makes a local decision that maximizes its own contribution to the modularity score.
  # When this procedure gets stuck (i.e. none of the nodes change their membership), then all the communities are collapsed into single nodes and the process continues (thus the name multilevel).
  if (method == "multilevel") {
    communities <- multilevel.community(graph)
    V(graph)$group <- communities$membership
  }
  # edge.betweenness.community: a Hierarchical decomposition process where edges are removed in the decreasing order of their edge betweenness scores. This is motivated by the fact that edges connecting different groups are more likely to be contained in multiple shortest paths simply because in many cases they are the only option to go from one group to another.
  # This method yields good results but is very slow because of the computational complexity of edge betweenness calculations and because the betweenness scores have to be re-calculated after every edge removal.
  # Another disadvantage is that it builds a full dendrogram and does not tell where to cut the dendrogram to obtain the final groups, so we use some other measure (e.g. modularity score of the partitions) to decide that
  else if (method == "edgebetweenness") {
    communities <- edge.betweenness.community(graph)
    V(graph)$group <- communities$membership
  }
  # leading.eigenvector.community: a top-down hierarchical approach that optimizes the modularity function. In each step, the graph is split into two parts in a way that the separation itself yields a significant increase in the modularity. The split is determined by evaluating the leading eigenvector of modularity matrix, and there is also a stopping condition which prevents tightly connected groups to be split further.
  # Due to the eigenvector calculations involved, it might not work on degenerate graphs where the ARPACK eigenvector solver is unstable. On non-degenerate graphs, it is likely to yield a higher modularity score than the fast greedy method, although it is a bit slower.
  else if (method == "eigenvector") {
    communities <- leading.eigenvector.community(graph)
    V(graph)$group <- communities$membership
  }
  # fastgreedy.community: a bottom-up hierarchical approach. It tries to optimize function modularity function in greedy manner. Initially, every vertex belongs to a separate community, and communities are merged iteratively such that each merge is locally optimal (i.e. has high increase in modularity value).
  # The algorithm stops when it is not possible to increase the modularity any more, so it gives you a grouping as well as a dendrogram. The method is fast and it is the method that is usually tried as a first approximation because it has no parameters to tune.
  # However, it has a limitation that communities below a given size threshold will always be merged with neighboring communities
  else if (method == "fastgreedy") {
    communities <- fastgreedy.community(graph)
    V(graph)$group <- membership(communities)
  }
  # spinglass.community: an approach based on Potts model. Each node can be in one of 'c' spin states, and the edges specify which pairs of vertices would prefer to stay in the same spin state and which ones prefer to have different spin states.
  # The model is then simulated for a given number of steps, and the spin states of the nodes in the end define the communities.
  # The consequences are that 1) There will never be more than 'c' communities in the end, although you can set c to as high as 200; 2) There may be less than 'c' communities in the end as some of the spin states may become empty; 3) In disconnected networks, it is not guaranteed that nodes in disconencted parts of the networks have different spin states.
  # The method is not particularly fast and not deterministic, but has a tunable resolution parameter that determines the cluster sizes.
  else if (method == "spinglass") {
    communities <- spinglass.community(graph, spins=10)
    V(graph)$group <- communities$membership
  }
  # walktrap.community: the general idea is that if you perform random walks on the graph, then the walks are more likely to stay within the same community because there are only a few edges that lead outside a given community.
  # Walktrap runs short random walks of 3-4-5 steps (depending on parameters) and uses the results of these random walks to merge separate communities in a bottom-up manner like fastgreedy.community. We can use the modularity score to select where to cut the dendrogram.
  # It is a bit slower than the fast greedy approach but also a bit more accurate.
  else if (method == "walktrap") {
    communities <- walktrap.community(graph)
    V(graph)$group <- membership(communities)
  }
  # label.propagation.community: a simple approach in which every node is assigned one of 'k' labels. The method then proceeds iteratively and re-assigns labels to nodes in a way that each node takes the most frequent label of its neighbors in a synchronous manner. The method stops when the label of each node is one of the most frequent labels in its neighborhood.
  # It is very fast but yields different results based on the initial configuration (which is decided randomly), therefore it should be run a large number of times before labeling.
  else if (method == "labelpropagation") {
    V(graph)$group <- label.propagation.community(graph)$membership
  }
  else if (method == "clique") {
    # TODO
  }
  else if (method == "largescale") {
    large_scale_community(graph)
  }
  # Plot the graph, showing communities
  if (plot) {
    graph$layout <- layout.kamada.kawai
    size <- length(unique(V(graph)$group))
    V(graph)$color <- rainbow(size)[V(graph)$group]
    plot(graph)
  }
  graph
}

#' @title Detects communities using cliques
#' @name clique_community
#' @param graph is igraph object
#' @param k is the number of communitites to find
#' @references Palla, Gergely, et al. "Uncovering the overlapping community structure of complex networks in nature and society." Nature 435.7043 (2005): 814-818.
#' @import igraph
#' @export
clique_community <- function(graph, k) {
  clq <- cliques(graph, min=k, max=k)
  edges <- c()
  for (i in seq_along(clq)) {
    for (j in seq_along(clq)) {
      if ( length(unique(c(clq[[i]], clq[[j]]))) == k+1 ) {
        edges <- c(edges, c(i,j))
      }
    }
  }
  clq.graph <- simplify(graph(edges))
  V(clq.graph)$name <- seq_len(vcount(clq.graph))
  comps <- decompose.graph(clq.graph)
  lapply(comps, function(x) { unique(unlist(clq[V(x)$name])) })
}

#' @title Detects communities in large-scale graphs
#' @name large_scale_community
#' @param graph is igraph object
#' @param mode is detection mode. Values can be "all", "in" or "out". Default is "all"
#' @references Raghavan, Usha Nandini, Réka Albert, and Soundar Kumara. "Near linear time algorithm to detect community structures in large-scale networks." Physical Review E 76.3 (2007): 036106.
#' @import igraph
#' @export
large_scale_community <- function(graph, mode="all") {
  group <- as.numeric(V(graph))
  current_order <- sample(vcount(graph), vcount(graph))
  t <- 0
  done <- FALSE
  while(!done){
    t <- t + 1
    done <- TRUE ## change to FALSE whenever a node changes groups
    for(v in current_order){
      ## get the neighbor group frequencies:
      freq <- table(group[V(graph)[neighbors(graph, v, mode=mode)]])
      ## pick one of the most frequent:
      new_group <- sample(names(freq) [freq==max(freq)],1)
      if(done) {
        done <- new_group == group[V(graph)[v]]
      }
      group[V(graph)[v]] <- new_group
    }
  }
  # V(graph)$group <- group
  ## now fix any distinct groups with same labels:
  # for(i in unique(V(graph)$group)) {
  #   ## only bother for connected groups
  #   if(!is.connected(induced_subgraph(graph, V(graph)[group==i]))) {
  #     theseNodes <- V(graph)[group==i]
  #     theseClusters <- clusters(subgraph(graph, theseNodes))
  #     ## iterate through the clusters and append their names
  #     for(j in unique(theseClusters$membership)) {
  #       V(graph)[theseNodes[theseClusters$membership==j]]$group <- paste(i,j,sep=".")
  #     }
  #   }
  # }
  graph
}

#' @title Detects function performs a Wilcoxon rank-sum test on the "internal" and "external" degrees of a community in order to quantify its significance.
#' @name community_significance_test
#' @description The edges within a community are "internal" and the edges connecting the vertices of a community with the rest of the graph are "external". More internal than external edges show that the community is significant; the otherwise suggests that the community is in fact an "anti-community".
#' @param graph is igraph object
#' @param vs is the vector of vertices of the original graph which will form the subgraph
#' @param ... arguments to pass to \link{wilcox.test} function
#' @import igraph stats
#' @export
community_significance_test <- function(graph, vs, ...) {
  if (is.directed(graph)) {
    stop("This method requires an undirected graph")
  }
  subgraph <- induced.subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}

#' @title Detects function performs modularity test of communities.
#' @name community_modularity_test
#' @param graph is igraph object
#' @import igraph
community_modularity_test <- function(graph) {
  # TODO:
}

