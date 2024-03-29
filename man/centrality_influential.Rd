% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/influence_maximization.R
\name{centrality_influential}
\alias{centrality_influential}
\title{Returns the most influential nodes in a graph using centrality-based heuristics}
\usage{
centrality_influential(
  graph,
  budget = 1,
  test_method = c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"),
  centrality_method = c("DEGREE", "ECCENTRICITY", "AVERAGE_DISTANCE", "BARYCENTER",
    "BETWEENNESS", "BOTTLENECK", "CENTROID", "CLOSENESS", "CLUSTERRANK",
    "COMMUNITY_BETWEENNESS", "COMMUNITY_CENTRALITY", "CROSS_CLIQUE",
    "CURRENTFLOW_CLOSENESS", "DECAY", "EDGE_PERCOLATION", "EIGENVECTOR", "ENTROPY",
    "FREEMAN_CLOSENESS", "GEODESIC_K_PATH", "HUBBELL", "KATZ", "LAPLACIAN",
    "LATORA_CLOSENESS", "LEADERRANK", "LEVERAGE", "LINCENT", "LOBBY", "MARKOV",
    "MAX_NEIGHBORHOOD_COMPONENT", "MAX_NEIGHBORHOOD_DENSITY", "PAIRWISE_DISCONNECTIVITY",
         "RADIALITY", "RESIDUAL_CLOSENESS", "SALSA", "SEMILOCAL",
    "TOPOLOGICAL_COEFFICIENT", "VITALITY_CLOSENESS")
)
}
\arguments{
\item{graph}{is the igraph object}

\item{budget}{number of influential nodes to be fetched. Default value is 1}

\item{test_method}{specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"}

\item{centrality_method}{defines the centrality method to be used. Value must be:}
}
\value{
object containing: 1. Vector of influential nodes. 2. Measure of influence. 3. Elapsed time in seconds.
}
\description{
Returns the most influential nodes in a graph using centrality-based heuristics
}
\references{
{
Harary, F., Norman, R. Z., & Cartwright, D. (1965). Structural models: An introduction to the theory of directed graphs. Wiley.;
Freeman, L. C. (1977). A set of measures of centrality based on betweenness. Sociometry, 35-41.;
Freeman, L. C. (1978). Centrality in social networks conceptual clarification. Social networks, 1(3), 215-239.
}
}
