#' @title A collaboration graph crawled from arXiv.org, High Energy Physics – Theory, from year 1991 to year 2003
#' @name load_arxiv_collaboration_graph
#' @return igraph object
#' @import igraph
#' @export
load_arxiv_collaboration_graph <- function() {
  largest_component(read.graph("dataset/arxiv_collaboration", directed=FALSE))
}

#' @title The undirected network of autonomous systems of the Internet connected with each other from the CAIDA project, collected in 2007
#' @name load_as_caida_graph
#' @return igraph object
#' @import igraph
#' @export
load_as_caida_graph <- function() {
  largest_component(read.graph("dataset/as_caida", directed=FALSE))
}

#' @title Coauthorship network of scientists working on network theory and experiment, as compiled by M. Newman in May, 2006
#' @name load_author_netscience_graph
#' @return igraph object
#' @import igraph
#' @export
load_author_netscience_graph <- function() {
  largest_component(read.graph("dataset/author_netscience", directed=FALSE))
}

#' @title Citation network of papers on influence mining from 1990 to 2015
#' @name load_influence_citation_graph
#' @return igraph object
#' @import igraph
#' @export
load_influence_citation_graph <- function() {
  largest_component(read.graph("dataset/influence_citation_network", directed=FALSE))
}

#' @title Network of American football games between Division IA colleges during regular season Fall 2000, as compiled by M. Girvan and M. Newman
#' @name load_football_graph
#' @return igraph object
#' @import igraph
#' @export
load_football_graph <- function() {
  largest_component(read.graph("dataset/football", directed=FALSE))
}

#' @title Air Transport Network is an undirected graph where nodes represent airports and edges represent a direct flight from one airport to the other
#' @name load_ita2000_graph
#' @return igraph object
#' @import igraph
#' @export
load_ita2000_graph <- function() {
  largest_component(read.graph("dataset/ita2000", directed=FALSE))
}

#' @title Class dependency network of the Java Development Kit 1.6.0.7 framework
#' @name load_jdk6_dependencies_graph
#' @return igraph object
#' @import igraph
#' @export
load_jdk6_dependencies_graph <- function() {
  largest_component(read.graph("dataset/jdk6_dependencies", directed=FALSE))
}

#' @title A social network of members of a karate club at an American University. The network represents ties between members of the club based on their social interactions both within and away from the club
#' @name load_karate_club_graph
#' @return igraph object
#' @import igraph
#' @export
load_karate_club_graph <- function() {
  largest_component(read.graph("dataset/karate_club", directed=FALSE))
}

#' @title Network of twitter links collected via snowball sampling, started from root account seekme_94
#' @name load_twitter_graph
#' @return igraph object
#' @import igraph
#' @export
load_twitter_graph <- function() {
  largest_component(read.graph("dataset/twitter_network", directed=FALSE))
}

#' @title Neural network of the nematode C. Elegans compiled by Duncan Watts and Steven Strogatz from original experimental data by White et al.
#' @name load_nematode_neural_network_graph
#' @return igraph object
#' @import igraph
#' @export
load_nematode_neural_network_graph <- function() {
  largest_component(read.graph("dataset/nematode_neural_network", directed=FALSE))
}

#' @title A directed network of hyperlinks between weblogs on US politics, recorded in 2005 by Adamic and Glance
#' @name load_political_blog_graph
#' @return igraph object
#' @import igraph
#' @export
load_political_blog_graph <- function() {
  largest_component(read.graph("dataset/political_blog", directed=TRUE))
}

#' @title Protein interaction network data for yeast (by Barabasi)
#' @name load_protein_barabasi_graph
#' @return igraph object
#' @import igraph
#' @export
load_protein_barabasi_graph <- function() {
  largest_component(read.graph("dataset/protein_barabasi", directed=FALSE))
}

#' @title Lexical network of words from the WordNet dataset. Nodes in the network are English words, and links are relationships between them, such as synonymy, antonymy, meronymy, etc.
#' @name load_wordnet_graph
#' @return igraph object
#' @import igraph
#' @export
load_wordnet_graph <- function() {
  largest_component(read.graph("dataset/wordnet", directed=FALSE))
}

#' @title Network of trade between manufacturers from around the world
#' @name load_world_trade_graph
#' @return igraph object
#' @import igraph
#' @export
load_world_trade_graph <- function() {
  largest_component(read.graph("dataset/world_trade", directed=FALSE))
}
