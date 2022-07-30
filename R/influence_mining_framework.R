# This framework provides flexible parameters to tune for influence maximization problem

library(igraph)
library(dplyr)
library(lubridate)
source('R/influence_maximization.R')
source('R/data_util.R')


#' Influence mining framework implementation based on CELF greedy method
#'
#' @name celf_greedy_framework
#' @param graph is the igraph object
#' @param criteria is the stopping criteria for influence framework. Supported values are TIME, BUDGET and SPREAD
#' @param threshold is the limit value for criteria function. The execution steps as soon as criteria meets the threshold value. If the criteria is TIME, then threshold should be defined in seconds
#' @param test_method specifies the method to measure influence. Value "RESILIENCE" (number of total nodes REMOVED (NOT THE REMAINING ones as in original resilience function) from the graph); "INFLUENCE_IC" (see simulate_ic method); "INFLUENCE_LT" (see simulate_lt method). Default is "RESILIENCE"
#' @return output object containing influential nodes, time taken and the influence of nodes according the test_method
#' @import igraph
#' @export
celf_greedy_framework <- function(graph, criteria=c("TIME","BUDGET","INFLUENCE"), threshold, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC")) {
  seed <- c()
  # For each node, compute its influence independently and add maintain in a table, while setting the flag=0
  V(graph)$name <- 1:vcount(graph)
  df <- data.frame(node=V(graph)$name, gain=0, flag=0)
  # Save the start time
  start <- as.numeric(Sys.time())
  for (i in 1:vcount(graph)) {
    df$gain[i] <- get_influence(graph, V(graph)[i], test_method=test_method, lt_threshold=prob)
  }
  # Arrange the data frame by marginal gains
  df <- arrange(df, desc(gain))

  # Modification for the framework
  if (criteria == "TIME") {
    # If the criteria is TIME, then the threshold should be set according to time
    criteria_value <- Sys.time()
    threshold <- criteria_value + seconds(threshold)
  } else {
    criteria_value <- 0
  }

  # Until the budget is met
  while (criteria_value < threshold) {
    top_row <- df[1,]
    u <- V(graph)[top_row$node]
    # If the flag is the size of the seed set so far, then add this node to the seed and remove from data frame
    if (top_row$flag == length(seed)) {
      seed <- c(seed, top_row$node)
      df <- df[-1,]
    }
    else {
      # Otherwise compute the marginal gain with this node
      current_influence <- get_influence(graph, V(graph)[seed], test_method=test_method, lt_threshold=prob)
      top_row$gain <- get_influence(graph, V(graph)[c(seed, u)], test_method=test_method, lt_threshold=prob) - current_influence
      # Store the length of seed in the flag
      top_row$flag <- length(seed)
      # Update the values for this row in data frame
      df[1,] <- top_row
      # Sort the data frame again by gain
      df <- arrange(df, desc(gain))
    }

    # Recalculate the new value for criteria to supply to loop condition
    if (criteria == "TIME") {
      # In case of TIME, the loop should end if the time threshold is up
      criteria_value <- Sys.time()
    } else if (criteria == "BUDGET") {
      # In case of BUDGET, the loop should end if the required budget threshold is met
      criteria_value <- length(seed)
    } else if (criteria == "INFLUENCE") {
      # In case of INFLUENCE, the loop should end if the current seed is sufficient according to threshold
      criteria_value <- get_influence(graph, V(graph)[seed], test_method=test_method, lt_threshold=prob)
    }
  }
  end <- as.numeric (Sys.time())
  output <- NULL
  output$influential_nodes <- V(graph)[seed]
  output$time <- (end - start)
  output$influence <- get_influence(graph, output$influential_nodes, test_method=test_method, lt_threshold=prob)
  output
}


## TEST
graph <- erdos.renyi.game(100, 0.2)
graph <- load_arxiv_collaboration_graph()
graph <- load_author_netscience_graph()
graph <- load_protein_barabasi_graph()

# TIME bound influence mining based on CELF
x <- celf_greedy_framework(graph, criteria="TIME", threshold=5, test_method="INFLUENCE_IC")
# BUDGET bound influence mining based on CELF
y <- celf_greedy_framework(graph, criteria="BUDGET", threshold=5, test_method="INFLUENCE_IC")
# SPREAD bound influence mining based on CELF
z <- celf_greedy_framework(graph, criteria="INFLUENCE", threshold=100, test_method="INFLUENCE_IC")

