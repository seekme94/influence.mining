# Project: influence.mining

This project provides various functions to perform influence mining operations on graphs.

## Introduction
The project before you is part of a doctorial study by the author of this repository on Influence mining in complex networks. This package is written and maintained in R language and uses igraph library as a mandatory dependency.
Influence mining is a branch of network science, in which the objective is to identify the few most influential nodes in a network, which are capable of propagating information. The applications of influence mining (a.k.a influence maximization) can be found in viral marketing, social networks, traffic control, epidemiology and more.

The aim of creating this package is to help other researchers and network scientists. If this package helped you in publishing your work, then please cite this repository.

This document serves as official guide to use this package. However, some detail might be missed or intentionally skipped. You are free to fork, suggest and make changes and create pull requests to make this package more useful.

## Installation
In order to install this package, you first need to install **devtools** package:
```
install.packages("devtools")
```
Then install this package from GitHub:
```
devtools::install_github("seekme94/influence.mining")
```


## Package tutorial

```
influence_maximization.R
```
All influence mining functions are provided in this source file. The primary interfacing function, **influence** is a wrapper function, which calls other functions. For most of the users, this should be enough for applications.
```
influence (graph, budget, prob, steps, optimal_solution, test_method, heuristic, centrality_method, parallel, logging)
```
- **graph** is the igraph object
- budget number of influential nodes to be fetched. Default value is 1
- **prob** the probability at which a node influences its neighbours
- **steps** (currently not implemented) is the time steps for which, the diffusion process should run. Provide NULL for exhaustive run. Default value is 1.
- **test_method** specifies the method to measure influence. Value MUST be "RESILIENCE", "INFLUENCE_IC" or "INFLUENCE_LT"
- **heuristic** specifies the heuristic method used for influence calculation. Required only when optimal_solution is FALSE
- **centrality_method** (optional) is the centrality algorithm to use when heuristic is "CENTRALITY" or "ADAPTIVE_CENTRALITY". Value must be "DEGREE", "BETWEENNESS", "CLOSENESS" or "EIGENVECTOR"
- **parallel** when true, executes the funtion using multiple CPU cores. Default value is TRUE
- **optimal_solution** should be TRUE if influential nodes are to be derived using optimal algorithm. Caution! This is the slowest apporach
- **logging** when true, a complete log is stored in output.log file

Returns an object containing:
1.  Vector of influential nodes.
2.  Measure of influence. 
3.  Elapsed time in milliseconds.

Example 1.
On a random graph of 200 nodes, find most influential nodes with Greedy method using Resilience test to quantify the influence.
```
require(igraph)
require(logging)
graph <- erdos.renyi.game(200, 0.2)
influence(graph, budget=5, prob=0.5, steps=1, optimal_solution=FALSE, test_method="RESILIENCE", heuristic="GREEDY")
```
Output
```
$`influential_nodes`
+ 5/200 vertices, from ab47fd1:
[1] 1 2 3 4 5

$time
[1] 1.96186

$influence
[1] 5
```

Example 2:
On a small world graph of 200 nodes, find most influential nodes with Pagerank heuristic using Resilience test to quantify the influence.
```
graph <- generate_scale_free(200)
influence(graph, budget=5, prob=0.5, steps=1, optimal_solution=FALSE, test_method="RESILIENCE", heuristic="PAGERANK")
```
Output
```
$`influential_nodes`
+ 5/200 vertices, named, from 29843b0:
[1] 13 15 21 3  9 

$influence
[1] 119

$time
[1] 0.0009799004
```

Example 3:
On the famous karate club graph, find the globally optimal influential nodes using Resilience test to quantify the influence.
```
graph <- generate_scale_free(200)
influence(graph, budget=5, prob=0.5, steps=1, optimal_solution=FALSE, test_method="RESILIENCE", heuristic="PAGERANK")
```
Output
```
$`influential_nodes`
+ 5/200 vertices, named, from 29843b0:
[1] 13 15 21 3  9 

$influence
[1] 119

$time
[1] 0.0009799004
```


### Libraries used
See the [DESCRIPTION file](https://github.com/seekme94/influence.mining/blob/master/DESCRIPTION)


### References:
[1] Kempe, D., Kleinberg, J., & Tardos, É. (2003). Maximizing the Spread of Influence through a Social Network. In Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining - KDD ’03 (p. 137). New York, New York, USA: ACM Press. doi:10.1145/956755.956769


<!-- badges: start -->
[![R build status](https://github.com/seekme94/influence.mining/workflows/R-CMD-check/badge.svg)](https://github.com/seekme94/influence.mining/actions)
<!-- badges: end -->
