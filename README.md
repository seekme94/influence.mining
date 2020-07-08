# Project: influence.mining
### Version: 0.2.0

This project provides various functions to perform influence mining operations on graphs.


```
Influence.R
```
The file contains source code for implementation of two basic influence mining models: Independent Cascade model and Linear Threshold model[1]
```
influence (graph, seed, budget, steps, model, maximize, seed_method, prob)
```
This is a wrapper function to call influence_LT and influence_IC functions
- graph: is the igraph object
- budget: defines what percentage of most influential nodes out of all nodes is required as output. Default value is 1
- seed: (optional) is a set of seed (initial nodes). If this parameter is NULL, then seed_method parameter should be given
- steps: is the time steps for which, the diffusion process should run. If exhaustive run is required, provide a high value (like 100). Default value is 1
- model: is influence model to run the dataset on. Value MUST either be "LT" or "IC"
- maximize: should be TRUE if influential nodes are to be derived using Greedy algorithm
- seed_method: is the selection method for seed (initial nodes). Value can be "random", "degree", "closeness", "betweenness", "coreness", "eigenvector", "a-degree", "a-closeness", "a-betweenness", "a-coreness", "a-eigenvector". Default value is "random"
- prob: is the probability of activation of a neighbour node. This is applicable only to IC model currently. Default value is 0.5

> Output: summary of influence process, including no. of nodes, edges, seed set size, nodes influenced and time taken

```
influence_LT
```
This function calculates influence (number of nodes in the network expected to be activated) under Linear Threshold model. For parameters, see influence function.

```
influence_IC
```
This function calculates influence (number of nodes in the network expected to be activated) under Independent Cascade model. For parameters, see influence function.

```
select_seed
```
This function returns a set of nodes, to be used as seed in influence functions on the basis of given seed selection method
- G: a graph object of library *igraph*
- k: percentage of seed nodes from the network to be chosen
- seed_method: see influence function

> Output: subset vector of nodes in a graph

```
select_adaptive_seed
```
This function returns a set of nodes, to be used as seed in influence functions on the basis of given adaptive method for seed selection
- G: a graph object of library *igraph*
- k: percentage of seed nodes from the network to be chosen
- seed_method: see influence function

> Output: subset vector of nodes in a graph

```
find_communities
```
This method finds communities in the given graph and returns the graph after adding a vector "group" to its vertices
- G: a graph object of library *igraph*
- method: is the method to generate communities. Available algorithms are "multilevel", "edgebetweenness", "fastgreedy", "eigenvector", "spinglass", "walktrap", "labelpropagation", "clique", "largescale"

> Output: graph object with additional vector "group" to vertices

```
community.significance.test
```
This function performs a Wilcoxon rank-sum test on the "internal" and "external" degrees of a community in order to quantify its significance.


### Examples:
1. Calculate influence under defaults (model="LT", budget=5, steps=1 and seed_method="random")
```
influence(edgesFile="C:/Datasets/twitter_edges.csv")
```
2. Calculate influence under IC model, budget=10% for 2 time steps and seed_method="random"
```
influence(edgesFile="C:/Datasets/twitter_edges.csv", budget=10, steps=2, model="IC")
```
3. Calculate influence under IC model to select 10% nodes for 2 time steps and seed selection criteria to be nodes with highest degree
```
influence(edgesFile="C:/Datasets/twitter_edges.csv", budget=10, steps=2, model="IC", seed_method="degree")
```
4. Calculate influence under LT model to select 5% nodes for 1 time steps and seed selection criteria to be nodes with highest betweenness
```
influence(edgesFile="C:/Datasets/twitter_edges.csv", seed_method="betweenness")
```
### Libraries used
jsonlite, uuid, sampling, digest, RWeka, doMC, snow, doSNOW, iterpc, foreach, igraph, caret, e1071, party, rpart, rpart.plot, randomForest, RColorBrewer, nnet, rattle, ggplot2, Rcpp


### References:
[1] Kempe, D., Kleinberg, J., & Tardos, É. (2003). Maximizing the Spread of Influence through a Social Network. In Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining - KDD ’03 (p. 137). New York, New York, USA: ACM Press. doi:10.1145/956755.956769
