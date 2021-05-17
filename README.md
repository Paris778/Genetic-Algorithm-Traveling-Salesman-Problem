# Genetic-Algorithm-Traveling-Salesman-Problem

## The Traveling Salesman problem

The Traveling Salesman Problem (TSP) is one of the most famous problems in
computer science. 
TSP is basically attempting to find the shortest complete tour through a series of points
(cities), starting and ending at the same point. Finding the shortest route that visits a
set of locations is an exponentially difficult problem. 
An brute-force search of all possible paths
would be guaranteed to find the shortest, but is computationally intractable for all but
small sets of locations. For larger problems, optimization techniques, such as GA, are
needed to intelligently search the solution space and find near-optimal solutions.
Mathematically, traveling salesman problem can be represented as a graph, where the
locations are the nodes and the edges (or arcs) represent direct routes between the
nodes. The weight of each edge is the distance between the nodes. It is a minimization
problem starting and finishing at a specified vertex after having visited each other vertex
exactly once. The goal is to find the path with the shortest sum of weights. Below, we
see a simple five-node graph:



![Alt text](https://github.com/Paris778/Genetic-Algorithm-Traveling-Salesman-Problem/blob/main/scnrSHots/Capture4.JPG "Title")

In reality, and for this probelm specifically, the nodes can be many-many more:

![Alt text](https://github.com/Paris778/Genetic-Algorithm-Traveling-Salesman-Problem/blob/main/scnrSHots/Capture2.jpg "Bes")

And using this implementation we are able to find the optimal path for this set of nodes.

![Alt text](https://github.com/Paris778/Genetic-Algorithm-Traveling-Salesman-Problem/blob/main/scnrSHots/Capture3.jpg "Best Solutio")

## Optimal parameters

For this implementation the Optimal parameters were:

### Population size : 50
For the Population Size parameter, it was decided that 50 was the optimal number in regards to a result to time ratio.
![Alt text](https://github.com/Paris778/Genetic-Algorithm-Traveling-Salesman-Problem/blob/main/scnrSHots/Capture5.JPG "Best Solution")
  
### Number of Iterations : 6000
This is a graph displaying the distance metric (y-axis) against the number of iterations (x-axis)
<br /><br />
![Alt text](https://github.com/Paris778/Genetic-Algorithm-Traveling-Salesman-Problem/blob/main/scnrSHots/Capture6.JPG "Best Solution")

### Mutation type: Flip mutation

Below are the results obtained for every type of mutation strategy implemented and tested
![Alt text](https://github.com/Paris778/Genetic-Algorithm-Traveling-Salesman-Problem/blob/main/scnrSHots/Capture7.JPG "Best Solution")
![Alt text](https://github.com/Paris778/Genetic-Algorithm-Traveling-Salesman-Problem/blob/main/scnrSHots/Capture8.JPG "Best Solution")
![Alt text](https://github.com/Paris778/Genetic-Algorithm-Traveling-Salesman-Problem/blob/main/scnrSHots/Capture9.JPG "Best Solution")
  
### Selection method : Elitist Selection

### Cross-over method : Single-Point Crossover
  
### Cross-over Probability:
It was determined that the best results were obtained with a crossover probability of *0.5-0.7*
### Mutation Probability:
It was determined that the best results were obtained with a mutation probability of *0.7-0.8*

## BONUS

During development, at some point the algorithm came up with a **harmonograph** as an optimal solution which I thought was very interesting.
<br /><br />
![Alt text](https://github.com/Paris778/Genetic-Algorithm-Traveling-Salesman-Problem/blob/main/scnrSHots/Harmonograph.JPG "Best Solution")


