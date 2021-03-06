# R package 'pdiv'

## Purpose

This R package contains functions related broadly to "biodiversity" topics,
as opposed to general functions (which are in the package `pascal`) and "geography"-
related function (which are in the package `pgeo`)

## Groups of functions 

### Additive and tripartite partitioning

- `addpart` partitions biodiversity effects following Loreau and Hectors scheme.
- `tripart` partitions biodiversity effects following Fox' scheme.

Both functions support parallelization.

### Biodiversity indices

These functions calculate indices which are part of Hill's series of biodiversity numbers, 
plus a few related indiced.

- `hill.Na(x,a=1)`
- `hill.Ha(x,a=1)`
- `shannon.H(x)`
- `simpson.D(x)`
- `evennes.1overD(x)`

### Interconversion between data formats representing community composition

- `communityFromMatrix(...)`
- `communityToMatrix(...)`

### Niche breadth and overlap

- `levinsB` calculates niche breadth as defined by Levins
- `prop.sim` calculates Colwell and Futuyma's proportional similarity index 

### Functional and phylogenetic diversity

These functions calculate functional and phylogenetic diversity as total branch length
of a functional trait dendrogram or a phylogenetic tree, following Petchey \& Gastons 2002 
scheme (including their 2006 update). Also, a function to calculate mean phylogenetic distances,
with and without abundance weighing, is provided (Warwick and Clark 1995).

- `calcFD(...)`
- `calcPD(...)`
- `calcMPD(...)`

## Installation

* download the ready-built package from the [pkgs directory](https://github.com/pascal-niklaus/pdiv/tree/master/pkgs)
* use `install_github`:  
`library(devtools)`  
`install_github("pascal-niklaus/pdiv/pdiv")`


