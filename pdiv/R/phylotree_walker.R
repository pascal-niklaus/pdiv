############################################################
##
## Functions to process "phylo" structures containing
## a tree structure as provided by the "ape" package
## Note that these are also used to calculate functional
## trait diversity sensu Petchey & Gaston.

############################################################
##
## walk_to_root: 
##
##    walk from tip to root, returning edge indices

.walk_to_root <- function(tr, tip) {
  r <- c();
  repeat{
    edge <- which(tr$edge[,2]==tip)
    if(length(edge)==0) 
      break;
    r <- append(r,edge);      
    tip <- tr$edge[edge,1];    
  }
  r;  
}

############################################################
##
## on_way_to_root:
##
##    test if edge given is part of the way to the root
##    this is done starting from all tree tips

.on_way_to_root <- function(br, edge) {
  sapply(br,function(x) { edge %in% x });
}

############################################################
##
## find_common_edges:
##
##    return edges that are common to all tips 
##    (when travelling to root)

.find_common_edges <- function(br) {
  shortest <- order(sapply(br,length))[1];
  to.remove <- integer(0);
  for(section in br[[shortest]]) {
    if(all(.on_way_to_root(br,section)))
      to.remove <- append(to.remove,section);
  }
  to.remove;
}

############################################################
##
## remove_common:
##
##    remove common edges from branch list
##    this function is used to remove the unused parts from
##    a tree when working with a subtree. 

.remove_common <- function(br) {
  to.remove <- .find_common_edges(br);
  lapply(br,function(x) setdiff(x,to.remove));
}

