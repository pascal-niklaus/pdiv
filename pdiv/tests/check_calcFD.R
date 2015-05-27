library(pdiv)
library(ape)

Getlength <- function(xtree, comp=NA) {
    if(!is.data.frame(comp))
        result <- Getlength.inner(xtree)
    if(is.data.frame(comp)){
        S <- tapply(comp[,2], comp[,1], function(x) length(x))
        FD <- tapply(comp[,2], comp[,1],
                     function(x)
                         Getlength.inner(list(xtree[[1]],
                                              xtree[[2]][!is.na(match(dimnames(xtree[[2]])[[1]], x)),]))) 
      FD <- FD/Getlength.inner(xtree)
      result <- data.frame(S=S, FD.new=FD)
  }
  result
}

Getlength.inner <- function(xtree) {   
    if(!is.matrix(xtree[[2]]))
        result <- 0
    if(is.matrix(xtree[[2]]))
        result = sum(xtree[[1]][colSums(xtree[[2]]) != 0 & colSums(xtree[[2]]) < length(xtree[[2]][,1])])
    result
}

Xtree <- function(h) {
    species.names <- h$labels
      
    H1 <- matrix(0, length(h$order), 2 * length(h$order) - 2)
    l <- vector("numeric", 2 * length(h$order) - 2)
    for(i in 1:(length(h$order) - 1)) {
        if(h$merge[i, 1] < 0) {
            l[2 * i - 1] <- h$height[order(h$height)[i]]
            H1[ - h$merge[i, 1], 2 * i - 1] <- 1
        }
        else {
            l[2 * i - 1] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 1]]]
            H1[, 2 * i - 1] <- H1[, 2 * h$merge[i, 1] - 1] + H1[, 2 * h$merge[i, 1]]
        }
        if(h$merge[i, 2] < 0) {
            l[2 * i] <- h$height[order(h$height)[i]]
            H1[ - h$merge[i, 2], 2 * i] <- 1
        }
        else {
            l[2 * i] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 2]]]
            H1[, 2 * i] <- H1[, 2 * h$merge[i, 2] - 1] + H1[, 2 *h$merge[i, 2]]
        }
    }
    dimnames(H1) <- list(species.names,NULL)  
    list(h2.prime=l, H1=H1)
}

data(traits.example)
data(composition.example)

print(traits.example)
print(composition.example)

dimnames(traits.example) <- list("species"=as.character(traits.example[,1]),
                                 "traits.example"=dimnames(traits.example)[[2]])
distances <- dist(traits.example[,-1],"euclidean")
tree <- hclust(distances,"complete")
xtree <- Xtree(tree)

FD1 <- tapply(composition.example[,2], composition.example[,1],
               function(x) Getlength(list(branchlengths=xtree$h2.prime,
                                          sppxbranch=xtree$H1[!is.na(
                                            match(dimnames(xtree$H1)[[1]],x)),])))

# > FD1
#        1        2        3        4        5        6        7        8 
# 3.750272 2.773446 3.418143 4.690748 0.000000 0.000000 0.000000 0.000000 

FD2 <- Getlength(xtree, composition.example)$FD.new


y <- calcFD(composition.example,traits.example,which="stdFD",shrink.tree = FALSE);

# print(FD2);

stopifnot(all(identical(as.numeric(FD2),as.numeric(y$stdFD))));

