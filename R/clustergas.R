################################################
# "clustergas.R" is a Genetic Algorithm  for
# building Hierarchical Clustering from datasets.
#
# Version:     1.0                                                  
# Author:      José A. Castellanos G.
#              Department of Computer Science
#              University of Valladolid
#  
# Last update: 13/06/2012  
#
# Note: this version is to upload to CRAN.
# This version does not create or read files 
# from the current directory. 
################################################


# Creating a new environment for constants and global variables
  myenv <- new.env()

# Constants  and variables
############

  myenv$OBJ.COUNT       <- 10 # Dataset objects number. 
  myenv$PART            <- 1/2 # Part of dendrogram to remove...
  myenv$PCENT.MUT       <- 0.50 # (0.15)Percent of time number to
                          # merge two different clusters (mutation)...
  myenv$DIFF            <- 0.03 # Differences max. among clusters (mutation)...  
  myenv$ITE.CROSS       <- 200 # Number of iteration on the clustering
                          # chosen in the two dendrograms of
                          # the crossover3.  
  #---------------------------------------------------------
  
  myenv$POPULATION_SIZE <- 5
  myenv$CHROM_LENGTH    <- myenv$OBJ.COUNT - trunc(myenv$OBJ.COUNT * myenv$PART) - 2
  myenv$PCROSS          <- 0.40
  myenv$PMUT            <- 0.1 #0.60   #0.05
  myenv$MAX_GEN         <- 10
  myenv$MODIFI          <- c()
  
  # R Package call
  ###############
  
  #library(cluster)
  
  
  # Global Variables
  ############
    
  clustervector  <- c(0) # Represents a cluster (a subset).
  clustering     <- list(clustervector) # Represents a clustering.
  dendogram      <- list(clustering) # A tree with each clustering 
                                     # (one individual).
  population     <- list(dendogram) # Population of the GA.
  
  myenv$fitness.values <- as.double(0) # Fitness function values vector.
  myenv$prob.values    <- as.double(0) # Probability vector by "fitness.values".
  
  # Disimilarity matrix.
  obj.data        <- matrix(sample(1:6, 10*10, TRUE), 10, 10, TRUE, NULL)
  myenv$DISTANCE  <- as.matrix(dist(obj.data, method = "euclidean", TRUE, TRUE))
  myenv$MAX_DIST  <- max(as.dist(myenv$DISTANCE))
  myenv$MIN_DIST  <- min(as.vector(as.dist(myenv$DISTANCE)))
  myenv$DIST_ONE  <- 95 * (myenv$MAX_DIST - myenv$MIN_DIST) / 100 + myenv$MIN_DIST   

     
# R Functions
############

####################################################
# "init.variables" function.
# This function creates and initializes constants 
# and global variables used by the functions of this 
# package.
####################################################

init.variables <- function(dist.matrix, pop.size = 10, gen.numb = 10, 
                           p.cross = 0.40, p.mut = 0.1, part = 1/2)
{
 # Updating global variable and constant values.
 myenv$OBJ.COUNT <- nrow(dist.matrix)
 myenv$PART <- part
 temp <- myenv$OBJ.COUNT - trunc(myenv$OBJ.COUNT * myenv$PART) - 2
 myenv$CHROM_LENGTH <- temp 
 myenv$PCROSS <- p.cross
 myenv$PMUT <- p.mut #0.60   #0.05
 myenv$POPULATION_SIZE <- pop.size
 myenv$MAX_GEN <- gen.numb
 myenv$DISTANCE  <- dist.matrix
 myenv$MAX_DIST <- max(as.dist(myenv$DISTANCE)) 
 myenv$MIN_DIST <- min(as.vector(as.dist(myenv$DISTANCE)))
 temp <- 90 * (myenv$MAX_DIST - myenv$MIN_DIST) / 100 + myenv$MIN_DIST
 myenv$DIST_ONE <- temp
 myenv$MODIFI <- integer(pop.size)
 myenv$DIFF <- 0.03 # Differences max. among clusters (mutation ops)...  
 myenv$ITE.CROSS <- 200  # Number of iterations in the crossover3
 myenv$PCENT.MUT <- 0.30 # only for mutation2

 # Other variables
 myenv$fitness.values <- as.double(0) # Fitness function values vector.
 myenv$prob.values <- as.double(0) # Probability vector for "fitness.values".
}


init.varsearch <- function(dist.matrix, pop.size, gen.numb, part){
   myenv$OBJ.COUNT <- nrow(dist.matrix)
   myenv$POPULATION_SIZE <- pop.size
   myenv$MAX_GEN <- gen.numb
   myenv$DISTANCE  <- dist.matrix
   myenv$MAX_DIST <- max(as.dist(myenv$DISTANCE)) 
   myenv$MIN_DIST <- min(as.vector(as.dist(myenv$DISTANCE)))
   temp <- 90 * (myenv$MAX_DIST - myenv$MIN_DIST) / 100 + myenv$MIN_DIST
   myenv$DIST_ONE <- temp
   myenv$PART <- part
   temp <- myenv$OBJ.COUNT - trunc(myenv$OBJ.COUNT * myenv$PART) - 2
   myenv$CHROM_LENGTH <- temp  
}
 
#########################################################
# "random.int" function
# random of vector numbers "vector.x".
# "long" is numbers amount selected. 
# if "repl=TRUE" then selection has repeated numbers.
#########################################################

random.int <- function(vector.x, long, repl = FALSE) 
 sample(vector.x, long, repl)

#######################################################
# "init.population" function.
# Creates the initial population of dendograms (individuals).
#######################################################

init.population <- function(population) 
{
 dendo <- dendogram; init <- trunc(myenv$OBJ.COUNT * myenv$PART)
 for (i in 1:myenv$POPULATION_SIZE) 
  {
   # Creates the first clustering and adds at the end it
   # a vector with the numbers of merged clusters to 
   # form next clustering. 
   clust <- c(as.list(1:myenv$OBJ.COUNT), list(c(0, 0)))
   # Dendograms constructing. 
   for (j in 1:(myenv$OBJ.COUNT - 2))
    {
     size <- myenv$OBJ.COUNT - j + 1  # Size of clustering before do it. 
     randm <- random.int(1:size, 2) # Generates two cluster positions
                                    # to merge it (it's not number repeating).  
     clust[[randm[1]]] <- c(clust[[randm[1]]], 
                            clust[[randm[2]]]) # Merging clusters...
     # Reducing by one, the clustering size.
     clust[[randm[2]]] <- clust[[size]]
     clust[[size]] <- randm # Update merged clusters numbers. 
     length(clust) <- size 
     # Adding the clustering. 
     if (j > init) dendo[[j - init]] <- clust 
    }
   population[[i]] <- dendo 
  } 
 return (population)
}

#-----------------------------------------------------------------------

################################################
# "mutation" function.
# Mutation of a dendrogram in a fixed level. 
################################################
 
mutation <- function(dendogram, level = 0, comb = F,
                     homog.f = "cluster.fitness")
{
 if (level == 0)
  level <- random.int(1:(myenv$CHROM_LENGTH - 1), 1) # Level choice for mutation.
 clust <- dendogram[[level]] # Taking clustering of level "level".
 for (i in (level + 1):myenv$CHROM_LENGTH)
   {
    size <- myenv$CHROM_LENGTH - i + 3 # Clustering Size (numb. of clusters).
    randm1 <- random.int(1:size, 2) # Choices two clusters to merge. 
    clust[[randm1[1]]] <- c(clust[[randm1[1]]], 
                            clust[[randm1[2]]]) # Merging clusters choice.  
    clust[[randm1[2]]] <- clust[[size]]  # Moving last one cluster toward 
                                         # "randm1[2]" position.
    clust[[size]] <- randm1 # Update merged clusters numbers.                                     
    length(clust) <- size  # Reducing list (clustering) size in 1.
    dendogram[[i]] <- clust    # Updates dendogram level.     
   }
 return (dendogram)
}

################################################
# "mutation2" function.
# It makes the mutation in a fixed level, 
# choising two clusters a myenv$PCENT.MUT of times. 
################################################
 
mutation2 <- function(dendogram, level = 0, comb = F, 
                      homog.f = "cluster.fitness")
{ 
 # Choicing the level for the mutation.
 if (level == 0)
  level <- random.int(1:(myenv$CHROM_LENGTH - 1), 1) 
 # Taking clustering level "level".
 clust <- dendogram[[level]] 
 size <- myenv$CHROM_LENGTH - (level + 1) + 3
 pos <- dendogram[[level + 1]][[size]]
 valuef <- do.call(homog.f, list(c(clust[[pos[1]]], clust[[pos[2]]])))
 #valuef <- cluster.fitness(c(clust[[pos[1]]], clust[[pos[2]]]))
 for (i in (level + 1):myenv$CHROM_LENGTH)
   {
    # Clustering Size (numb. of clusters).
    size <- myenv$CHROM_LENGTH - i + 3 
    if (!comb) temp <- round(size * myenv$PCENT.MUT, 0)
     else temp <- round((size * (size - 1)/2) * myenv$PCENT.MUT, 0)
    if (temp < 1) temp <- 1
    # Choices two clusters to merge.
    for (j in 1:temp)
     {
      randm1 <- random.int(1:size, 2)
      tclust <- c(clust[[randm1[1]]], 
                  clust[[randm1[2]]])
      valuef2 <- do.call(homog.f, list(tclust))
      #valuef2 <- cluster.fitness(tclust)                  
      if ((valuef2 < valuef) | 
         (round(abs(valuef2 - valuef), 2) <=  myenv$DIFF)) 
       {
        valuef <- valuef2
        pos <- randm1 
       }           
     } 
    valuef <- +Inf     
    # Merging clusters choice. 
    clust[[pos[1]]] <- c(clust[[pos[1]]], 
                         clust[[pos[2]]])   
    # Moving last one cluster toward 
    # "randm1[2]" position.                       
    clust[[pos[2]]] <- clust[[size]]  
    # Update merged clusters numbers.
    clust[[size]] <- pos                                      
    # Reducing list (clustering) size in 1.
    length(clust) <- size  
    # Updates dendogram level.
    dendogram[[i]] <- clust         
   }
 return (dendogram)
}

################################################
# "mutation3" function.
# It makes a mutation by switching two objects
# between the two clusters of the last level.
# Afterwards, that change is propagated to the low
# levels.
################################################

mutation3 <- function(dendogram, level = 0, comb = F,
                      homog.f = "cluster.fitness")
{
 # Look for an element in a clustering.
 search.clust <- function(elem, clustering, size)
 {
  i <- 1
  while (i <= size)
   if (elem %in% clustering[[i]]) break()
    else i <- i + 1
  return (c(i, which(clustering[[i]] == elem)))
 }
 
 clust <- dendogram[[myenv$CHROM_LENGTH]]
 pos1 <- random.int(1:(temp <- length(clust[[1]])), 1)
 pos2 <- random.int(1:(myenv$OBJ.COUNT - temp), 1)
 elem1 <- clust[[1]][pos1]
 elem2 <- clust[[2]][pos2]
 clust[[1]][pos1] <- elem2
 clust[[2]][pos2] <- elem1
 dendogram[[myenv$CHROM_LENGTH]] <- clust
 size <- 2 # Clustering Size (numb. of clusters of the last level).
 for (i in (myenv$CHROM_LENGTH - 1):1)
   {
    size <- size + 1
    clust <- dendogram[[i]]
    vpos1 <- search.clust(elem1, clust, size) # cluster and inside cluster position
    vpos2 <- search.clust(elem2, clust, size) # of elem1 and elem2 in "clust".
    clust[[vpos1[1]]][vpos1[2]] <- elem2
    clust[[vpos2[1]]][vpos2[2]] <- elem1
    dendogram[[i]] <- clust    # Updates the current dendogram level.
   }
 return (dendogram)
}

#----------------------------------------------------------------------

#################################################
# "obj.dist" function.
# Return an element of class "dist"  
# belongs to (i, j) element of 
# disimilarity simetrix matrix.
# "n" is rows number.
#################################################
 
obj.dist <- function(i, j, n)
{
 if (i < j) {temp <- i; i <- j; j <- temp}
 return (myenv$DISTANCE[(j - 1) * n - j * (j - 1) / 2 + i])
}

#################################################
# "cluster.fitness" function.
# Fitness computing of a cluster, it's
# mean distance of cluster elements.   
#################################################

cluster.fitness <- function(cluster1)
{
 if (is.na(cluster1[1])) return (myenv$MAX_DIST) # if cluster is empty.
  else if (is.na(cluster1[2])) return (myenv$DIST_ONE) # if cluster has 
   else                                          # size 1.
    {
     c <- length(cluster1)
     return (sum(myenv$DISTANCE[cluster1, cluster1]) / (c * (c - 1))) 
    }
}

#################################################
# "cluster.dist" function.
# Distance between two clusters,
# the smallest distance between theirs.
#################################################

cluster.dist <- function(cluster1, cluster2)
 min(myenv$DISTANCE[cluster1, cluster2]) # disjoin clusters. 

cluster.dist.mean <- function(cluster1, cluster2)
 mean(myenv$DISTANCE[cluster1, cluster2]) # disjoin clusters.
 
###################################################
# "clustering.fitness" function.
# Clustering fitness, it's the diference between 
# distances mean and clusters fitness mean.
###################################################  
 
clustering.fitness <- function(clustering)
{                                    
 mean1 <- 0; k <- length(clustering) - 1; betha <- 0
 # mean of cluster homogeneity in the clustering.
 mean0 <- mean(sapply(clustering[-k-1], "cluster.fitness"))
 # mean of distances between all clusters pair.
 for (i in 1:(k - 1))
  for (j in (i + 1):k)
   mean1 <- mean1 + cluster.dist(clustering[[i]], 
                                 clustering[[j]])
 mean1 <- 2 * mean1 / (k * (k - 1))
 # Aptitud parametrizada, función sigma.
 #if (sigma) {
 # td <- round(myenv$OBJ.COUNT / k, 0)
 # m <- sapply(clustering[-k-1], "length") - td
 # sigma <- sqrt(sum(m^2) / k) 
 # betha <- ((k - 1) * td^2 + (myenv$OBJ.COUNT - td)^2) / k
 #}
 #print(paste("Beta - Sigma: ", sqrt(betha) - sigma))
 # "myenv$MAX_DIST" avoides negatives fitness values.
 #return ((myenv$MAX_DIST + mean1 - mean0) + (sqrt(betha) - sigma)) 
 return (myenv$MAX_DIST + mean1 - mean0)
}

#####################################################
# "fitness.mean2" function.
# dendogram fitness, it's fitness mean of clustering.
#####################################################
 
fitness.mean2 <- function(dendogram)
 mean(sapply(dendogram, "clustering.fitness"))

#####################################################
# "fitness.list" function.
# Return a vector with dendogram clustering fitness. 
#####################################################

fitness.list <- function(dendogram, homog.f = "cluster.fitness")
{
 sum.fitness <- 0
 # Computing first fitness of dendogram.
 fc <- clustering.fitness(dendogram[[myenv$CHROM_LENGTH]])
 # Mean of all fitness of all clusters.
 H <- mean(sapply(dendogram[[myenv$CHROM_LENGTH]][-3], homog.f))
 # Mean of all distances of all clusters.                
 S <- fc - myenv$MAX_DIST + H
 fitness.vector <- as.vector(fc)
 # Size of update clustering. 
 k <- 3; 
 for (i in (myenv$CHROM_LENGTH - 1):1)
  {
   # Position of two merged clusters. 
   j <- dendogram[[i + 1]][[k]]
   # Get the two clusters splited
   c1 <- dendogram[[i]][[j[1]]]; c2 <- dendogram[[i]][[j[2]]] 
   # Size of two splited clusters.
   l <- c(length(c1), length(c2))
   # Count of distances of used clusters.
   kj <- c(l[1]*(l[1] - 1) / 2, l[2] * (l[2] - 1) / 2,
           (l[1] + l[2]) * (l[1] + l[2] - 1) / 2)
   # Fitness of clusters t1 and t2
   f <- c(do.call(homog.f, list(c1)), do.call(homog.f, list(c2)))
   #f <- c(cluster.fitness(c1), cluster.fitness(c2))
   # Distance mean between clusters j1 and j2.
   d.mean <- cluster.dist.mean(c1, c2)
   # Clustering homogenity fitness. 
   H <- ((k - 1) * H + f[1] + f[2] - (kj[1] * f[1] + kj[2] * f[2] + 
        l[1] * l[2] * d.mean) / kj[3]) / k
   # Count of distance betwen clusters.       
   g <- k * (k - 1) / 2
   # Distance min. between clusters j1 and j2.
   d.min <- cluster.dist(c1, c2)
   # Removes two clusters and position...
   clust <- dendogram[[i]][-c(j, k + 1)]
   # Distances of both c1 and c2 to other clusters.
   t1 <- as.vector(0)
   for (d in  1:(k - 2))
    {
     t1[d] <- cluster.dist(c1, clust[[d]])
     t1[k - 2 + d] <- cluster.dist(c2, clust[[d]])
    }
   plus <- sum(t1)
   # Distance to take away...
   minus <- 0
   for (d in  1:(k - 2)) 
    minus <- minus + (t1[d] + t1[k - 2 + d] - 
                      abs(t1[d] - t1[k - 2 + d])) / 2
   # Clustering separability fitness.
   S <- ((g - k + 1) * S - minus + d.min + plus) / g
   fitness.vector <- c(myenv$MAX_DIST + S - H, fitness.vector)
   k <- k + 1                    
  }
 return (fitness.vector)
}

#####################################################
# "fitness.mean" function.
# dendogram fitness, it's fitness mean of clustering.
#####################################################

fitness.mean <- function(dendogram)
 mean(fitness.list(dendogram))

#####################################################
# "fitness.max" function.
# dendogram fitness, it's fitness max of clustering.
#####################################################
 
fitness.max <- function(dendogram)
 max(fitness.list(dendogram))

####################################################
# "S1" function.
# Gives the separation betwen two clusters, it not
# perform centroid.
####################################################

S1 <- function(cluster1, cluster2)
 mean(myenv$DISTANCE[cluster1, cluster2]) # disjoin clusters.
 
####################################################
# "a1" function.
# Gives the object average distance "obj" to other  
# objects in the same cluster.
####################################################

a1 <- function(obj, clustering)
{
 i <- 1; c <- length(clustering) - 1
 while (i <= c)
  if (obj %in% clustering[[i]]) break()
   else i <- i + 1
 cluster1 <- setdiff(clustering[[i]], obj)
 if (is.na(cluster1[1])) return (myenv$DIST_ONE)
 
 return (S1(c(obj), cluster1))  
}

####################################################
# "b1" function.
# It gives the object average distance "obj" to objects 
# in its nearest neighbor cluster.
####################################################

b1 <- function(obj, clustering)
{
 least <- pos <- Inf
 for (i in 1:(length(clustering) - 1))
  if (!(obj %in% clustering[[i]]))
   {
    d <- S1(c(obj), clustering[[i]])
    if (d < least) least <- d
   }
   
 return (least)  
} 

####################################################
# "silhouette.obj" function.
# Gives the Silhouette of a object in a clustering.
####################################################

silhouette.obj <- function(obj, clustering)
{
 a.g <- a1(obj, clustering)
 b.g <- b1(obj, clustering)
 # Computing the max{a.g, b.g}
 major <- (a.g + b.g + abs(a.g - b.g)) / 2 
 # Return the Silhouette of "obj"
 return (myenv$MAX_DIST + (b.g - a.g) / major)
}

####################################################
# "silhouette.clust" function.
# Gives the Silhouette mean of all objects in a 
# clustering.
####################################################

silhouette.clust <- function(...)
{
 add <- 0
 for (i in 1:myenv$OBJ.COUNT)
  add <- add + silhouette.obj(i, ...)
 # add <- sum(sapply(1:myenv$OBJ.COUNT, silhouette.obj, clustering)) 
 return (add / myenv$OBJ.COUNT) 
}

#####################################################
# "fitness.mean.silhouette" function.
# dendogram fitness, it's the mean fitness  of all 
# clustering.
#####################################################

fitness.mean.silhouette <- function(dendogram)
 mean(sapply(dendogram, "silhouette.clust"))

#####################################################
# "fitness.max.silhouette" function.
# dendogram fitness, it's fitness max of all 
# clustering.
#####################################################

fitness.max.silhouette <- function(dendogram)
 max(sapply(dendogram, "silhouette.clust"))
 
#####################################################
# Fitness functions from the cluster validation 
# measures of homogeneity and separation
#####################################################
 
fitness.meanC.HS <- function(clustering)
 max(myenv$DISTANCE) + S.ave(clustering) - H.ave(clustering)
 
sigma.clust.HS <- function(clustering)
 max(myenv$DISTANCE) + S.ave(clustering) - H.ave(clustering) + 
 clust.sigma(clustering)

fitness.meanD.HS <- function(dendogram)
 mean(sapply(dendogram, fitness.meanC.HS))

sigma.meanD.HS <- function(dendogram)
 mean(sapply(dendogram, fitness.meanC.HS) +
      sapply(dendogram, "clust.sigma"))
       
#####################################################
# "sigma.clust.fitness" function.
# it is a objetive function to optimize the element  
# number inside each cluster in a clustering.
#####################################################

clust.betha <- function(size, tt)
 sqrt(((size - 1) * tt^2 + (myenv$OBJ.COUNT - tt)^2) / size)

clust.sigma <- function(clustering)
{
 size <- length(clustering) - 1
 # estimated of the cluster number in each cluster
 tt <- round(myenv$OBJ.COUNT / size, digits = 0) 
 # betha is the maximal value that sigma can reach.
 betha <- clust.betha(size, tt) 
 #print(betha)
 # cluster number of each cluster
 vect <- sapply(clustering[-(size + 1)], length)
 sigma <- sqrt(mean((vect - tt)^2))
 #print(sigma) 
 return (betha - sigma)
}   

sigma.clust.fitness <- function(clustering)
 clustering.fitness(clustering) + clust.sigma(clustering)

#####################################################
# "dendo.sigma1" function.
# objetive function to estimate cluster number  
# of the whole dendrogram, it is more quick than 
# the before.
##################################################### 
 
sigma.list <- function(dendogram)
{
 vbetha <- vsigma <- c()
 size <- 2; clust <- dendogram[[myenv$CHROM_LENGTH]]
 # estimated of the cluster number in each cluster
 t1 <- round(myenv$OBJ.COUNT / size, digits = 0) 
 # betha is the maximal value that sigma can reach.
 vbetha[myenv$CHROM_LENGTH] <- clust.betha(size, t1) 
 vsigma[myenv$CHROM_LENGTH] <- (clust.sigma(clust) - vbetha[myenv$CHROM_LENGTH])^2
 vpos <- clust[[size + 1]] # position of the clusters splited.
 for (i in (myenv$CHROM_LENGTH - 1):1)
  {
   clust <- dendogram[[i]]; size <- size + 1
   ms <- length(clust[[vpos[1]]]); mp <- length(clust[[vpos[2]]])
   t2 <- round(myenv$OBJ.COUNT / size, digits = 0)
   #print(vsigma)
   vsigma[i] <- vsigma[i + 1] * (size - 1) - (ms + mp)^2 + ms^2 + mp^2 +
                2 * myenv$OBJ.COUNT *(t1 - t2) + size * t2^2 - (size - 1) * t1^2
   vsigma[i] <- vsigma[i] / size
   vbetha[i] <- clust.betha(size, t2)
   # updating for the next iteration
   t1 <- t2             
   vpos <- clust[[size + 1]]
  }
 return (vbetha - sqrt(vsigma)) 
} 
 
#####################################################
# "sigma.dendo.fitness" function.
# dendogram fitness, it's fitness mean of clustering
# more another function to check the element number.
#####################################################
 
sigma.dendo.fitness <- function(dendogram)
 mean(sapply(dendogram, "clustering.fitness") +
      sapply(dendogram, "clust.sigma"))

#####################################################
# "sigma.dendo.fitness1" function.
# dendogram fitness, it's fitness mean of clustering.
#####################################################

sigma.dendo.fitness1 <- function(dendogram)
 mean(fitness.list(dendogram) + sigma.list(dendogram))

######################################################
# "allfitness" function.
# Gives all individuals fitness in the population.
# there are four available, "fitness.mean, fitness.max, 
# fitness.mean.silhouette, fitness.max.silhouette". 
######################################################
                   
allfitness <- function(population, vfitness, 
                       fitness.f = "fitness.mean")
{
 if (any(myenv$MODIFI == 1)){
  pos <- which(myenv$MODIFI == 1)
  vfitness[pos] <- sapply(population[pos], fitness.f)
 } 
 return (vfitness)
}

#-----------------------------------------------------------------------

##################################################
# "remove.rep.elemt" function.
# Removes repeated elements in diferent clusters.
##################################################
                               
remove.rep.elemt <- function(clustering, level)
{
 pos <- 1; 
 c <- myenv$CHROM_LENGTH - level + 2 # Length of clustering.
 temp <- k <- 0; v <- i <- integer(0)
 
 for (elem in 1:myenv$OBJ.COUNT)
 {
  i <- 0
  while (pos <= c)
  {
   cl <- clustering[[pos]]
   if (elem %in% cl) # if "elem" belongs to "cluster"...
    { 
     if (i > 0) # Removes repeated elements of worse cluster...
      {#removes "elem".
       k1 <- clustering[[i]][-match(elem, clustering[[i]])]
       k2 <- cl[-match(elem, cl)]
       if (is.na(cl[2])) 
        clustering[[i]] <- k1
       else if (is.na(clustering[[i]][2]))
        clustering[[pos]] <- k2
       else
        {
         temp <- cluster.dist(c(elem), k1) < cluster.dist(c(elem), k2)
         if (temp) clustering[[pos]] <- k2  
          else clustering[[i]] <- k1
        }  
       break() 
      }
     i <- pos
    }
   pos <- pos + 1  
  }
  # "v" will be missing elements vector.
  if (i == 0) v[k <- k + 1] <- elem
  pos <- 1
 }
 return (list(v, clustering))
}

#######################################################
# "long.cluster" function.
# Returns position of an empty cluster when "long=0",
# in other case (long>0), it returns position of first 
# cluster with size more than "long".
#######################################################

long.cluster <- function(clustering, level, long = 0)
{
 i <- 1; c <- myenv$CHROM_LENGTH - level + 2 # Length of clustering.
 if (long > 0)  
   while ( i <= c)
    {# ¿Are there clusters with size more than "long"? 
     if (!is.na(clustering[[i]][long + 1])) break()
     i <- i + 1
    }
  else  
   {# ¿Are there empty clusters?
    i <- match(list(integer(0)), clustering)
    if (is.na(i)) i <- c + 1
   } 
 if (i > c) return (0)
  else return (i)
}
                                    
#################################################
# "choice.elem" function.
# Choices randomly a number of "set", 
# removes choiced element.
#################################################

choice.elem <- function(set)
{
 pos <- random.int(1:length(set), 1)
 nomb <- set[pos]
 set <- set[-pos]
 return (list(nomb, set)) 
}

#####################################################
# "worse.cluster" function.
# Returns position of worse cluster, so that, the 
# cluster is differents of either empty or one size.
#####################################################

worse.cluster <- function(clustering, level, homog.f = "cluster.fitness")              
{
 c <- myenv$CHROM_LENGTH - level + 2 # Length of clustering.
 v <- sapply(clustering, homog.f) 
 elem <- v[pos <- long.cluster(clustering, level, 1)]
 for (i in pos:c)
  # If cluster is not either empty or has size one...
  if (!is.na(clustering[[i]][1]) && !is.na(clustering[[i]][2]))
   if (v[i] > elem) {elem <- v[i]; pos <- i}# Taking worse cluster position.             
 return (pos) 
}

#########################################################
# "min.dist" function.
# Returns the cluster position with smallest distance 
# from "elem" to "clustering".
#########################################################

min.dist <- function(elem, clustering, level)
{
 i <- pos <- 1; minimum <- myenv$MAX_DIST + 1
 len <- myenv$CHROM_LENGTH - level + 2 # Length of clustering. 
 for (i in 1:len)
  if ((j <- cluster.dist(c(elem), clustering[[i]])) < minimum)
   { # Take the smallest distance... 
    minimum <- j
    pos <- i
   } 
 return (pos) 
}

########################################################
# "repair.cluster" function.
# Repear a clustering when it's doing the crossover,
# returns a partition, so that, it's in "clustering".
########################################################

repair.cluster <- function(clustering, level, 
                            homog.f = "cluster.fitness")
{
 v <- remove.rep.elemt(clustering, level) # Removes repeated elements...
 clustering <- v[[2]]; v <- v[[1]] # Gets the both  new clustering and
                                   # missing elements. 
 # if v is not empty...
 while (!is.na(v[1]) &&  
        ((pos <- long.cluster(clustering, level)) > 0)) # Are there 
  {                                                     # empty sets?
   l <- choice.elem(v); v <- l[[2]]
   clustering[[pos]][1] <- l[[1]]
  }
 while ((pos <- long.cluster(clustering, level)) > 0) # If there are still
  {                                            # empty sets...
   l <- worse.cluster(clustering, level, homog.f)
   i <- random.int(1:length(clustering[[l]]), 1)
   clustering[[pos]][1] <- clustering[[l]][i]
   clustering[[l]] <- clustering[[l]][-i]
  }
 while (!is.na(v[1])) # If there are still elements into v...
  { 
   i <- choice.elem(v); v <- i[[2]]
   pos <- min.dist(i[[1]], clustering, level)  
   clustering[[pos]] <- c(clustering[[pos]], i[[1]])
  }  
 return (clustering)
}

####################################################
# "clusters.select" function.
# Choices parents clusters to form one cluster child.
####################################################

clusters.select <- function(clustering1, clustering2, level, 
                            homog.f = "cluster.fitness")
{
 l <- myenv$CHROM_LENGTH - level + 2 # Length of clustering.
 newclustering <- list(c(0))
 v1 <- sapply(clustering1[-l-1], homog.f)
 v2 <- sapply(clustering2[-l-1], homog.f)
 c <- l %/% 2
 
 for (i in 1:c) # Choicing of clusters to crossover.
  {
   pos <- which.min(v1)
   newclustering[[i]] <- clustering1[[pos]]
   v1[pos] <- myenv$MAX_DIST + 1
   pos <- which.min(v2)
   newclustering[[c + i]] <- clustering2[[pos]]
   v2[pos] <- myenv$MAX_DIST + 1
  }
 if ((l %% 2) > 0) # If "l" is odd, then, it fails a cluster.
  {
   if (random.int(1:2, 1) == 1) 
    newclustering[[c + i + 1]] <- clustering1[[which.min(v1)]]
   else
    newclustering[[c + i + 1]] <- clustering2[[which.min(v2)]]
  } 
  
 return (newclustering)
}

####################################################
# "div.cluster" function.
# Splits randomy a cluster in two clusters,
# the "cluster" size most be more than one.
####################################################

div.cluster <- function(cluster1)
{  
 c <- length(cluster1); l <- list(0,0)
 pos <- random.int(1:c, 1)
 if (pos == 1) {l[[1]] <- cluster1[1]; l[[2]] <- cluster1[2:c]}
  else if (pos == c) {l[[1]] <- cluster1[1:(c - 1)]; l[[2]] <- cluster1[c]}
    else {l[[1]] <- cluster1[1:(pos - 1)]; l[[2]] <- cluster1[pos:c]}
    
 return (l) 
}

####################################################
# "div.clustering" function.
# Constructs clustering from level "level.up" 
# toward up levels, it splits clusters in two 
# clusters and so on. 
####################################################

div.clustering <- function(dendogram, level.up)
{
 len <- myenv$CHROM_LENGTH - level.up + 2 # clustering size.
 while (len < (myenv$CHROM_LENGTH + 1))
  {
   pos <- random.int(1:len, 1)
   # while cluster size lets 1 be.
   while (is.na(dendogram[[level.up]][[pos]][2])) 
    if (pos == len) pos <- 1 # if "pos" is the last one position.
     else pos <- pos + 1
   l <- div.cluster(dendogram[[level.up]][[pos]])
   dendogram[[level.up - 1]] <- dendogram[[level.up]]   
   dendogram[[level.up - 1]][[pos]] <- l[[1]]
   dendogram[[level.up - 1]][[len <- len + 1]] <- l[[2]]
   dendogram[[level.up]][[len]] <- c(pos, len)
   level.up <- level.up - 1 
  }
 dendogram[[level.up]][[len + 1]] <- c(0, 0) 
 
 return (dendogram)
}

####################################################
# "div.clustering" function.
# Constructs clustering from level "level.up" 
# toward up levels, it splits clusters in two 
# clusters and so on. Implemented by "mutation2".
####################################################

div.clustering2 <- function(dendogram, level.up, homog.f = "cluster.fitness")
{
 len <- myenv$CHROM_LENGTH - level.up + 2 # clustering size.
 while (len < (myenv$CHROM_LENGTH + 1))
  {
   #pos <- random.int(1:len, 1)
   # Searching position of cluster less compact
   # for split.
   #pos <- which.max(sapply(dendogram[[level.up]], 
   #                 "H.ave"));
   
   pos <- which.max(sapply(dendogram[[level.up]], homog.f))
   # while cluster size be 1...
   while (is.na(dendogram[[level.up]][[pos]][2])) 
    if (pos == len) pos <- 1 # if "pos" is the last one position.
     else pos <- pos + 1
   l <- div.cluster(dendogram[[level.up]][[pos]])
   dendogram[[level.up - 1]] <- dendogram[[level.up]]   
   dendogram[[level.up - 1]][[pos]] <- l[[1]]
   dendogram[[level.up - 1]][[len <- len + 1]] <- l[[2]]
   dendogram[[level.up]][[len]] <- c(pos, len)
   level.up <- level.up - 1 
  }
 dendogram[[level.up]][[len + 1]] <- c(0, 0) 
 return (dendogram)
}


####################################################
# "crossover" function.
# It carries out a crossover between two dendograms.
####################################################

crossover <- function(dendogram1, dendogram2, 
                      homog.fit = "cluster.fitness") 
{
 level <- random.int(1:myenv$CHROM_LENGTH, 1)
 cl <- clusters.select(dendogram1[[level]], 
                       dendogram2[[level]], level)
 cl <- repair.cluster(cl, level, homog.fit)
 dendogram1[[level]] <- cl
 dendogram1 <- div.clustering(dendogram1, level)
 if (level < myenv$CHROM_LENGTH) 
  dendogram1 <- mutation(dendogram1, level, homog.f = homog.fit)
  
 return (dendogram1)
}

####################################################
# "crossover2" function.
# It carries out a crossover between two dendograms
# with mutation2.
####################################################

crossover2 <- function(dendogram1, dendogram2, homog.fit = "cluster.fitness") 
{   
 level <- random.int(1:myenv$CHROM_LENGTH, 1)
 cl <- clusters.select(dendogram1[[level]], 
                       dendogram2[[level]], level)
 cl <- repair.cluster(cl, level, homog.fit)
 dendogram1[[level]] <- cl
 dendogram1 <- div.clustering2(dendogram1, level, homog.f = homog.fit)
 
 if (level < myenv$CHROM_LENGTH)
  dendogram1 <- mutation2(dendogram1, level, homog.f = homog.fit) 
  #dendogram1 <- do.call(mutation.op, list(dendogram1, level))
   
 return (dendogram1)
}

####################################################
# "crossover3" function.
# It carries out a crossover between two dendograms
# with mutation2. "exec.strag" is used to improve the
# new clustering.
####################################################

crossover3 <- function(dendogram1, dendogram2, homog.fit = "cluster.fitness") 
{   
 level <- random.int(1:myenv$CHROM_LENGTH, 1)
 len <- myenv$CHROM_LENGTH - level + 2
 cl <- clusters.select(dendogram1[[level]], 
                       dendogram2[[level]], level)
 cl <- repair.cluster(cl, level, homog.fit)
 cl <- evol.cluster(c(cl, list(c(0, 0))), myenv$ITE.CROSS %/% 2) # Improving the clustering.
 cl <- evol.cluster(cl, myenv$ITE.CROSS %/% 2, mop = "gene.change2") # Improving the clustering.
 dendogram1[[level]] <- cl[-(len + 1)]
 dendogram1 <- div.clustering2(dendogram1, level, homog.f = homog.fit)
 
 if (level < myenv$CHROM_LENGTH)
   dendogram1 <- mutation2(dendogram1, level, homog.f = homog.fit)
  #dendogram1 <- do.call(mutation.op, list(dendogram1, level))
   
 return (dendogram1)
}

####################################################
# "selection" function.
# This is the selection phase, rulette wheel 
# method is applied to select the individuals...
####################################################

selection <- function(population, fitness.values)
{
 modify <- myenv$MODIFI
 add <- sum(fitness.values)
 rulette <- fitness.values / add # Relative probability vector.
 rulette <- round(cumsum(rulette), 2) # Rulette wheel method.
 # Vector of random numbers betwen 0 and 1, 
 # it's for selection of the most fit individuals.
 flip <- round(runif(myenv$POPULATION_SIZE), 2) 
 k <- 0; newpopulation <- list(rep(0, myenv$POPULATION_SIZE))
 #fitnessvector <- rep(0, myenv$POPULATION_SIZE)
 for (i in 1:myenv$POPULATION_SIZE)
  {
   j <- 1
   while (flip[i] > rulette[j]){ 
    j <- j + 1
   } 
   newpopulation[[k <- k + 1]] <- population[[j]] # Update new population.
   modify[k] <- fitness.values[j]
  }
 myenv$MODIFI <- modify
 return (newpopulation)
}

####################################################
# "mutation.all" function.
# dendogram mutation with probability <= PMUT. 
####################################################
                      
mutation.all <- function(population, mutate = "mutation2", comb1 = F,
                         homog = "cluster.fitness")                        
{
 modify <- myenv$MODIFI
 # Vector of mutation probabilities.                                
 flip <- round(runif(myenv$POPULATION_SIZE), 1) #3 
 for (i in 1:myenv$POPULATION_SIZE)
  if (flip[i] <= myenv$PMUT) 
  {
   #print(i)
   population[[i]] <- do.call(mutate, list(population[[i]], comb = comb1,
                                           homog.f = homog))
   #population[[i]] <- mutation2(population[[i]], comb = comb1)
   modify[i] <- 1
  } 
 myenv$MODIFI <- modify   
 return (population) 
}

####################################################
# "mutation.all2" function.
# dendrogram mutation without using probability values,
# it mutates each individual.  
####################################################

mutation.all2 <- function(population, vfitness, mutate = "mutation2",
                          fitness.f = "fitness.mean", 
                          homog = "cluster.fitness", comb1 = F)                      
{
 pop <- lapply(population, mutate, comb = comb1, homog.f = homog)
 
 vect <- sapply(pop, fitness.f, USE.NAMES = F)
 pos <- which(vect > vfitness)
 population[pos] <- pop[pos]
 vfitness[pos] <- vect[pos] 
 
 return (list(population, vfitness)) 
}

####################################################
# "crossover.all" function.
# dendograms crossover with probability <= PCROSS. 
####################################################
                 
crossover.all <- function(population, fitness.f = "fitness.mean",
                          homog = "cluster.fitness", 
                          cross.op = "crossover2")
{
 modify <- myenv$MODIFI
 k <- pos <- 0
 # Vector of crossover probabilities.                                
 flip <- round(runif(myenv$POPULATION_SIZE), 1)
 for (i in 1:myenv$POPULATION_SIZE)
  if (flip[i] <= myenv$PCROSS){ 
    k <- k + 1
    if (k == 2){ # if it founds two dendograms to crossover...
     temp <- do.call(fitness.f, list(population[[pos]])) < 
             do.call(fitness.f, list(population[[i]]))
     if (temp)
      {
       population[[i]] <- do.call(cross.op, list(population[[pos]], 
                                                 population[[i]], 
                                                 homog.fit = homog)) 
       #population[[i]] <- crossover2(population[[pos]], 
       #                              population[[i]])  
       modify[i] <- 1
      } 
     else
      {
       population[[pos]] <- do.call(cross.op, list(population[[pos]], 
                                                   population[[i]], 
                                                   homog.fit = homog))
       #population[[pos]] <- crossover2(population[[pos]], 
       #                                population[[i]])    
       modify[pos] <- 1
      }
     k <- pos <- 0
    }
    else pos <- i 
  }
     
 if (k == 1){ # if there exists a dendogram without a crossover... 
  i <- random.int(1:myenv$POPULATION_SIZE, 1)
  if (i == pos){
   if (i < myenv$POPULATION_SIZE) i <- i + 1
    else i <- 1
  }
  population[[pos]] <- do.call(cross.op, list(population[[pos]], 
                                              population[[i]], 
                                              homog.fit = homog))
  #population[[pos]] <- crossover2(population[[pos]], 
  #                                population[[i]])
  modify[pos] <- 1
 } 
 myenv$MODIFI <- modify
 return (population) 
}                                      
              
####################################################
# "crossover.all2" function.
# dendogram crossover with probability <= PCROSS. 
####################################################

crossover.all2 <- function(population, fitness.f = "fitness.mean")
{
 modify <- myenv$MODIFI
 k <- pos <- 0
 # Vector of crossover probabilities.                                
 flip <- round(runif(myenv$POPULATION_SIZE), 1)
 for (i in 1:myenv$POPULATION_SIZE)
  if (flip[i] <= myenv$PCROSS){ 
    k <- k + 1
    if (k == 2){ # if it found two dendograms to crossover...
     temp <- do.call(fitness.f, list(population[[pos]])) < 
             do.call(fitness.f, list(population[[i]]))
     if (temp)
      {
       population[[i]] <- crossover3(population[[pos]], 
                                     population[[i]])  
       modify[i] <- 1
      } 
     else
      {
       population[[pos]] <- crossover3(population[[pos]], 
                                       population[[i]])    
       modify[pos] <- 1
      }
     k <- pos <- 0
    }
    else pos <- i 
  }  
   
 if (k == 1){ # if there exist a dendogram without a crossover... 
  i <- random.int(1:myenv$POPULATION_SIZE, 1)
  if (i == pos){
   if (i < myenv$POPULATION_SIZE) i <- i + 1
    else i <- 1
  }
  population[[pos]] <- crossover3(population[[pos]], 
                                  population[[i]])
  modify[pos] <- 1
 } 
 myenv$MODIFI <- modify
 return (population) 
}                                      

##########################################
# "reordering.dendro" reorders the clusters  
# of a dendrogram according to the order 
# of the objects in the last level of it.
# This is useful when we graphic
# dendrograms with function "plot".
##########################################

reordering.dendro <- function(dendro){
  # Internal function
  creating.newdendro <- function(i){
    pos <- dendro[[i]][[size]][1]
    cluster1 <- dendro[[i-1]][[pos]]
    pos <- searchcluster(cluster1, dendro.result[[i-1]], size)  
    clustering <- dendro.result[[i-1]] 
    cluster1 <- c(clustering[[pos]], clustering[[pos+1]])
    # Creating the next upper level
    if (pos > 1)
      clustering <- c(clustering[1:(pos-1)], list(cluster1), 
                      clustering[(pos+2):(size+1)])
    else
      clustering <- c(list(cluster1), clustering[3:(size+1)])
    
    clustering[[size]] <- c(pos, pos+1)  
    dendro.result[[i]] <<- clustering
    size <<- size - 1         
  }
  # Internal function
  searchcluster <- function(cluster1, clustering, long){
    subcluster <- function(cluster1, cluster2)
      all(cluster1 %in% cluster2) && all(cluster2 %in% cluster1)
 
    i <- 1                                
    # Searching a cluster in a clustering.
    while (i <= long){
      if (subcluster(cluster1, clustering[[i]])) break()
      i <- i + 1
    } 
    if (i <= long) return (i)
    return (0)            
  }

  len <- length(dendro)
  ordered.objs <- c(dendro[[len]][[1]], dendro[[len]][[2]])
  setsize <- length(ordered.objs)
  dendro.result <- as.list(1:len)
  dendro.result[[1]] <- create.firstlevel(dendro[[1]], ordered.objs, 
                                          setsize)                                                                  
  size <- length(dendro[[1]]) - 1
  sapply(2:len, creating.newdendro, USE.NAMES = F)
  return (dendro.result)
}

#################################################
# Auxiliary function of "reordering.dendro" that
# builts the first clustering the new dendrogram
# rearranged by "reordering.dendro".
#################################################

create.firstlevel <- function(firstclustering, ordered.objs, setsize){
  # Internal function
  search.firstobj <- function(cluster1)
    match(cluster1[1], ordered.objs)

  # Internal function
  search.lessnumb <- function(){
    findless <- function(i){
      if (posvector[i] < less){
         less <<- posvector[i]
         pos <<- i
      }      
    }
    less <- posvector[1]; pos <- 1
    sapply(2:len, findless, USE.NAMES = F)
    posvector[pos] <<- Inf
    return (c(less, pos)) 
  }
  # Internal function
  creatingthelevel <- function(i){
    pos <- search.lessnumb()
    clust <- ordered.objs[pos[1]:(pos[1] + lengvector[pos[2]] - 1)]
    clustresult <<- c(clustresult, list(clust))  
  }
  
  len <- length(firstclustering) - 1
  clustering <- firstclustering[1:len] 
  posvector <- sapply(clustering, search.firstobj, USE.NAMES = F)
  lengvector <- sapply(clustering, length)
  clustresult <- list()
  sapply(1:len, creatingthelevel, USE.NAMES = F)
  
  return (c(clustresult, list(c(0,0))))        
}

####################################################
# "complete.tree" function.
# It randomly builts the first half of a dendogram. 
####################################################

complete.tree <- function(dendogram)
{
 # reordering the dendrogram.
 dendogram <- reordering.dendro(dendogram)
 # building the fist part of the dendrogram.
 level.up <- trunc(myenv$OBJ.COUNT * myenv$PART) + 1
 len <- myenv$OBJ.COUNT - level.up # clustering size.
 dendogram <- c(as.list(rep(NA, level.up - 1)), dendogram)
 while (len < (myenv$OBJ.COUNT - 1))
  {
   pos <- random.int(1:len, 1)
   # while cluster size lets 1 be.
   while (is.na(dendogram[[level.up]][[pos]][2])) 
    if (pos == len) pos <- 1 # if "pos" is the last one position.
     else pos <- pos + 1
   l <- div.cluster(dendogram[[level.up]][[pos]])
   dendogram[[level.up - 1]] <- dendogram[[level.up]]   
   dendogram[[level.up - 1]][[pos]] <- l[[1]]
   dendogram[[level.up - 1]][[len <- len + 1]] <- l[[2]]
   dendogram[[level.up]][[len]] <- c(pos, len)
   level.up <- level.up - 1 
  }
 dendogram[[level.up]][[len + 1]] <- c(0, 0) 
 return (dendogram)
}

####################################################
# "replacebyname" function.
# Changes object numbers by their real names. 
####################################################

replacebyname <- function(clustering, dimname)
{
 for (i in 1:(length(clustering) - 1))
  {
   cluster1 <- character(0)
   for (j in 1:length(clustering[[i]]))
    cluster1[j] <- dimname[clustering[[i]][j]]
   clustering[[i]] <- cluster1 
  }  
  
 return (clustering)  
}

#--------------------------------------------------------------------

####################################################
# "ag.coef" function.
# Computing A.C, it is the clustering with best fitness
# from the half to end of dendogram, considering all the
# levels of the dendrogram.
####################################################

ag.coef <- function(dendogram, 
                    fitness.clust = "clustering.fitness")
{
 part <- trunc(myenv$OBJ.COUNT * myenv$PART)
 return (which.max(sapply(dendogram[-(1:part)],
                   fitness.clust)) + part)
}                  

ag.coef1 <- function(dendogram, 
                    fitness.clust = "clustering.fitness")
  which.max(sapply(dendogram, fitness.clust))
                  

####################################################
# "ag.coef.mean" function.
# Computing A.C, it is clustering with best fitness
# from half to end of dendogram. Considering all 
# levels of the dendrogram.
####################################################

ag.coef.mean <- function(dendogram)
{
 part <- trunc(myenv$OBJ.COUNT * myenv$PART)
 return (which.max(fitness.list(dendogram[-(1:part)])) + part)
}         

####################################################
# "ag.coef.mean1" function.
# Computing A.C, it is clustering with best fitness
# from half to end of dendogram. It is not Considers 
# all levels of the dendrogram, just the meaning
# levels.
####################################################

ag.coef.mean1 <- function(dendogram)
 which.max(fitness.list(dendogram))

####################################################
# "ag.coef.sigma" function.
# Computing A.C, it is clustering with best fitness
# from half to end of dendogram. It is not Considers 
# all levels of the dendrogram, just the meaning
# levels.
####################################################

ag.coef.sigma <- function(dendogram)
 which.max(fitness.list(dendogram) + sigma.list(dendogram))
    
#---------------------------------------------------------------------    
          
####################################################
# "subcluster" function.
# It gives TRUE if cluster1 is subset of cluster2.
####################################################

subcluster <- function(cluster1, cluster2)
 all(cluster1 %in% cluster2)  
 
####################################################
# "subset.cluster" function.
# It gives TRUE if cluster is into clustering.
####################################################

subset.cluster <- function(cluster1, clustering, long)
{
 i <- 1                                
 # Searching "cluster" into "clustering".
 while (i <= long){
  if (subcluster(cluster1, clustering[[i]])) break()
  i <- i + 1
 } 
 return (i <= long)         
}

####################################################
# "merged.cluster" function.
# Gives a clusters list with all clusters merged 
# into dendogram.
####################################################

merged.cluster <- function(dendogram)
{
 k <- 0; clust <- clustering
 tempclust <- as.list(1:myenv$OBJ.COUNT)
 for (i in 0:(myenv$OBJ.COUNT - 3))
  {
   long <- myenv$OBJ.COUNT - i; j <- 1
   while (j <= (long - 1))
    if (!subset.cluster(dendogram[[i + 1]][[j]], tempclust, long))
     { # Gets the new cluster merged in this level.
      clust[[k <- k + 1]] <- dendogram[[i + 1]][[j]]
      break()
     }
    else j <- j + 1
   tempclust <- dendogram[[i + 1]]   
  }
 clust[[k + 1]] <- c(dendogram[[myenv$OBJ.COUNT - 2]][[1]],
                     dendogram[[myenv$OBJ.COUNT - 2]][[2]])
 return (clust)
}

###########################################################
# Compute the distances where two clusters are merged for
# each level of the dendrogram.
############################################################
 
union.distances <- function(dendro, homog.f = "cluster.fitness")
{
 create.dist <- function(i){
   posmerge <- dendro[[i]][[numb.clusts]]
   clust <- c(dendro[[i - 1]][[posmerge[1]]], dendro[[i - 1]][[posmerge[2]]])
   #print(numb.clusts);print(clust)
   distances <<- c(distances, do.call(homog.f, list(clust)))
   numb.clusts <<- numb.clusts - 1
 }
 
 numb.clusts <- length(dendro[[2]])
 distances <- c()            
 
 sapply(2:length(dendro), create.dist, USE.NAMES = FALSE)
  
 distances <- sort(distances) 
 distances <- c(distances[1] / 4, distances[1] / 2, distances)
 
 return(distances)
}
 
####################################################
# "height.cluster" function.
# Gives a vector with distances between merging
# clusters at the successive stages.
####################################################

height.cluster <- function(list.cluster, distances)
{
 permut <- list.cluster[[myenv$OBJ.COUNT - 1]]
 vect <- vect1 <- integer(0)
 
 for (i in 1:myenv$OBJ.COUNT) # Inverse permutation.
  vect1[permut[i]] <- i
 k <- 0
 for (i in 1:(myenv$OBJ.COUNT - 1))
  {
   j <- 1; c <- length(list.cluster[[i]]) - 1
   while (j <= c)
    { # Taking position of each object en "pos".
     pos <- vect1[list.cluster[[i]][j]]
     if (is.na(vect[pos]))
      { # Assigning incremental values to each
        # merged cluster into dendogram.
       vect[pos] <- distances[k <- k + 1] 
       #vect[pos] <- (k <- k + 1)
       break()
      }
     j <- j + 1
    }
  }
 return (vect)
}
####################################################
# "merge.matrix" function.
# Given a matrix, row i of merge describes the  
# merging of clusters at step i of clustering.
####################################################

merge.matrix <- function(list.cluster)
{
 tmatrix <- matrix(nrow = (myenv$OBJ.COUNT - 1), ncol = 2, byrow = TRUE)
 tmatrix[1, ] <- -1 * list.cluster[[1]]
 j <- 1
 for (k in 2:(myenv$OBJ.COUNT - 1))
  { # Searching the subset clusters of cluster in 
    # position k.
   i <- 1; j <- 0; pos <- integer(0)
   while (i < k)
    {
     if (!any(is.na(list.cluster[[i]])))
      if (subcluster(list.cluster[[i]], list.cluster[[k]]))
       {
        j <- j + 1
        pos[j] <- i
        if (j == 2) break() 
       }
     i <- i + 1  
    }
   # filling the matrix with either objects or levels numbers.   
   if (j == 0) tmatrix[k, ] <- -1 * list.cluster[[k]]
    else if (j == 2) 
     {
      tmatrix[k, ] <- pos
      list.cluster[[pos[1]]] <- NA
      list.cluster[[pos[2]]] <- NA
     } 
   else # when j == 1 ...
    {
     if (list.cluster[[pos]][1] == list.cluster[[k]][1])
      {
       tmatrix[k, 2] <- list.cluster[[k]][length(list.cluster[[k]])]
       tmatrix[k, 2] <- -1 * tmatrix[k, 2]
       tmatrix[k, 1] <- pos
      }
     else 
      {
       tmatrix[k, 1] <- -1 * list.cluster[[k]][1]
       tmatrix[k, 2] <- pos
      } 
     list.cluster[[pos]] <- NA 
    }
  }
 return (tmatrix) 
}

####################################################
# "dendrogram.graph" function.
# Given a full dendrogram on a dataset, it builds  
# a object of class "agnes", which can display the 
# graphic of the dendrogram by using R function "plot".
####################################################

dendrogram.graph <- function(dendogram, datam, 
                             fitness.clust = "clustering.fitness",
                             homog.f = "cluster.fitness")
{
 sort.dimnames <- function(i)
    names.obj[i] <<- temp[graph$order[i]]
    
 # vector with the distances where are merged two clusters in each level.
 distances <- union.distances(dendogram, homog.f)
 
 mc <- merged.cluster(dendogram)
 
 graph <- list(order = c(dendogram[[myenv$OBJ.COUNT - 2]][[1]],
                         dendogram[[myenv$OBJ.COUNT - 2]][[2]]))
                         
 names.obj <- character(0); temp <- (dimnames(myenv$DISTANCE))[[1]]
 
 # Ordering the object names (dimnames).
 sapply(1:myenv$OBJ.COUNT, sort.dimnames, USE.NAMES = FALSE)
 #for (i in 1:myenv$OBJ.COUNT)
 # names.obj[i] <- temp[graph$order[i]]
  
 # Computing the A.C, that is, the clustering of best fitness.
 a.c_dist <- distances[ag.coef(dendogram, fitness.clust) - 1]
 
 # Constructing the "graph" list.
 graph <- c(graph, list(order.lab = names.obj,              
                        height =  height.cluster(mc, distances), 
                        ac = round(a.c_dist / distances[myenv$OBJ.COUNT - 1], 2),
                        #ac = round(a.c / (myenv$OBJ.COUNT - 1), 2),
                        merge = merge.matrix(mc), 
                        diss = as.dist(myenv$DISTANCE),
                        data = datam))
                        
 class(graph) <- c("agnes", "twins")
 
 return (graph)
}             

#-------------------------------------------------------------------
                 
####################################################
# "init.variable" function.
# Initials the global variables...   
####################################################

reduce.tree <- function(dmatrix, popsize, endpart,
                        filename = "tree.Rdata", 
                        fitness = "fitness.mean")
{
 # load saved tree.
 load(filename)
 init.variables(dmatrix, pop.size=popsize, gen.numb=5, part=endpart) 
 temp <- length(elite)
 temp <- temp - myenv$CHROM_LENGTH
 elite <- elite[-(1:temp)]
#print(elite.fitness) 
 elite.fitness <- do.call(fitness, list(elite))
 for (i in 1:popsize)
  {
   oldp[[i]] <- oldp[[i]][-(1:temp)]
   vector.fitness[i] <- do.call(fitness, list(oldp[[i]]))
  } 
 save(elite, elite.fitness, oldp, vector.fitness, 
      file = filename)
 #print(length(elite))
 
 return (elite.fitness)     
}

                                              
local.search <- function(dist.matrix, pop.size = 10, gen.numb, 
                      part = 1/2, mutation = "mutation2",
                      fitness.f = "fitness.mean",
                      homog.f = "cluster.fitness", 
                      return.pop = FALSE, init.pop = NULL)
{# Initializing constants and variables...

 init.varsearch(dist.matrix, pop.size, gen.numb, part)
 
 oldp <- newp <- population 
 
 myenv$MODIFI <- rep(1, myenv$POPULATION_SIZE)
 vector.fitness <- integer(myenv$POPULATION_SIZE)
 # Initial population.
 if (is.null(init.pop)) oldp <- init.population(oldp)
  else oldp <- init.pop 
 # Gets a fitness vector of all individuals.
 vector.fitness <- allfitness(oldp, vector.fitness, fitness.f) 
 
 # Starts the local evolutionary strategy.
 
 t <- 0
 
 start=Sys.time() # get initial time...
 
 while (t < myenv$MAX_GEN)            
 {  
  t <- t + 1  # Time increment.
 
  temp <- mutation.all2(oldp, vector.fitness, mutation, fitness.f, 
                        comb1 = T, homog = homog.f)
  # Get the both dendogram and vector fitness.
  newp <- temp[[1]]; vector.fitness <- temp[[2]]
 
  oldp <- newp # Update old population with the current population.
  rm(newp, temp)
  #print(round(vector.fitness, 2))
 }     
 
 temp <- which.max(vector.fitness)
 
 elite <- oldp[[temp]]
 elite.fitness <- vector.fitness[temp]
 
 #if (savefile) save(elite, elite.fitness, oldp, vector.fitness, 
 #                   file = filename)
 
 #elite <- complete.tree(elite) 
 
 #if (read.file) print(round(vector.fitness, 2))
 
 print(Sys.time()-start)
 flush.console()     # get finish time...
 print(paste("Elite fitness: ", as.character(elite.fitness)))
 
 if (!return.pop) return (elite)
 
 return( list(elite = elite, lastpop = oldp) )
}
 
 
gene.change1 <- function(clustering, size, minimal = 1)
{
 clust <- clustering
 rand <- random.int(1:size, 2)
 leng <- length(clust[[rand[1]]])
 if (leng <= minimal) 
  {
   rand <- rev(rand)
   leng <- length(clust[[rand[1]]])
  } 
 if (leng > minimal)
  {
   pos <- random.int(1:leng, 1)
   elem <- clust[[rand[1]]][pos]
   clust[[rand[1]]] <- clust[[rand[1]]][-pos]
   clust[[rand[2]]] <- c(clust[[rand[2]]], elem)    
  }
 return (clust) 
}     

gene.change2 <- function(clustering, size)
{
 clust <- clustering
 rand <- random.int(1:size, 2)
 leng1 <- length(clust[[rand[1]]])
 leng2 <- length(clust[[rand[2]]])
 pos1 <- random.int(1:leng1, 1)
 pos2 <- random.int(1:leng2, 1)
 elem <- clust[[rand[1]]][pos1]
 clust[[rand[1]]][pos1] <- clust[[rand[2]]][pos2]
 clust[[rand[2]]][pos2] <- elem
 return (clust) 
}     

evol.cluster <- function(clustering, gen.numb, 
                         fitness.f = "clustering.fitness",
                         mop = "gene.change1",
                         print.fitness = FALSE)
{
 t <- 1; count <- 0; clust <- clustering 
 leng <- length(clustering) - 1
 ff <- do.call(fitness.f, list(clust))
 
 if (print.fitness){
   print(paste("Init_Fitness: ", as.character(ff)))
   start=Sys.time() # get initial time...
 }  
   
 while (t <= gen.numb)
  {
   t <- t + 1
   clust1 <- do.call(mop, list(clust, leng))
   ff1 <- do.call(fitness.f, list(clust1))
   if (ff1 > ff){  
    count <- count + 1
    ff <- ff1
    clust <- clust1
   }
  } 
  
 if (print.fitness){
   print(paste("End_Fitness: ", as.character(ff)))
   print(Sys.time()-start)
   flush.console()     # get finish time...
 }
    
 return (clust) 
} 

####################################################
# "exec.strag" function.
# It improves a clustering by swapping and moving 
# inter-cluster elements.
####################################################

exec.strag <- function(clustering, times = 1000, count = 10, 
                       fitness.f = "clustering.fitness",
                       print.fitness = FALSE)
{
 #start=Sys.time() # get initial time...
 repeating <- function(clustering, times){
  clust <<- evol.cluster(clustering, times, fitness.f, 
                         mop = "gene.change1")
  clust <<- evol.cluster(clust, times, fitness.f, 
                         mop = "gene.change2")
 }
 
 clust <- clustering
 
 if (print.fitness){
   ff <- do.call(fitness.f, list(clust))
   print(paste("Init_Fitness: ", as.character(ff)))
   start=Sys.time() # get initial time...
 }
 
 replicate(count, repeating(clust, times))
 
 if (print.fitness){
   ff <- do.call(fitness.f, list(clust))
   print(paste("End_Fitness: ", as.character(ff)))
   print(Sys.time()-start)
   flush.console()     # get finish time...
 }
 
 #print(Sys.time()-start)
 #flush.console()     # get finish time...
 
 return (clust)           
}           

####################################################
# "clustering.dendo" function.
# It improves a clustering and builts a dendrogram
# from the improved clustering.
####################################################

clustering.dendo <- function(clustering, level, times = 500, count = 10, 
                             fitness.f = "clustering.fitness",
                             homog.f = "cluster.fitness",
                             print.fitness = FALSE) 
{
 cl <- clustering; len <- myenv$CHROM_LENGTH - level + 2
 
 if (print.fitness){
   ff <- do.call(fitness.f, list(cl))
   print(paste("Clustering Init_Fitness: ", as.character(ff)))
   start=Sys.time() # get initial time...
 }
 
 cl <- exec.strag(cl, times, count, fitness.f)
 
 if (print.fitness){
   ff <- do.call(fitness.f, list(cl))
   print(paste("Clustering End_Fitness: ", as.character(ff)))
 }
 
 dendo <- as.list(rep(NA, myenv$CHROM_LENGTH))
 dendo[[level]] <- cl[-(myenv$CHROM_LENGTH - level + 3)]
 dendo <- div.clustering2(dendo, level, homog.f)
 
 if (level < myenv$CHROM_LENGTH){ 
  pos <- random.int(1:len, 2)
  dendo[[level + 1]] <- list() 
  dendo[[level + 1]][[len]] <- pos
  dendo <- mutation2(dendo, level, comb = T, homog.f)
 } 
 
 if (print.fitness){
   print(Sys.time()-start)
   flush.console()     # get finish time...
 }
 
 return (dendo)
}


improve.dendo <- function(dendo, ite = 10, operator = "mutation2", 
                          fitness.f = "fitness.mean",
                          homog.f = "cluster.fitness", 
                          print.fitness = FALSE)
{
if (print.fitness) start=Sys.time() # get initial time...
 
 t <- 0 
 new.dendo <- dendo 
 
 new.fitness <- do.call(fitness.f, list(new.dendo))
 
 if (print.fitness)
   print(paste("Init_Fitness: ", as.character(new.fitness)))
 
 while (t < ite)
  {
   t <- t + 1
   #dendogram, level = 0, comb = F
   temp.dendo <- do.call(operator, list(new.dendo, 0, comb = T, homog.f))
   temp.fitness <- do.call(fitness.f, list(temp.dendo))
   
   if (temp.fitness > new.fitness){
    new.fitness <- temp.fitness
    new.dendo <- temp.dendo
   } 
  }
  
 if (print.fitness){ 
   print(paste("End_Fitness: ", as.character(new.fitness))) 
   print(Sys.time()-start)
   flush.console()     # get finish time...  
 }
 
 return (new.dendo) 
}

###############################################################
# CreateFiles function, creates three files from the output dendrogram
# returned by one of the hierarchical clustering methods built in R.
# That is, a dendrogram of classes "twins", "hclust", "agnes", etc.
# Then, these files have a format that can be loaded by the
# 3D-VisualCluster tool.
#
# INPUT:
#
# datamat: a data matrix representing the used dataset.
# agn: an dendrogram in R of class twins, the output of a hierarchical
#      clustering method.
# namedataset: name given to the dataset without the extension. This
#              name will be used to identify the three created files.
# part: it is an optional parameter that if used, then it means,
#       the part of the dendrogram that we want to remove (a valor in
#       form of fraction). This is usefull when we do not want to show
#       all dendrogram levels. For example, if we do not want to show
#       the first half of a dendrogram, then "part" have to be assigned
#       to 1/2 (part = 1/2). Hence, the first half of the dendrogram
#       will be removed.
#
# EXAMPLE IN R OF "CreateFiles":
#
# > data(votes.repub)
# > dendrogram <- diana(votes.repub, metric = "manhattan", stand = TRUE)
# > CreateFiles(votes.repub, dendrogram, "VotesRepub_Diana_Manhattan", part = 3/4)
#
# This creates the files "VotesRepub_Diana_Manhattan50x31.txt",
# "VotesRepub_Diana_Manhattan50x31_dendo.txt" and
# "VotesRepub_Diana_Manhattan50x31_dendo_gr.txt" in the current directory.
#
# Note that "part = 3/4" removes the first 3/4 part of the dendrograma,
# that is, since "dendrogram" has 50 levels, the first 39 levels are
# removed (50 * 3/4 + 2 = 39) and just the 11 remaining levels (50 - 39)
# are used.  On the other hand, if we want to explore the whole dendrogram,
# it is enough do not use "part", namely:
# CreateFiles(votes.repub, dendrogram, "VotesRepub_Diana_Manhattan")
###############################################################

CreateFiles <- function(datamat, agn, dataset.name, part = 2){
  numrow <- nrow(datamat); numcol <- ncol(datamat)
  dataset.name <- paste(dataset.name, numrow, "x", numcol, sep="")
  # It creates a file with "datamat" as a dataset, includes column
  # and row heads.
  write.table(datamat, file= paste(dataset.name,".txt", sep=""), sep =";",
              row.names = paste(c("g"), 1:numrow, sep=""),
              col.names = paste(c("c"), 1:numcol, sep=""), quote = F)
  # It transforms "agn" to the 3D-VisualCluster structure.
  dendo <- transf.tree(agn, part)
  #dendo <- pos.mergedcluster(dendo)
  # It prints to a file all levels of the "dendo" dendrogram.
  write(paste(dataset.name,".txt", sep=""),
        file = paste(dataset.name, "_dendo", ".txt", sep=""))
  sink(paste(dataset.name, "_dendo", ".txt", sep=""), append = TRUE)
  print(dendo)
  sink()
  # It prints to a file the heights and merges of the "agn" dendrogram.
  sink(paste(dataset.name, "_dendo_gr", ".txt", sep=""))
  print(agn$height)
  print(agn$merge)
  sink()
  # Print information for the user.
  print("The following files have been created in the current directory: ")
  print(paste(dataset.name,".txt", sep=""))
  print(paste(dataset.name, "_dendo", ".txt", sep=""))
  print(paste(dataset.name, "_dendo_gr", ".txt", sep=""))
}

#############################################################
# The Create.RefPartitionFile function converts a clustering (which 
# is a reference partition) given by a vector structure (internal 
# structure of R)  into a clustering in form of list (3D-VisualCluster 
# structure) and saves the list to a file.
#
# INPUT:
#
# vector.clustering: reference partition (or clustering) as a vector of R.
# filename: name of the file (without extension) where the reference 
#           partition will be saved.
#
# OUTPUT: it saves and returns a clustering in form of list of R. 
#
# EXAMPLE:
#
# > vector.clustering <- c(1,2,1,3,2,3,4,4,5,1,2,5)
# > clust <- Create.RefPartitionFile(vector.clustering, "RefPart_Test")
#
# [1] "File RefPart_Test.txt was created in the current directory."
# [[1]]
# [1]  1  3 10
# [[2]]
# [1]  2  5 11
# [[3]]
# [1] 4 6
# [[4]]
# [1] 7 8
# [[5]]
# [1]  9 12
# [[6]]
# [1] 0 0
#############################################################

Create.RefPartitionFile <- function(vector.clustering, filename){
  extract.cluster <- function(i)
    list.clustering[[i]] <<- which(vector.clustering == i)
    
  cluster.indexes <- unique(vector.clustering)
  list.clustering <- list()
  
  sapply(cluster.indexes, extract.cluster, USE.NAMES = F)
  list.clustering <- c(list.clustering, list(c(0,0)))
  
  # It prints to a file the reference partition in list.clustering.
  sink(paste(filename, ".txt", sep=""))
  print(list.clustering)
  sink()
  
  # Print information for the user.
  print(paste("File ", filename, ".txt was created in ",
              "the current directory.", sep=""))
              
  return(list.clustering)
}


# Replace missing values of a data set by the values computed using
# the technique of k-nearest neigbors. 
missvalue.knn <- function(datamat, K = 20)
{
 missvalue <- function(vect.col, K, x){
  tv <- which(is.na(vect.col)) 
  tk <- K %/% 2
  for (i in tv){
   if ((i-1) >= tk) s1 <- sum(na.omit(vect.col[(i-tk):(i-1)]))
    else s1 <- sum(na.omit(vect.col[1:(i-1)])) + 
               sum(na.omit(vect.col[(x-tk+i):x]))
   if ((x - i) >= tk) s1 <- s1 + sum(na.omit(vect.col[(i+1):(i+tk)]))
    else s1 <- s1 + sum(na.omit(vect.col[(i+1):x])) + 
                    sum(na.omit(vect.col[1:(tk-x+i)]))
   vect.col[i] <- s1/K                 
  }
  return (vect.col) 
 }
 
 x <- nrow(datamat); y <- ncol(datamat)
 for (i in 1:y)
  datamat[,i] <- missvalue(datamat[,i], K, x)
 return (datamat) 
}

# "standard", normalizes a data set with mean 0 and variance 1. 
standard <- function(mat){
 for (i in 1:ncol(mat))
  mat[,i] <- (mat[,i] - mean(mat[,i])) / sqrt(var(mat[,i]))
 return (mat)
}

##################################################
# Module of putting numbers to dendrogram of 
# other methods.
###################################################
 
# Function indentifying the positions of the clusters 
# that are merged to built the next level in the dendrogram.
pos.mergedcluster <- function(dendo)
{
 update.clustering <- function(i){
   k <- 0; j <- 1; sizeclust <<- sizeclust - 1
   while ((j <= sizeclust) && (k < 2)){
     if (!cluster.in.clust(dendo[[i]][[j]], 
                           dendo[[i + 1]], sizeclust - 1))
       {
        k <- k + 1
        dendo[[i + 1]][[sizeclust]][k] <<- j 
       }
     j <- j + 1 
  }  
 }
 
 sizedendo <- length(dendo)
 sizeclust <- length(dendo[[1]]) 
  
 sapply(1:(sizedendo - 1), update.clustering, USE.NAMES = F)
 
 return (dendo)   
}
 
cluster.in.clust <- function(cluster1, clustering, long)
{
 subcluster <- function(cluster1, cluster2)
  all(cluster1 %in% cluster2) && all(cluster2 %in% cluster1)
 
 i <- 1                                
 # Searching "cluster" into "clustering".
 while (i <= long){
  if (subcluster(cluster1, clustering[[i]])) break()
  i <- i + 1
 } 
 return (i <= long)         
}
# Fin del módulo.

agnes.gas <- function(dist.matrix, pop.size = 10, gen.numb = 10, 
                      part = 1/2, fitness.f = "fitness.mean",
                      homog.f = "cluster.fitness",
                      crossover = "crossover2",  
                      mutation = "mutation2", p.cross = 0.40, 
                      p.mut = 0.10, return.pop = FALSE, 
                      init.pop = NULL)
{
 # Initialing constant variables...
 init.variables(dist.matrix, pop.size, gen.numb, p.cross, p.mut, part)
 
 oldp <- newp <- population 
 
 myenv$MODIFI <- rep(1, myenv$POPULATION_SIZE)  
 vector.fitness <- integer(myenv$POPULATION_SIZE)
 # Initial population.
 if (is.null(init.pop))  oldp <- init.population(oldp) 
  else oldp <- init.pop
 # Gets a fitness vector from all individuals.
 vector.fitness <- allfitness(oldp, vector.fitness, fitness.f) 
 # Gets the best individual and its fitness. 
 elite <- oldp[[k <- which.max(vector.fitness)]] 
 elite.fitness <- vector.fitness[k] 
 
 # Starts the genetic algorithm
 t <- 0
 start=Sys.time() # get initial time...
 while (t < myenv$MAX_GEN)            
 {  
 
  t <- t + 1  # Time increment.

  myenv$MODIFI <- integer(pop.size)
  newp <- selection(oldp, vector.fitness)# Selection for mating...

  # "selection" updates myenv$MODIFI with the fitness...
  vector.fitness <- myenv$MODIFI

  myenv$MODIFI <- integer(pop.size)
  rm(oldp)
  
  # Crossover of the individuals with probability <= PCROSS. 
  newp <- crossover.all(newp, fitness.f, homog = homog.f, crossover) 
  
  # Mutation of the individuals with probability <= PMUT.  
  newp <- mutation.all(newp, mutation, comb1 = T, homog = homog.f)  
  
  vector.fitness <- allfitness(newp, vector.fitness, fitness.f)
             
  if (elite.fitness < (k <- max(vector.fitness)))# Updates the elite individual.
   {
    elite <- newp[[which.max(vector.fitness)]] # Best individual.
    elite.fitness <- k
   } 
  else  
   {
    # Elite...
    newp[[1]] <- elite
    vector.fitness[1] <- elite.fitness 
   }

  oldp <- newp # Update old population with the current population.
  rm(newp)#; gc() #- Running the garbage collection
 }

 #if (savefile) save(elite, elite.fitness, oldp, vector.fitness, part, 
 #                   file = filename)
 #elite <- complete.tree(elite) 
 print(Sys.time()-start)
 flush.console()     # get finish time...
 print(paste("Elite fitness: ", as.character(elite.fitness)))

 if (!return.pop) return (elite)
  
 return ( list(elite = elite, lastpop = oldp) ) 
}

#################################################################
# This function saves to a R data file a population of dendrograms
# that would be read later as an initial population by functions
# agnes.gas and local.search.
#################################################################

create.initialpop <- function(population, fitness.f = "fitness.mean",
                              part = 1/2, filename = "Tree.RData"){
 temp <- myenv$MODIFI; lenpop <- length(population)
 myenv$MODIFI <- rep(1, lenpop)
 
 vector.fitness <- integer(lenpop)
 vector.fitness <- allfitness(population, vector.fitness, fitness.f)
 
 elite <- population[[pos <- which.max(vector.fitness)]]
 elite.fitness <- vector.fitness[pos]
 
 oldp <- population
 save(elite, elite.fitness, oldp, vector.fitness, part, file = filename)                              
 
 myenv$MODIFI <- temp
}
                            
#################################################################
# Functions to generate different probality values for 
# the mutation and crossover operators
##################################################################

#p.mutation <- function()
#{
# start=Sys.time() # get initial time...
# LEN <- 4
# LENX <- 20
# LENY <- 20
# pmutation <- seq(length = LENX, from = 0.0, by =0.05)
# #pcrossover <- c(0, seq(length = LENY - 1, from = 0.1, by = 0.045))
# #fitness <- matrix(nrow = LENX, ncol = LENY)
# m.fitness <- matrix(nrow = LEN, ncol = LENX)
# load("cellcycle384x17.RData")
# m1 <- as.matrix(daisy(cellcycle384x17,  metric = "euclidean", stand = TRUE))
# 
# PCROSS <- 0.00 
# for (i in 1:LEN){ k <- 0; print(i)
#  for (j in pmutation)
#   {
#    k <- k + 1
#    PMUT <- j 
#    m.fitness[i, k] <- agnes.gas(m1, 3, 100, part = 12/13)
#   }   
# } 
# dimnames(m.fitness) <- list(dimnames(m.fitness)[[1]], as.character(pmutation)) 
# #dimnames(m.fitness) <- list(as.character(pcrossover), as.character(pmutation))
# #save(pmutation, pcrossover, fitness, file = "fitness.Rdata")  
# print(Sys.time()-start)
# flush.console()     # get finish time...
# return (m.fitness)
#}
#
#p.crossover <- function()
#{
# start=Sys.time() # get initial time...
# LEN <- 4
# LENX <- 20
# LENY <- 20
# pcrossover <- seq(length = LENX, from = 0.0, by =0.05)
# #pcrossover <- c(0, seq(length = LENY - 1, from = 0.1, by = 0.045))
# #fitness <- matrix(nrow = LENX, ncol = LENY)
# m.fitness <- matrix(nrow = LEN, ncol = LENX)
# load("cellcycle384x17.RData")
# m1 <- as.matrix(daisy(cellcycle384x17,  metric = "euclidean", stand = TRUE))
# 
# PMUT <- 0.00 
# for (i in 1:LEN){ k <- 0; print(i)
#  for (j in pcrossover)
#   {
#    k <- k + 1 
#    PCROSS <- j
#    m.fitness[i, k] <- agnes.gas(m1, 3, 100, part = 12/13)
#   }   
# } 
# dimnames(m.fitness) <- list(dimnames(m.fitness)[[1]], as.character(pcrossover)) 
# #dimnames(m.fitness) <- list(as.character(pcrossover), as.character(pmutation))
# #save(pmutation, pcrossover, fitness, file = "fitness.Rdata")  
# print(Sys.time()-start)
# flush.console()     # get finish time...
# return (m.fitness)
#}       
#

#-------------------------------------------------------------------------------
#                           ***********************


###############################################################
# From here begins the contents of file "evaluation.R".
# Functions of transformation and cluster validation.
################################################################ 


####################################################
# "H1" function.
# Gives the homogenity of a cluster, it not perform
# centroid.
####################################################

H1 <- function(cluster1)
{
 if (is.na(cluster1[2]))
  return (max(myenv$DISTANCE)) # if cluster has
  else                   # size 1.
   {
    c <- length(cluster1)
    return (sum(myenv$DISTANCE[cluster1, cluster1]) / (c * (c - 1)))
   }
}

####################################################
# "S1" function.
# Gives the separation betwen two clusters, it not
# perform centroid.
####################################################

S1 <- function(cluster1, cluster2)
 mean(myenv$DISTANCE[cluster1, cluster2]) # disjoin clusters.
 
####################################################
# "lencluster" function.
# Gives the pairwise count into a cluster.
####################################################

lencluster <- function(cluster1)
{
 c <- length(cluster1)
 if (c == 1) return (1)
 return ((c * (c - 1)) %/% 2)
}

####################################################
# "H.ave" function.
# Gives the Homogenity of a clustering, how Shamir
# and Sharan definition but without centroid.
####################################################

H.ave <- function(clustering)
{
 l <- length(clustering); clust <- clustering[-l]
 h1 <- sapply(clust, "H1")
 len <- sapply(clust, "lencluster")
 h1 <- h1 * len
 return (sum(h1) / sum(len))
}

####################################################
# "S.ave" function.
# Gives the Separation of a clustering, how Shamir
# and Sharan definition but without centroid.
####################################################

S.ave <- function(clustering)
{
 mean1 <- 0; sum1 <- 0
 c <- length(clustering) - 1
 # mean of distances between all clusters pair.
 for (i in 1:(c - 1))
  {
   len0 <- length(clustering[[i]])   
   for (j in (i + 1):c)
    {
     len <- len0 * length(clustering[[j]])
     mean1 <- mean1 + len * S1(clustering[[i]], clustering[[j]])
     sum1 <- sum1 + len
    }
  }  
 return (mean1 / sum1) 
}

    
####################################################
# "a" function.
# Gives the gene average distance "gene" to other  
# genes in the same cluster.
####################################################

a <- function(gene, clustering)
{
 i <- 1; c <- length(clustering) - 1
 while (i <= c)
  if (gene %in% clustering[[i]]) break()
   else i <- i + 1
 cluster1 <- setdiff(clustering[[i]], gene)
 if (is.na(cluster1[1])) return (max(myenv$DISTANCE))
 return (S1(c(gene), cluster1))  
}

####################################################
# "b" function.
# Gives the gene average distance "gene" to genes 
# in its nearest neighbor cluster.
####################################################

b <- function(gene, clustering)
{
 least <- pos <- Inf
 for (i in 1:(length(clustering) - 1))
  if (!(gene %in% clustering[[i]]))
   {
    d <- S1(c(gene), clustering[[i]])
    if (d < least) least <- d
   }
 return (least)  
}

####################################################
# "silhouette.g" function.
# Gives the Silhouette of a gene in a clustering.
####################################################

silhouette.g <- function(gene, clustering)
{
 a.g <- a(gene, clustering)
 b.g <- b(gene, clustering)
 # Computing the max{a.g, b.g}
 major <- (a.g + b.g + abs(a.g - b.g)) / 2 
 # Return the Silhouette of "gene"
 return ((b.g - a.g) / major)
}

####################################################
# "silhouette.mean" function.
# Gives the Silhouette mean of all genes in a 
# clustering.
####################################################

silhouette.mean <- function(...)
{
 add <- 0
 for (i in 1:myenv$OBJ.COUNT)
  add <- add + silhouette.g(i, ...)
 return (add / myenv$OBJ.COUNT) 
} 

####################################################
# "transf.clusters" function.
# Transform a clusters vector in a clusters list.
####################################################

transftoclusters <- function(vector.clust, orderdata, k){
  ordercluster <- function(clustering, orderdata, k){
     setclust <- function(i){
        elem <- orderdata[pos]; j <- 1
        while (j <= k)
          if (elem %in% clustering[[j]]) break()
            else j <- j + 1
        len <- length(clustering[[j]])
        tempclust[[i]] <<- orderdata[pos:(pos + len - 1)]
        pos <<- pos + len
     }
 
     tempclust <- list(NA); pos <- 1
     sapply(1:k, setclust, USE.NAMES = F)
     
     return (tempclust)
  }

  createcluster <- function(i)
    clustering[[i]] <<- which(vector.clust == i)
     
  clustering <- list(NA)
  vector.clust <- as.vector(vector.clust)
  sapply(1:k, createcluster, USE.NAMES = F)
  clustering <- ordercluster(clustering, orderdata, k)
  clustering[[k + 1]] <- c(0, 0)
  
  return (clustering)
}

transftovector <- function(clustering)
{
 k <- length(clustering) - 1
 tvector <- rep(0, myenv$OBJ.COUNT)
 for (i in 1:k)
  tvector[clustering[[i]]] <- i
 return (tvector) 
}

transf.tree <- function(agn, part = 2){
  setlevels <- function(i){
    tree[[n - i]] <<- transftoclusters(cutree(agn, k = size), orderv, size)
    size <<- size - 1
  }  
 
  n <- length(agn$order)
  orderv <- agn$order
  tree <- list(NULL)
  if (part != 2)  part <- trunc(n * part) + 2
  size <- n - part + 1
  sapply((n - 1):part, setlevels, USE.NAMES = F)

  return (tree)
}

####################################################
# "degree.simalarity" function.
# Agreement with a reference partition.
####################################################

cluster.similarity <- function(cluster, cmat)
{
 k <- length(cluster)
 for (i in 1:k-1)
  for (j in (i+1):k){
    cmat[cluster[i], cluster[j]] <- 1
    cmat[cluster[j], cluster[i]] <- 1
  } 
 return (cmat)  
}
 
degree.similarity <- function(clustering, groundtruth) 
{
 cmat <- (pmat <- matrix(0, nrow = myenv$OBJ.COUNT, ncol = myenv$OBJ.COUNT, byrow = TRUE))
 k <- length(clustering) - 1
 for (i in 1:k)
  cmat <- cluster.similarity(clustering[[i]], cmat)   # V
 k <- length(groundtruth) - 1 
 for (i in 1:k)
  pmat <- cluster.similarity(groundtruth[[i]], pmat)   # U
 nij <- integer(4) 
 for (i in 1:myenv$OBJ.COUNT)
  for (j in 1:myenv$OBJ.COUNT)
   {
    if (cmat[i, j] == 1 && pmat[i, j] == 1) nij[1] <- nij[1] + 1   # a (agreement)
    if (cmat[i, j] == 1 && pmat[i, j] == 0) nij[2] <- nij[2] + 1   # c (disagreements)
    if (cmat[i, j] == 0 && pmat[i, j] == 1) nij[3] <- nij[3] + 1   # b (disagreements)
    if (cmat[i, j] == 0 && pmat[i, j] == 0) nij[4] <- nij[4] + 1   # d (agreement)
   } 
 return (nij) 
}

# The Rand index lies between 0 and 1.
# When the two partitions agree perfectly, the Rand index is 1.
rand.index <- function(nij)
 (nij[1] + nij[4]) / sum(nij)

# W. M. Rand (1971). "Objective criteria for the evaluation of 
# clustering methods". Journal of the American Statistical 
# Association 66: 846–850. doi:10.2307/2284239.  
#a.rand.index <- function(nij)
# 2 * (nij[1] * nij[4] - nij[2] * nij[3]) / 
# ((nij[1] + nij[3]) * (nij[3] + nij[4]) + 
#  (nij[1] + nij[2]) * (nij[2] + nij[4]))

a.rand.index <- function(nij){
 a <- nij[1]; b <- nij[3]; c1 <- nij[2]; d <- nij[4]
 temp <- (a + b) * (a + c1)
 
 return ( (a * choose(myenv$OBJ.COUNT, 2) - temp) / 
          (1/2 * choose(myenv$OBJ.COUNT, 2) * (2 * a + b + c1) - temp) )
}  

a.rand.index1 <- function(nij){
 a <- nij[1]; b <- nij[3]; c1 <- nij[2]; d <- nij[4]
 return ( ( choose(myenv$OBJ.COUNT, 2) * (a + d) - ((a + b) * (a + c1) + (c1 + d) * (b + d)) ) /
        ( choose(myenv$OBJ.COUNT, 2)^2 - ((a + b) * (a + c1) + (c1 + d) * (b + d)) ) )
}

# Jacard and Minkowski coefficient may be more effective in 
# gene-based clustering than Rand index (Jiang2004).
jacard.coef <- function(nij)  
 nij[1] / (sum(nij) - nij[4])                   
 
minkowski <- function(nij)
 sqrt((nij[2] + nij[3]) / (nij[1] + nij[3]))


#####################################################################
# Ajusted Rand Index with contingency matrix
#####################################################################

# Intersecta dos clusters y devuelve la longitud del resultado
# de la intersección.
intersect.clusters <- function(cluster1, cluster2 = NULL)
 length(intersect(cluster1, cluster2))

# Intersecta una partición con todos los clusters de un clustering.
# Toma la logitud de cada intersección y devuelve la suma total de
# esas longitudes.
intersectclusters.sum <- function(partition, clustering = NULL)
 sum(sapply(clustering[1:(length(clustering) - 1)], intersect.clusters,
            cluster2 = partition, USE.NAMES = F))
            
# Intersecta una partición con todos los clusters de un clustering.
# Toma la logitud de cada intersección y devuelve un vector con
# las longitudes.
intersectclusters.vector <- function(partition, clustering = NULL)
  sapply(clustering[1:(length(clustering) - 1)], intersect.clusters,
         cluster2 = partition, USE.NAMES = F)

# Devuelve un vector o una matriz, según el valor de functname.
# Si se utiliza "functname = intersectclusters.sum", se devuelve
# un vector columna de totales por filas, de las longitudes de la
# intersecciones de cada partición de "refpartition" con cada
# cluster de "clustering".
#
# En caso contrario, es decir con "functname = intersectclusters.vector", se
# devuelve una matriz (matriz de contingencia) con las longitudes de las
# intersecciones de cada partición de "refpartition" con cada cluster de
# "clustering".
intersectclusters.total <- function(refpartition, clustering,
                                    functname = "intersectclusters.sum")
   sapply(refpartition[1:(length(refpartition) - 1)], functname,
          clustering = clustering, USE.NAMES = F)
          
#####################################################################
# Ajusted Rand Index
#####################################################################

a.rand.index2 <- function(refpartition, clustering){
  n <- myenv$OBJ.COUNT 
  # Totales de la matriz de contingencia por filas.
  rowtotal <- intersectclusters.total(refpartition, clustering)
  # Totales de la matriz de contingencia por columnas.
  coltotal <- intersectclusters.total(clustering, refpartition)
  # Total de los binomios de los valores de la matriz de contingencia,
  # tomados dos a dos "binomio(matriz(i, j), 2)".
  a <- sum(choose(intersectclusters.total(clustering, refpartition,
                                    functname = "intersectclusters.vector"), 2))
  # Suma total de los totales por filas de la matriz de contingencia.
  tr <- sum(choose(rowtotal, 2))
  # Suma total de los totales por columnas de la matriz de contingencia.
  tc <- sum(choose(coltotal, 2))
  # Valor esperado de "Ajusted Rand Index"
  expectedvalue <- tr * tc
  
  numerator <- choose(n, 2) *  a - expectedvalue

  denominator <- 1/2 * choose(n, 2) * (tr + tc) - expectedvalue
  
  return (numerator / denominator)
}
