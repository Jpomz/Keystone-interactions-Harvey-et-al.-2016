
#... Script: code to rank interactions in an ecological network by their impact on structural stability
#......................................................................................................
#... Supplementary information of for Harvey et al. (2016). Bridging ecology and conservation: a road map to ecological network protection. Journal of applied ecology.
#... Script author: Isabelle Gounand, e-mail: isabelle.gounand@eawag.ch 
#... Date: March 18, 2016


#... References: 
# Sauve et al. 2016 Ecology 
# Neutel et al. 2002 Science, 2007 Nature
# Allesina and Tang 2012 Nature

#.......................................................................................................
#............ Example of a food web, using the community stability index as in Sauve et al. 2016 Ecology
#.......................................................................................................
# The community stability index represent the minimum intraspecific competition for the system to be stable (the smaller the more stable)
# This community stability index can be used if there is no cannibalism (only zero on the diagonal of the adjacency matrix)
# If there are some of these interactions, use directly the resilience as index of stability 
# (= real part of the dominant eigenvalue of the jacobian matrix): the more negative, the more resilient

#... Analysis steps: 
# Step 1: Calculate the change in structural stability with the removal of each interaction separately
# Step 2: Calculate the number of extinctions following the removal of each interaction separately



#.......................................................................................................
#... Clear workspace
#...................
rm(list=ls())


#... Load the stability functions (Online Supplementary material of Sauve et al. 2016 Ecology)
#................................
source("C:\\Users\\Justin\\Documents\\Literature\\Sauve et al 2016\\stability_functions.R")


#... Load data
#.............
#... Load the adjacency matrix
# the adjacency matrix m is a square matrix of size S, with S the total number of species
# undirected interactions (or bidirectional) are such that m_ij = m_ji = 1
# directected interactions (such as antagonistic ones) are such that m_ij = 1 and m_ji = 0 if i feeds on j (j -> i)
# i are the rows, with species in their role of consumer in antagonistic interactions
# j are the columns, with species in their role of resource in antagonistic interactions
m = as.matrix(read.csv("adjacency_matrix_foodweb.csv",h=T, row.names=1, sep=";"))
s = nrow(m) # the number of species

#... Load the functional role of species basal (primary producers) or not
# this is necessary to calculate robustess and secondary extinctions following removal of interactions
# b is a array representing species in the same order than in the adjacency matrix
# 1 is for basal species (persist without any resource), 0 for consumer (need an interaction with a resource to persist)   
b = t(as.matrix(read.csv("basal_species_foodweb.csv",h=F, row.names=1, sep=";")))


#... Variable to initiate the seeds (arbitrary positive integer)
#..................................
move_seed = 13000


#... Define the Interaction links
#................................

links = which(m!=0) # continuous position in the matrix

#... function to get the row - column coordinates from continuous positions in a matrix
# x is the position
# s is the size of the square matrix (number of species)
get_coordinates <-function(x){
  row = ifelse(x%%s!=0,x%%s,s) #T == x%%s; F == s
  col = ifelse(x%%s!=0,trunc(x/s)+1,trunc(x/s))
  return(c(row,col))
}

#... names of the links
names.links = as.character(links)
for(i in 1:length(links)) names.links [i] = paste(rownames(m)[get_coordinates(links[i])[1]],"->",colnames(m)[get_coordinates(links[i])[2]],sep="")

#... Get the row - column coordinates of the links
links_ij = t(sapply(links,FUN = get_coordinates))
colnames(links_ij) = c("row","col")
rownames(links_ij) = names.links



#....................................................................................................... 
#.......................................................................................................
#... Step 1:  Change in structural stability with the removal of each interaction
#................................................................................

#... Calculate the initial stability index
#.........................................
# the seed is set in order to keep the same set of interaction strength parameterization for each interaction removal
# when only one interaction strength is changed
nbiter = 10000  # number of iteration for random parameterization of interaction strengths (the higher the more accurate the final ranking of interactions)
stab = numeric(nbiter) # array to record community stability index for each random parameterization

#... Loop for random parameterizations of interaction strength
for(i in 1:nbiter)
{
  #... Set a new (known) seed
  set.seed(i+move_seed) 
  
  #... Parameterize the Jacobian matrix (function from Sauve et al. 2016's script)
  # interaction strengths drawn from a normal distribution with mean 1 and sd 0.1
  # m is the adjacency matrix with 0 and 1
  J = jacobian_binary(m=m) 
  
  #... Calculate the community stability index (function from Sauve et al. 2016's script)
  # s2 is an arbitrary number; if it is too low, then increase it
  stab[i] = stability(J,s2 = 1) # 
}

#... Central metric of initial stability
# median of the distribution of stability, calculate from many random parameterizations
stab_init = median(stab)


#... Sequential interaction removal and effect on structural stability
#.....................................................................
delta_stab <-numeric(length(links)) # array to record change in stability when an interaction is removed

#... Loop for link removals
#..........................
for (j in 1:length(links)){
  m2 = m
  m2 [links_ij[j,1],links_ij[j,2]] = 0
  stab_changed = numeric(nbiter) # array to record community stability indices when 1 interaction is removed
  for(i in 1:nbiter)
  {
    set.seed(i+move_seed) # Set the same seeds as for initial stability calculation
    J2 <- jacobian_binary(m=m2) # Parameterize the interaction strengths
    stab_changed[i] <- stability(J2,1) # Calculate the community stability index
  }
  
  #... change in average stability
  # Percent of change compared to initial stability
  # Negative numbers mean that stability decreases
  delta_stab[j] <- ((stab_init - median(stab_changed))/stab_init*100)
}


#... Rank the interactions depending on stability change
#.......................................................
names(delta_stab) <- names.links
delta_stab_sorted <- sort(delta_stab,decreasing = T)


#... Plot the stability change against ranked interactions
#.........................................................
colos = rainbow(length(links),start=0,end=0.35) # colors of the points: red (stability loss) to green (stability gain)
pdf("Structural_stability_changes_seed13000_nbiter10000.pdf")
plot(1:length(links),delta_stab_sorted,xlab = "Links",ylab="% Change in structural stability",xaxt="n",pch=19,col=rev(colos))
axis(1,at=1:length(links),labels=names(delta_stab_sorted),las=2,cex.axis=0.8)
abline(h=0,lty=3)
dev.off()



#....................................................................................................... 
#.......................................................................................................
#... Step 2:  number of extinctions following the removal of each interaction
#............................................................................

#... Function to calculate the number of extinct species following the removal of interactions
# am is the adjacency matrix, to which some links have been removed (1's turned into 0's)
# basal is the array of species function where basal species are indicated by a 1 (otherwise 0)
# the species have to be in the same order than in the square adjacency matrix (rows or columns)
# basal species are those which doesn't need an interaction with a resource species to persist
nb.extinct = function(am,basal){
  am2 = am
  basal2 = basal
  is.extinct = function(x){res = ifelse(sum(x[-1])==0 & x[1]==0,1,0);return(res)} # function determining if a species in row is extinct
  while(sum(apply(cbind(c(basal2),am2),1,is.extinct))!=0) # Loop until no more species is extinct
  {
    extinct = which(apply(cbind(c(basal2),am2),1,is.extinct)==1) # Positions of species extinct
    am2 = am2[-extinct,-extinct] # removal of extinct species in the adjacency matrix
    basal2 = basal2[-extinct] # removal of extinct species in the array of species function
  }
  return(length(basal) - nrow(am2)) # return the number of extinct species
}

#... Define the array to record the number of extinctions following one interaction removal
#..........................................................................................
nb_ext = numeric(length(links)) 
names(nb_ext) = names.links


#... Loop for link removals
#..........................
for (j in 1:length(links)){
  m2 = m
  m2 [links_ij[j,1],links_ij[j,2]] = 0
  nb_ext[j] = nb.extinct(m2,b)
}

#... Print the result
#....................
print(nb_ext)




