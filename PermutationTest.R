# ON POOLS
#load data
pools <- read.table(file ="/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/PermutationCalculations.csv", header = T, sep = ",")

# calculate the difference in sample MEANS 
mean(pools$count449[pools$validbinary == "1"]) # mean for validmuts 
mean(pools$count449[pools$validbinary == "0"]) # mean for validWT 

# lets calculate the absolute diff in means 
test.stat1 <- abs(mean(pools$count449[pools$validbinary == "1"]) -  
                    mean(pools$count449[pools$validbinary == "0"]))  
test.stat1 

# seed for reproducability 
set.seed(1979)   
# the number of pools to sample 
n <- length(pools$sample)   
# the number of permutations 
P <- 100000  
# the variable we will resample from  
variable <- pools$count449
# initialize a matrix to store the permutation data 
PermSamples <- matrix(0, nrow = n, ncol = P) 
# each column is a permutation sample of data 
# now, get those permutation samples, using a loop 
# let's take a moment to discuss what that code is doing 
for(i in 1:P) 
{ 
  PermSamples[, i] <- sample(variable,  
                             size = n,  
                             replace = FALSE) 
} 

# we can take a quick look at the first 5 columns of PermSamples 
PermSamples[, 1:5] 

#check same number of 1s and 0s as data 
sum(pools$count449)
sum(PermSamples[,23])
length(PermSamples[,23])

# initialize vectors to store all of the Test-stats 
Perm.test.stat1 <- rep(0, P) 

# loop thru, and calculate the test-stats 
for (i in 1:P) { 
  # calculate the perm-test-stat1 and save it 
  Perm.test.stat1[i] <- abs(mean(PermSamples[pools$validbinary == "1",i]) -  
                              mean(PermSamples[pools$validbinary == "0",i])) 
  
  # Perm.test.stat2[i] <-   abs(mean(PermSamples[pools$count449 == "1",i]) -  
  #       mean(PermSamples[pools$count449 == "0",i]))
} 
  
# the TEST STATS 
test.stat1
  
# permutation-TEST STATS for 1  
round(Perm.test.stat1[1:15], 1)

# and, let's calculate the permutation p-value 
# notice how we can ask R a true/false question 
(Perm.test.stat1 >= test.stat1)[1:180] 

# and if we ask for the mean of all of those, 
# it treats 0 = FALSE, 1 = TRUE 
mean((Perm.test.stat1 >= test.stat1)[1:15]) 

# Calculate the p-value, for all P = 100,000 
mean(Perm.test.stat1 >= test.stat1) 



#Test Strict Cand vs Valid
# calculate the difference in sample MEANS 
mean(pools$candbin[pools$validbinary == "1"]) # mean for validmuts 
mean(pools$candbin[pools$validbinary == "0"]) # mean for validWT 

# lets calculate the absolute diff in means 
test.stat1 <- abs(mean(pools$candbin[pools$validbinary == "1"]) -  
                    mean(pools$candbin[pools$validbinary == "0"]))  
test.stat1 

# seed for reproducability 
set.seed(1979)   
# the number of pools to sample 
n <- length(pools$sample)   
# the number of permutations 
P <- 100000  
# the variable we will resample from  
variable <- pools$candbin
# initialize a matrix to store the permutation data 
PermSamples <- matrix(0, nrow = n, ncol = P) 
# each column is a permutation sample of data 
# now, get those permutation samples, using a loop 
# let's take a moment to discuss what that code is doing 
for(i in 1:P) 
{ 
  PermSamples[, i] <- sample(variable,  
                             size = n,  
                             replace = FALSE) 
} 

# we can take a quick look at the first 5 columns of PermSamples 
PermSamples[, 1:5] 

#check same number of 1s and 0s as data 
sum(pools$candbin)
sum(PermSamples[,23])
length(PermSamples[,23])

# initialize vectors to store all of the Test-stats 
Perm.test.stat1 <- rep(0, P) 

# loop thru, and calculate the test-stats 
for (i in 1:P) { 
  # calculate the perm-test-stat1 and save it 
  Perm.test.stat1[i] <- abs(mean(PermSamples[pools$validbinary == "1",i]) -  
                              mean(PermSamples[pools$validbinary == "0",i])) 
  
  # Perm.test.stat2[i] <-   abs(mean(PermSamples[pools$count449 == "1",i]) -  
  #       mean(PermSamples[pools$count449 == "0",i]))
} 

# the TEST STATS 
test.stat1

# permutation-TEST STATS for 1  
round(Perm.test.stat1[1:15], 1)

# and, let's calculate the permutation p-value 
# notice how we can ask R a true/false question 
(Perm.test.stat1 >= test.stat1)[1:180] 

# and if we ask for the mean of all of those, 
# it treats 0 = FALSE, 1 = TRUE 
mean((Perm.test.stat1 >= test.stat1)[1:15]) 

mean(Perm.test.stat1 >= test.stat1) 

