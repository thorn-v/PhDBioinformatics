#requires poppr at the moment!

# workhorse function to check 2 loci vectors for how many combos there are - only one allele per element
FGT <- function(l1, l2){
  # table does not include NAs by default - maybe have to add checks for other chars (ie. !(ACGT)) to be converted to NA
  cross <- table(l1, l2)
  # table makes a punnet square basically - so each element over 0 is a combo seen in the set
  return(sum(cross > 0))
}


### Start Here ####

Pinfdf <- #Your data - needs to be loci in columns and samples to be rows; missing data needs to be encoded as NA

# gets all possible combos of input of size given (2 in this case for every combo of 2 loci); will not do with self so no diagonal and counts l1/l2 as the same as l2/l1
loci <- combn(colnames(Pinfdf), 2)

results.fgt <- matrix(data = 0, nrow = length(colnames(Pinfdf)), ncol = length(colnames(Pinfdf)))
colnames(results.fgt) <- colnames(Pinfdf)
rownames(results.fgt) <- colnames(Pinfdf)

## actually do the test ##
# results in an upper triangle
for(i in 1:ncol(loci)){
  x <- loci[1, i]
  y <- loci[2, i]
  results.fgt[x, y] <- FGT(Pinfdf[,x], Pinfdf[,y])
}


### Significance testing ####
# currently takes a bit, depending on the number of loci and sample size
# TODO make this into a parallel loop

fgt.boot <- vector(length=999) # 999 reps gives you 0.001 p val cut off - will hold the number of 4 gamete loci pairs arose in each cycle

for(i in 1:999){
  set.seed(i) # sets seed each time so we can track the randomness, but still get a random shuffle each loop iteration - if needed, can modify i for new seeds (ex. i*3 or i-(2*i) )
  
  temp.fgt <- shufflepop(   noRoot.cc   , method = 1) #put your genind object in the big space
  Pinfdf.boot <- genind2df(temp.fgt, sep = "")
  
  loci.boot <- combn(colnames(Pinfdf.boot)[], 2) #if your resultant data frame has other columns that are not loci (like population info) then make sure to not include those columns (can use the [] to subset)
  
  fgt.vec <- vector(length = ncol(loci)) # don't care about which pairs have what number of gametes so we only need a vector
  
  ## actually do the test ##
  for(j in 1:ncol(loci)){
    x <- loci[1, j]
    y <- loci[2, j]
    fgt.vec[j] <- FGT(Pinfdf.boot[,x], Pinfdf.boot[,y])
  }
  
  fgt.boot[i] <- tail(table(fgt.vec), n=1) # gets the number of pairs that had all 4 gametes in this boot cycle

}


s0 <- as.vector(tail(table(results.fgt), n=1)) # gets the number of loci pairs that had 4 from your original data

hist(fgt.boot) # may need to change the x-lims to see your value (in red) against the distribution
abline(v = s0, col = "red") # your value

p.val <- (1+sum(fgt.boot <= s0))/(length(fgt.boot)+1) # formula for p value from random sampling 

# if you did 999 permutations, then you should have a 99.9% certainty I belive
p.val




