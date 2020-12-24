install.packages("remotes")
remotes::install_github("ysamwang/highDNG",subdir = "highDlingam")
install.packages("RPtests")

library(doParallel)
library(highDLingam)
library(RPtests)

# Load  data
dir <- getwd()
dir <- paste0(dir, "/99_sandbox")
setwd(dir)

dataX <- as.matrix(read.table("sampledata_gasshuku.csv", sep = ",", header = T))
dim(dataX)

x2 <- runif(1000)
x1 <- 2.0 * x2 + runif(1000)
x3 <- 3.0 * x2 + 4.0 * x2 + runif(1000)

dataX <- as.matrix(data.frame(x1 , x2, x3)) 

### cores to use for parallel computation
ncores <- 3
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

### assumed max in-degree (J in paper)
md <- 3
### Cut-off scaling which corresponds (alpha parameter in the paper)
cs <- 1
### Degree of moment to check for test statistic (K in paper)
deg <- 4

### Estimate topological ordering (used in Section 4.4)
?highDLingam::findGraphMulti
output3 <- highDLingam::findGraphMulti(dataX, maxInDegree = md, cutOffScaling = cs, degree = deg, verbose = T)

output3

# Create the adjacency matrix B

# Draw the causal graph

# Do bootstrapping and collect the  bootstrap prob. in Bboot
nboot <- 200

# Save the results

