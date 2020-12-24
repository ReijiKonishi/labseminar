#install.packages("remotes")
#remotes::install_github("ysamwang/highDNG",subdir = "highDlingam")
#install.packages("RPtests")

rm(list = ls())

library(doParallel)
library(highDLingam)
library(RPtests)

# sample data
set.seed(123)
num_samples <- 1000
x4 <- 2*runif(n = num_samples)
x1 <- 2*runif(n = num_samples) + 3 * x4
x3 <- 2*runif(n = num_samples) + 6 * x4 + 3 * x1
# x2 <- 2*runif(n = num_samples) + 3 * x1 + 2 * x3
# x6 <- 2*runif(n = num_samples) + 4 * x1
# x5 <- 2*runif(n = num_samples) + 8 * x1 - 1 * x3
# X  <- data.frame(x1, x2, x3, x4, x5, x6)
X  <- data.frame(x1, x3, x4)
write.csv(X,"sampledata_dummy.csv",row.names=FALSE)

### Load  data
dataX <- as.matrix(read.table("sampledata_dummy.csv", sep = ",", header = T))
dim(dataX)

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
output3 <- highDLingam::findGraphMulti(dataX, maxInDegree = md, cutOffScaling = cs, degree = deg, verbose = F)

output3

### Create the adjacency matrix 
adj_mat <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #adjacency matrix #�L���ӂ̗L���i1 or 0�j�̌��ʂ��i�[
colnames(adj_mat) <- colnames(dataX);rownames(adj_mat) <- colnames(dataX)

for(i in 1:length(output3$topOrder)){
     k <- output3$topOrder[i]
     for (j in output3$parents[i]){
          adj_mat[k,j] <- 1
     }
}

adj_mat #adjacency matrix : �L���ӂ̗L���i1 or 0�j�̌��ʂ�\��
 

### Draw the causal graph
# �Q�l�R�[�h�Fhttps://github.com/gkikuchi/rlingam/blob/master/R/plot_adjacency_mat.R
node_labels <- colnames(adj_mat)
round_digit <- 2
nodes <- DiagrammeR::create_node_df(ncol(adj_mat), label = node_labels)
graph <- DiagrammeR::create_graph(nodes, directed = TRUE)
adj_mat <- round(adj_mat, round_digit)
for (i in 1:ncol(adj_mat)) {
     for (j in which(abs(adj_mat[, i]) > 0)) {
          graph <- DiagrammeR::add_edge(graph, from = i, to = j, edge_aes = DiagrammeR::edge_aes(label = " "))
     }
}
dot <- DiagrammeR::render_graph(graph, layout = "dot")
dot$x$diagram <- gsub("neato", "dot", dot$x$diagram)
print(DiagrammeR::grViz(dot$x$diagram)) #�O���t��`��


### Do bootstrapping and collect the bootstrap prob. in B
B <- 1000 #�u�[�g�X�g���b�v������
n <- nrow(dataX) #�W�{�̑傫��

adj_mat_boot <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #�L���ӂ̗L���i1 or 0�j�̍ŏI���ʂ̊i�[�p
colnames(adj_mat_boot) <- colnames(dataX);rownames(adj_mat_boot) <- colnames(dataX)

adj_mat <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #�L���ӂ̗L���i1 or 0�j�̓r�����ʂ̊i�[�p
colnames(adj_mat) <- colnames(dataX);rownames(adj_mat) <- colnames(dataX)

# for���������̐i���m�F
#install.packages("tcltk")
library(tcltk)
#min�Ń��[�v�ϐ��̍ŏ��l�Amax�Ń��[�v�ϐ��̍ő�l�Astyle�ŕ\���X�^�C��(��{��3)��ݒ�
pb <- txtProgressBar(min=1, max=B, style=3)

for (b in 1:B) { 
     bt <- sample(1:n,replace = TRUE) #�u�[�g�X�g���b�v�W�{�ԍ� #�������o
     mydata <- dataX[bt,] #�u�[�g�X�g���b�v�W�{
     output3 <- highDLingam::findGraphMulti(mydata, maxInDegree = md, cutOffScaling = cs, degree = deg, verbose = F)
     adj_mat <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #�L���ӂ̗L���i1 or 0�j�̓r�����ʊi�[�p��������
     for(i in 1:length(output3$topOrder)){
          k <- output3$topOrder[i]
          for (j in output3$parents[i]){
               adj_mat[k,j] <- 1
          }
     } #create the adjacency matrix
     adj_mat_boot <- adj_mat_boot + adj_mat #�ŏI���ʂ̊i�[�p�ɁA�L���ӂ̗L���i1 or 0�j�������X�V����
     setTxtProgressBar(pb, b) 
}
adj_mat_boot #B��J��Ԃ������́A�L���ӂ̗L���i1 or 0�j�̍ŏI���ʂ�\��

adj_mat_boot_prob <- adj_mat_boot / B 
adj_mat_boot_prob #B��J��Ԃ������̃u�[�g�X�g���b�v�m����\��

### Save the results
adj_mat_boot <- data.frame(adj_mat_boot)
write.csv(adj_mat_boot,"adj_mat_boot.csv",row.names=TRUE)

adj_mat_boot_prob <- data.frame(adj_mat_boot_prob) #�u�[�g�X�g���b�v�m��
write.csv(adj_mat_boot_prob,"adj_mat_boot_prob.csv",row.names=TRUE)