#install.packages("remotes")
#remotes::install_github("ysamwang/highDNG",subdir = "highDlingam")
#install.packages("RPtests")

library(doParallel)
library(highDLingam)
library(RPtests)

# sample data
set.seed(123)
num_samples <- 1000
x4 <- 2*runif(n = num_samples)
x1 <- 2*runif(n = num_samples) + 3 * x4
x3 <- 2*runif(n = num_samples) + 6 * x4
# x2 <- 2*runif(n = num_samples) + 3 * x1 + 2 * x3
# x6 <- 2*runif(n = num_samples) + 4 * x1
# x5 <- 2*runif(n = num_samples) + 8 * x1 - 1 * x3
# X  <- data.frame(x1, x2, x3, x4, x5, x6)
X  <- data.frame(x1, x3, x4)
write.csv(X,"99_sandbox/sampledata_dummy.csv",row.names=FALSE)

# Load  data
dataX <- as.matrix(read.table("99_sandbox/sampledata_dummy.csv", sep = ",", header = T))
dim(dataX)

### cores to use for parallel computation
ncores <- 3
cl <- makeCluster(ncores)
registerDoParallel(cl, setup_timeout = 0.5)

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


# Create the adjacency matrix B

adj_mat <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #adjacency matrix #�L���ӂ̗L���i0 or 1�j�̌��ʊi�[�p
colnames(adj_mat) <- colnames(dataX)
rownames(adj_mat) <- colnames(dataX)

for(i in 1:length(output3$topOrder)){
     k <- output3$topOrder[i]
     for (j in output3$parents[i]){
          adj_mat[k,j] <- 1
     }
} #

adj_mat;output3$topOrder;output3$parents #adj_mat : �L���ӂ̗L���i0 or 1�j�̌���(�s��)


# Draw the causal graph
# https://github.com/gkikuchi/rlingam/blob/master/R/plot_adjacency_mat.R
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
print(DiagrammeR::grViz(dot$x$diagram))


# Do bootstrapping and collect the bootstrap prob. in Bboot
B <- 1000 #�u�[�g�X�g���b�v������
n <- nrow(dataX) #�W�{�̑傫��

adj_mat <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #�L���ӂ̗L���i0 or 1�j�̍ŏI���ʂ̊i�[�p
colnames(adj_mat) <- colnames(dataX)
rownames(adj_mat) <- colnames(dataX)

adj_mat_boot <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #�L���ӂ̗L���i0 or 1�j�̓r�����ʂ̊i�[�p
colnames(adj_mat_boot) <- colnames(dataX)
rownames(adj_mat_boot) <- colnames(dataX)

#install.packages("tcltk")
library(tcltk)
#min�Ń��[�v�ϐ��̍ŏ��l�Amax�Ń��[�v�ϐ��̍ő�l�Astyle�ŕ\���X�^�C��(��{��3?)��ݒ�
pb <- txtProgressBar(min=1, max=B, style=3)

for (b in 1:B) { ### B�ɏC��
     bt <- sample(1:n,replace = TRUE) #�u�[�g�X�g���b�v�W�{�ԍ� #�������o
     mydata <- dataX[bt,] #�u�[�g�X�g���b�v�W�{
     output3 <- highDLingam::findGraphMulti(mydata, maxInDegree = md, cutOffScaling = cs, degree = deg, verbose = F)
     adj_mat_boot <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #�L���ӂ̗L���i0 or 1�j�̓r�����ʂ̊i�[�p��������
     for(i in 1:length(output3$topOrder)){
          k <- output3$topOrder[i]
          for (j in output3$parents[i]){
               adj_mat_boot[k,j] <- 1
          }
     } #create the adjacency matrix
     adj_mat <- adj_mat + adj_mat_boot #�s��̊e�v�f�𑫂��Z����
     setTxtProgressBar(pb, b) 
}
adj_mat #1000��J��Ԃ������A�u�[�g�X�g���b�v�L���ӂ̗L���i0 or 1�j�̍ŏI����

adj_mat_boot_prob <- adj_mat / B #�u�[�g�X�g���b�v�m��


# Save the results
adj_mat_boot_prob <- data.frame(adj_mat_boot_prob)
write.csv(adj_mat_boot_prob,"adj_mat_boot_prob.csv",row.names=FALSE)







########## Directlingam.R �Ƃ̐��茋�ʂ̔�r ##########

#remotes::install_github("gkikuchi/rlingam")
library(rlingam)

# Load  data
dataX <- as.matrix(read.table("sampledata_200922.csv", sep = ",", header = T))
dim(dataX)

# model fit
mdl <- rlingam::DirectLiNGAM$new()
mdl$fit(dataX)

print(mdl$causal_order) # order
print(mdl$adjacency_matrix) #adj mat
plot_adjacency_mat(mdl$adjacency_matrix, node_labels = names(X)) #plot

# Do bootstrapping and collect the bootstrap prob. in Bboot
B <- 1000 #�u�[�g�X�g���b�v������
n <- nrow(dataX) #�W�{�̑傫��

adj_mat <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #�L���ӂ̗L���i0 or 1�j�̍ŏI���ʂ̊i�[�p
colnames(adj_mat) <- colnames(dataX)
rownames(adj_mat) <- colnames(dataX)

adj_mat_boot <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #�L���ӂ̗L���i0 or 1�j�̓r�����ʂ̊i�[�p
colnames(adj_mat_boot) <- colnames(dataX)
rownames(adj_mat_boot) <- colnames(dataX)

#install.packages("tcltk")
library(tcltk)
#min�Ń��[�v�ϐ��̍ŏ��l�Amax�Ń��[�v�ϐ��̍ő�l�Astyle�ŕ\���X�^�C��(��{��3?)��ݒ�
pb <- txtProgressBar(min=1, max=B, style=3)

for (b in 1:B) { ### B�ɏC��
        bt <- sample(1:n,replace = TRUE) #�u�[�g�X�g���b�v�W�{�ԍ� #�������o
        mydata <- dataX[bt,] #�u�[�g�X�g���b�v�W�{
        #mdl <- rlingam::DirectLiNGAM$new()
        mdl$fit(mydata)
        
        adj_mat_boot <- mdl$adjacency_matrix #adj mat
        adj_mat <- adj_mat + adj_mat_boot
        
        #output3 <- highDLingam::findGraphMulti(mydata, maxInDegree = md, cutOffScaling = cs, degree = deg, verbose = F)
        #for(i in 1:length(output3$topOrder)){
        #        k <- output3$topOrder[i]
        #        for (j in output3$parents[i]){
        #                adj_mat_boot[k,j] <- 1
        #        }
        #} #create the adjacency matrix
        #adj_mat <- adj_mat + adj_mat_boot #�s��̊e�v�f�𑫂��Z����
        setTxtProgressBar(pb, b) 
}
adj_mat #1000��J��Ԃ������A�u�[�g�X�g���b�v�L���ӂ̗L���i0 or 1�j�̍ŏI����

adj_mat / B #�u�[�g�X�g���b�v�m��





















##########
# Directlingam.R
#remotes::install_github("gkikuchi/rlingam")
library(rlingam)

X <- rlingam::gen_dummy_data(random_state = 10) # data
mdl <- rlingam::DirectLiNGAM$new()
mdl$fit(X)

print(mdl$causal_order) # order
print(mdl$adjacency_matrix) #adj mat
plot_adjacency_mat(mdl$adjacency_matrix, node_labels = names(X)) #plot
##########

mdl$adjacency_matrix # adj mat

node_labels <- adj #label
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
print(DiagrammeR::grViz(dot$x$diagram))


##### for���������̐i���m�F
library(tcltk)
n <- 100
#min�Ń��[�v�ϐ��̍ŏ��l�Amax�Ń��[�v�ϐ��̍ő�l�Astyle�ŕ\���X�^�C��(��{��3?)��ݒ�
pb <- txtProgressBar(min=1, max=n, style=3)
for(i in 1:n){
        #���������v���O���X�o�[(pb)�ƃ��[�v�ϐ�(i)���w��
        setTxtProgressBar(pb, i) 
}
#####