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
adj_mat <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #adjacency matrix #有向辺の有無（1 or 0）の結果を格納
colnames(adj_mat) <- colnames(dataX);rownames(adj_mat) <- colnames(dataX)

for(i in 1:length(output3$topOrder)){
     k <- output3$topOrder[i]
     for (j in output3$parents[i]){
          adj_mat[k,j] <- 1
     }
}

adj_mat #adjacency matrix : 有向辺の有無（1 or 0）の結果を表示
 

### Draw the causal graph
# 参考コード：https://github.com/gkikuchi/rlingam/blob/master/R/plot_adjacency_mat.R
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
print(DiagrammeR::grViz(dot$x$diagram)) #グラフを描画


### Do bootstrapping and collect the bootstrap prob. in B
B <- 1000 #ブートストラップ反復回数
n <- nrow(dataX) #標本の大きさ

adj_mat_boot <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #有向辺の有無（1 or 0）の最終結果の格納用
colnames(adj_mat_boot) <- colnames(dataX);rownames(adj_mat_boot) <- colnames(dataX)

adj_mat <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #有向辺の有無（1 or 0）の途中結果の格納用
colnames(adj_mat) <- colnames(dataX);rownames(adj_mat) <- colnames(dataX)

# for文内処理の進捗確認
#install.packages("tcltk")
library(tcltk)
#minでループ変数の最小値、maxでループ変数の最大値、styleで表示スタイル(基本は3)を設定
pb <- txtProgressBar(min=1, max=B, style=3)

for (b in 1:B) { 
     bt <- sample(1:n,replace = TRUE) #ブートストラップ標本番号 #復元抽出
     mydata <- dataX[bt,] #ブートストラップ標本
     output3 <- highDLingam::findGraphMulti(mydata, maxInDegree = md, cutOffScaling = cs, degree = deg, verbose = F)
     adj_mat <- matrix(0,dim(dataX)[2],dim(dataX)[2]) #有向辺の有無（1 or 0）の途中結果格納用を初期化
     for(i in 1:length(output3$topOrder)){
          k <- output3$topOrder[i]
          for (j in output3$parents[i]){
               adj_mat[k,j] <- 1
          }
     } #create the adjacency matrix
     adj_mat_boot <- adj_mat_boot + adj_mat #最終結果の格納用に、有効辺の有無（1 or 0）を代入し更新する
     setTxtProgressBar(pb, b) 
}
adj_mat_boot #B回繰り返した時の、有向辺の有無（1 or 0）の最終結果を表示

adj_mat_boot_prob <- adj_mat_boot / B 
adj_mat_boot_prob #B回繰り返した時のブートストラップ確率を表示

### Save the results
adj_mat_boot <- data.frame(adj_mat_boot)
write.csv(adj_mat_boot,"adj_mat_boot.csv",row.names=TRUE)

adj_mat_boot_prob <- data.frame(adj_mat_boot_prob) #ブートストラップ確率
write.csv(adj_mat_boot_prob,"adj_mat_boot_prob.csv",row.names=TRUE)
