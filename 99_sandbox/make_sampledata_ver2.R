# 下三角行列から作る
adj_mat <- matrix(0,100,100)
for (i in 2:nrow(adj_mat)) {
     for (j in 1:(i-1)) {
          adj_mat[i,j] <- (runif(1,0,1) > 0.95)*1 #一様乱数を発生させる。5％の確率で1を得る。
     }
}
adj_mat #100*100の下三角行列

# 列名
colname <- c()
for (i in 1:ncol(adj_mat)) {
     colname <- c(colname,paste0("x",i))
}
colname
colnames(adj_mat) <- colname;rownames(adj_mat) <- colname

#adj_mat <- data.frame(adj_mat)
#write.csv(adj_mat,"adj_mat.csv",row.names=FALSE)

# Draw the causal graph
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


# adj_matから要素が1である行・列番号を確認
for (i in 2:nrow(adj_mat)) {
     for (j in 1:(i-1)) {
          #print(adj_mat[i,j])
          if (adj_mat[i,j] == 1){
               print(c(i,j))
          }
     }
}

　
# 係数行列の作成
coef_mat <- matrix(0,nrow(adj_mat),ncol(adj_mat))
for (j in 1:ncol(adj_mat)) {
     for (i in 2:nrow(adj_mat)) {
          coef_mat[i,j] <- sample(2:5,1) * adj_mat[i,j] #adj_mat:隣接行列の各要素に、2~5のいずれかの整数を掛ける。
     }
}
coef_mat
rownames(coef_mat) <- colname


# x_parent_sum を作る
# xに親が存在する場合、その親の番号の情報を抽出する。
# 例えば、x2 = 3 * x6 + 2 * x17 であれば、x_parent_num[2,6] = 6, x_parent_num[2,17] = 17 が入る。
x_parent_num <- matrix(0,nrow(coef_mat),ncol(coef_mat))
for (i in 1:nrow(coef_mat)) {
     for (j in 1:ncol(coef_mat)) {
          if (coef_mat[i,j] != 0){
               x_parent_num[i,j] <- j
          }
     }
}
x_parent_num


# D:Test data の作成（80*100）
set.seed(123)
n_samples <- 80
D <- matrix(0,n_samples,ncol(coef_mat))
colnames(D) <- colname


z <- numeric(0) #親変数の数を代入するz(正の整数)
parent_sum_mat <- matrix(0,n_samples,1) #親変数の内積和を格納する行列（80*1）
for (i in 1:nrow(coef_mat)) {
     #for (j in 1:ncol(coef_mat)) {
     if(sum(coef_mat[i,] != 0) == 0){ #係数行列のi行目がすべて0であれば
          D[,i] <- runif(n=n_samples) #xiに親変数は存在しないので,一様乱数を生成し、代入
     }
     else{ #係数行列のi行目がすべて0でなければ、xiに親変数が存在するので
          coef <- subset(coef_mat[i,],coef_mat[i,] != 0) #親変数の内積和を計算するために、係数を抽出
          parent_sum_mat <- matrix(0,n_samples,1)  #親変数の内積和を格納する行列（80*1）
          z <- numeric(0) #親変数の数を代入するz(正の整数)
          for (z in 1:length(coef)) { 
               parent_sum_mat <- parent_sum_mat + ( coef[z]*D[, subset(x_parent_num[i,],x_parent_num[i,] != 0)[z]] ) #D[,z]のzがちがう -> x_parent_num[i,j]
          } #parent_sum_matに、xiの親変数の内積和を代入し、更新していく
          D[,i] <- runif(n=n_samples) + parent_sum_mat #一様乱数に親変数の値を加えて、データを生成
     }
     #}
}

D
X <- D
write.csv(X,"sampledata_200922.csv",row.names=FALSE)
