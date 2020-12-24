# ���O�p�s�񂩂���
adj_mat <- matrix(0,100,100)
for (i in 2:nrow(adj_mat)) {
     for (j in 1:(i-1)) {
          adj_mat[i,j] <- (runif(1,0,1) > 0.95)*1 #��l�����𔭐�������B5���̊m����1�𓾂�B
     }
}
adj_mat #100*100�̉��O�p�s��

# ��
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


# adj_mat����v�f��1�ł���s�E��ԍ����m�F
for (i in 2:nrow(adj_mat)) {
     for (j in 1:(i-1)) {
          #print(adj_mat[i,j])
          if (adj_mat[i,j] == 1){
               print(c(i,j))
          }
     }
}

�@
# �W���s��̍쐬
coef_mat <- matrix(0,nrow(adj_mat),ncol(adj_mat))
for (j in 1:ncol(adj_mat)) {
     for (i in 2:nrow(adj_mat)) {
          coef_mat[i,j] <- sample(2:5,1) * adj_mat[i,j] #adj_mat:�אڍs��̊e�v�f�ɁA2~5�̂����ꂩ�̐������|����B
     }
}
coef_mat
rownames(coef_mat) <- colname


# x_parent_sum �����
# x�ɐe�����݂���ꍇ�A���̐e�̔ԍ��̏��𒊏o����B
# �Ⴆ�΁Ax2 = 3 * x6 + 2 * x17 �ł���΁Ax_parent_num[2,6] = 6, x_parent_num[2,17] = 17 ������B
x_parent_num <- matrix(0,nrow(coef_mat),ncol(coef_mat))
for (i in 1:nrow(coef_mat)) {
     for (j in 1:ncol(coef_mat)) {
          if (coef_mat[i,j] != 0){
               x_parent_num[i,j] <- j
          }
     }
}
x_parent_num


# D:Test data �̍쐬�i80*100�j
set.seed(123)
n_samples <- 80
D <- matrix(0,n_samples,ncol(coef_mat))
colnames(D) <- colname


z <- numeric(0) #�e�ϐ��̐���������z(���̐���)
parent_sum_mat <- matrix(0,n_samples,1) #�e�ϐ��̓��Ϙa���i�[����s��i80*1�j
for (i in 1:nrow(coef_mat)) {
     #for (j in 1:ncol(coef_mat)) {
     if(sum(coef_mat[i,] != 0) == 0){ #�W���s���i�s�ڂ����ׂ�0�ł����
          D[,i] <- runif(n=n_samples) #xi�ɐe�ϐ��͑��݂��Ȃ��̂�,��l�����𐶐����A���
     }
     else{ #�W���s���i�s�ڂ����ׂ�0�łȂ���΁Axi�ɐe�ϐ������݂���̂�
          coef <- subset(coef_mat[i,],coef_mat[i,] != 0) #�e�ϐ��̓��Ϙa���v�Z���邽�߂ɁA�W���𒊏o
          parent_sum_mat <- matrix(0,n_samples,1)  #�e�ϐ��̓��Ϙa���i�[����s��i80*1�j
          z <- numeric(0) #�e�ϐ��̐���������z(���̐���)
          for (z in 1:length(coef)) { 
               parent_sum_mat <- parent_sum_mat + ( coef[z]*D[, subset(x_parent_num[i,],x_parent_num[i,] != 0)[z]] ) #D[,z]��z�������� -> x_parent_num[i,j]
          } #parent_sum_mat�ɁAxi�̐e�ϐ��̓��Ϙa�������A�X�V���Ă���
          D[,i] <- runif(n=n_samples) + parent_sum_mat #��l�����ɐe�ϐ��̒l�������āA�f�[�^�𐶐�
     }
     #}
}

D
X <- D
write.csv(X,"sampledata_200922.csv",row.names=FALSE)