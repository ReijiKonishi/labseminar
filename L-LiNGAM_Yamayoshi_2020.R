
# 目的 ----------------------------------------------------------------------

# Yamayoshi_2020のL-LiNGAMのNumerical Studyで書かれているデータ生成過程を試してみる

# -------------------------------------------------------------------------

library(dplyr)

# 関数定義 --------------------------------------------------------------------

logistic <- function(x) {
  1 / (1 + exp(x))
}

indicator <- function(score) {
  dplyr::if_else(score < 0.5, 0, 1)
}


# データ生成 -------------------------------------------------------------------


set.seed(123)
sample_size <- 1000
e1 <- rt(sample_size, df = 10)
e2 <- rt(sample_size, df = 10)
e3 <- rt(sample_size, df = 10)
e4 <- rt(sample_size, df = 10)

f4 <- rt(sample_size, df = 10)
f3 <- 2 * f4
f2 <- 2*f4 + 2*f3
f1 <- 2*f4 + 2*f3 + 2*f2

x1 <- f1 + e1

score2 <- logistic(-f2 + e2)

x3 <- f3 + e3

score4 <- logistic(-f4 + e4)

dat <- tibble(x1 = x1,
              score2 = score2,
              x3 = x3,
              score4 = score4) %>% 
  mutate(x2 = indicator(score2),
         x4 = indicator(score4)) %>% 
  select(x1, x2, x3, x4)

