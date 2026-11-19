
set.seed(11)
shap <- (0.01/0.001)^2
shap

rate <- 0.01/(0.001)^2
rnd <- rgamma(761,shap,rate)
mean(rnd)
# [1] 0.01002665
sd(rnd)
# [1] 0.001037943

A <- generateRandomGLVadjMat(200,0.05,c(0.3,0.3,0.3,0.05,0.05)) 
set.seed(11)
A1 <- generateGLVparmsFromAdj(A,0.1,0.01,preserveInt = 0,predIntSd=0.001,rndType = 0)
set.seed(11)
A2 <- generateGLVparmsFromAdjR(A,0.1,0.01,preserveInt = 0,predIntSd=0.001,rndType = 0)

all.equal(A2,A1)

diag(A1$interM) <- 0
ii <- abs(A1$interM[A1$interM != 0])
mean(ii)
sd(ii)

ii <- read.csv("test.txt",header=TRUE)
ii <- abs(ii$empty.rndgamma)
