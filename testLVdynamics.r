#
# Test package against LV dynamics 
#

require(deSolve)
require(ggplot2)
require(tidyr)
require(tibble)
# 
#


parms <- tibble::tribble(
  ~a1,   ~a2,  ~a3,    ~a4,   ~a5,  ~a6,    ~r,
    0,  4e-6, 3e-5,   2e-5,     0,    0,  -0.01,
-1e-5, -1e-4,    0,   1e-4,  6e-5,    0, -0.001,
-1e-4,     0,    0,      0,  6e-5, 5e-5,  -0.01,
-1e-4,  1e-4,    0,  -5e-4, -1e-5,    0,   0.08,
    0, -1e-4,-1e-4,      0, -1e-6,    0,   0.04,
    0,     0,-1e-4,   2e-4,     0,    0,   0.00
)

# Solve for equilibrium populations
#
eq <- lm(r ~ . + 0, parms)
-coef(eq)



# Generalized Lotka-Volterra 6 species  equation 
#
GLV_det<-function(Time,State,Pars){
  with(as.list(c(State, Pars)), {
    A <- matrix(Pars,6,7)
    dP1 <- N1*(A[1,7]+A[1,1]*N1+A[1,2]*N2+A[1,3]*N3+A[1,4]*N4+A[1,5]*N5+A[1,6]*N6) 
    dP2 <- N2*(A[2,7]+A[2,1]*N1+A[2,2]*N2+A[2,3]*N3+A[2,4]*N4+A[2,5]*N5+A[2,6]*N6)
    dP3 <- N3*(A[3,7]+A[3,1]*N1+A[3,2]*N2+A[3,3]*N3+A[3,4]*N4+A[3,5]*N5+A[3,6]*N6)
    dP4 <- N4*(A[4,7]+A[4,1]*N1+A[4,2]*N2+A[4,3]*N3+A[4,4]*N4+A[4,5]*N5+A[4,6]*N6)
    dP5 <- N5*(A[5,7]+A[5,1]*N1+A[5,2]*N2+A[5,3]*N3+A[5,4]*N4+A[5,5]*N5+A[5,6]*N6)
    dP6 <- N6*(A[6,7]+A[6,1]*N1+A[6,2]*N2+A[6,3]*N3+A[6,4]*N4+A[6,5]*N5+A[6,6]*N6)
    
    return(list(c(dP1,dP2,dP3,dP4,dP5,dP6)))
  })
}

# Simulate the model by integration
#
pars  <- c(A <- as.matrix(parms[,1:6]),
           r<- as.numeric(parms$r)
)   

yini  <- c(N1 = 100,N2=100,N3=100,N4=100,N5=100,N6=100)
out   <- as.data.frame(ode(yini, 1:600, GLV_det, pars))
df <- gather(out,key="Species",value="N", -time)
ggplot(df, aes(time,N,colour=Species)) + geom_line(size=.5)

A <- as.matrix(parms[,1:6])
r<- as.numeric(parms$r)
m <- c(0,0,0,0,0,0)
#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A1 <- metaWebNetAssemblyGLV(A,m,r,yini,1000,0.1)
df1 <- as.data.frame(t(A1$STime))
df1$time <- 1:1000
A1 <- metaWebNetAssemblyGLV(A,m,r,yini,1000,0.1)
df2 <- as.data.frame(t(A1$STime))
df2$time <- 1:1000

df1 <- rbind(df1,df2 )


names(df1)[1:6] <- unique(df$Species) 
require(tidyr)
df1 <- gather(df1,key="Species",value="N", N1:N6)
require(ggplot2)
ggplot(df, aes(time,N,colour=Species)) + geom_line() + geom_point(data=df1,aes(time,N,colour=Species),size=0.1) + theme_bw() + scale_color_viridis_d()

ggplot(df1, aes(time,N,colour=Species)) + geom_point(size=0.1) +  theme_bw() + scale_color_viridis_d()

#
# LV with migration
#

parms <- tibble::tribble(
  ~a1,   ~a2,  ~a3,    ~a4,   ~a5,  ~a6,    ~r,   ~m,
    0,  4e-6, 3e-5,   2e-5,     0,    0,  -0.01, 0.1,
-1e-5, -1e-4,    0,   1e-4,  6e-5,    0, -0.001, 0.1,
-1e-4,     0,    0,      0,  6e-5, 5e-5,  -0.01, 0.1,
-1e-4,  1e-4,    0,  -5e-4, -1e-5,    0,   0.08, 0.1,
    0, -1e-4,-1e-4,      0, -1e-6,    0,   0.04, 0.1,
    0,     0,-1e-4,   2e-4,     0,    0,   0.00, 0.1
)

# Generalized Lotka-Volterra with migration 6 species  equation 
#
GLVM_det<-function(Time,State,Pars){
  with(as.list(c(State, Pars)), {
    A <- matrix(Pars,6,8)
    dP1 <- N1*(A[1,7]+A[1,1]*N1+A[1,2]*N2+A[1,3]*N3+A[1,4]*N4+A[1,5]*N5+A[1,6]*N6)+A[1,8] 
    dP2 <- N2*(A[2,7]+A[2,1]*N1+A[2,2]*N2+A[2,3]*N3+A[2,4]*N4+A[2,5]*N5+A[2,6]*N6)+A[2,8]
    dP3 <- N3*(A[3,7]+A[3,1]*N1+A[3,2]*N2+A[3,3]*N3+A[3,4]*N4+A[3,5]*N5+A[3,6]*N6)+A[3,8]
    dP4 <- N4*(A[4,7]+A[4,1]*N1+A[4,2]*N2+A[4,3]*N3+A[4,4]*N4+A[4,5]*N5+A[4,6]*N6)+A[4,8]
    dP5 <- N5*(A[5,7]+A[5,1]*N1+A[5,2]*N2+A[5,3]*N3+A[5,4]*N4+A[5,5]*N5+A[5,6]*N6)+A[5,8]
    dP6 <- N6*(A[6,7]+A[6,1]*N1+A[6,2]*N2+A[6,3]*N3+A[6,4]*N4+A[6,5]*N5+A[6,6]*N6)+A[6,8]
    
    return(list(c(dP1,dP2,dP3,dP4,dP5,dP6)))
  })
}

pars <- as.matrix(parms)
yini  <- c(N1 = 0,N2=0,N3=0,N4=0,N5=0,N6=0)
out   <- as.data.frame(ode(yini, 1:600, GLVM_det, pars))
df <- gather(out,key="Species",value="N", -time)
ggplot(df, aes(time,N,colour=Species)) + geom_line(size=.5)
calcPropInteractionsGLVadjMat(A1$interM,A2$STime[,600])


A <- as.matrix(parms[,1:6])
r<- as.numeric(parms$r)
m <-as.numeric(parms$m)
#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A1 <- metaWebNetAssemblyGLV(A,m,r,yini,600,0.1)
df1 <- as.data.frame(t(A1$STime))
df1$time <- 1:600
A1 <- metaWebNetAssemblyGLV(A,m,r,yini,600,0.1)
df2 <- as.data.frame(t(A1$STime))
df2$time <- 1:600

df1 <- rbind(df1,df2 )


names(df1)[1:6] <- unique(df$Species) 
require(tidyr)
df1 <- gather(df1,key="Species",value="N", N1:N6)
require(ggplot2)
ggplot(df, aes(time,N,colour=Species)) + geom_line() + geom_point(data=df1,aes(time,N,colour=Species),size=0.1) + theme_bw() + scale_color_viridis_d()

#
# Generate random adjacency matrix 
#
require(igraph)

#### ER random graph
build.random.structure <- function(n, connectance){
  S <- matrix(0, nrow = n, ncol = n)
  S[upper.tri(S)] <- rbinom(n * (n - 1)/2,size = 1, p = connectance)
  S[lower.tri(S)] <- rbinom(n * (n - 1)/2,size = 1, p = connectance)
  #S <- S + t(S)
  #diag(S) <- 0
  return(S)    
}


#
# Size 106, links 4623
#
require(multiweb)

C <- 4623/(106*106)

#
# Generate random GLV matrix
#

A <- generateRandomGLVadjMat(100,0.01,c(0.3,0.3,0.3,0.05,0.05)) 
A1 <- generateGLVparmsFromAdj(A,0.01,0.01)
yini <- rep(10,times=nrow(A1$interM))
calcPropInteractionsGLVadjMat(A ,yini)

yini <- rep(0,times=nrow(A1$interM))
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,yini,600,0.1)
A2$STime[,600]
df1 <- as.data.frame(t(A2$STime))
df1$time <- 1:600
require(tidyr)
require(ggplot2)
require(dplyr)
df1 <- gather(df1,key="Species",value="N", -time)
ggplot(filter(df1,time>500), aes(time,N,colour=Species)) + geom_line() +  theme_bw() + scale_color_viridis_d(guide=FALSE) 
ggplot(df1, aes(time,N,colour=Species)) + geom_line() +  theme_bw() + scale_color_viridis_d(guide=FALSE) 

df1 <- data_frame(time=1:600,S = A2$S,C = A2$L/(A2$S*A2$S))
ggplot(filter(df1,time>500), aes(time,S)) + geom_line() +  theme_bw() + scale_color_viridis_d() 
ggplot(df1, aes(time,S)) + geom_line() +  theme_bw() + scale_color_viridis_d() 

ggplot(filter(df1,time>500), aes(time,C)) + geom_line() +  theme_bw() + scale_color_viridis_d() 
ggplot(df1, aes(time,C)) + geom_line() +  theme_bw() + scale_color_viridis_d() 
ggplot(df1, aes(S,C)) + geom_line() +  theme_bw() + scale_color_viridis_d() 

calcPropInteractionsGLVadjMat(A1$interM,A2$STime[,600])

# With no migration 
yini <- rep(100,times=nrow(A1$interM))
A1$m <- rep(0.0, 100)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,yini,600,0.1)
A2$STime[,600]
df1 <- as.data.frame(t(A2$STime))
df1$time <- 1:600
require(tidyr)
require(ggplot2)
require(dplyr)
df1 <- gather(df1,key="Species",value="N", -time)
ggplot(df1, aes(time,N,colour=Species)) + geom_line() +  theme_bw() + scale_color_viridis_d(guide=FALSE) 

calcPropInteractionsGLVadjMat(A1$interM,A2$STime[,600])


#
# To igraph
#
require(multiweb)
require(igraph)
diag(A2$A) <- 0
g <- graph_from_adjacency_matrix(A2$A,mode="directed")
plotTrophLevel(g,modules=TRUE)
calcTopologicalIndices(g)

#
# More testing
#

A <- generateRandomGLVadjMat(100,0.05,c(0.2,0.5,0.1,0.0,0.0)) 
A
A1 <- generateGLVparmsFromAdj(A,0.1,0.01)

yini <- rep(1,times=nrow(A1$interM))

calcPropInteractionsGLVadjMat(A1$interM,yini)

A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,yini,100,0.1)

df <- as.data.frame(t(A2$STime))
df$Time <- 1:100
require(tidyr)
df <- gather(df,key="Species",value="N", V1:V100)
require(ggplot2)
ggplot(df, aes(Time,N,colour=Species)) + geom_line() +guides(colour=FALSE)

calcPropInteractionsGLVadjMat(A1$interM,A2$STime[,100])

# Delete diagonals
diag(A2$A) <- FALSE
# Delete columns and rows with no species
d <- A2$STime[,100]!=0
A <- A2$A[d,d]
# Size of d == Number of species
sum(d)==A2$S[100]

g <- graph_from_adjacency_matrix(A,mode="directed")
g
plotTrophLevel(g)
calcTopologicalIndices(g) # Trophic level do not work with multiple interactions types

require(NetIndices)
TrophInd(A2$A)
