#
# Test LV dynamics
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
#
#

