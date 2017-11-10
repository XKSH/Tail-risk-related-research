require(StableEstim)
theta <- c(1.45,0,1,0)
#theta <- c(2,0,1,0)
pm <- 0
# set.seed(2345)
# x <- rstable(200,theta[1],theta[2],theta[3],theta[4],pm)
# set.seed(1234)
# y<- rstable(200,theta[1],theta[2],theta[3],theta[4],pm)
set.seed(2345)
rxy=rstable(80,theta[1],theta[2],theta[3],theta[4],pm)
x <-rxy[1:40]
y=rxy[41:80]
a <-0.5*x+y;b=-2*x+y
#Kout bad; should with adjustment of initialisation
objKout <- Estim(EstimMethod="Kout",data=x,pm=pm, ComputeCov=FALSE,HandleError=FALSE,spacing="Kout")
objKout@par
#ML bad and long
ML <- MLParametersEstim(x=x,pm=pm,PrintTime=TRUE)
ML$Estim$par
#McCulloch good
#principle is simple use x~S(a,b,1,0);y=cx+mu;y~S(a,b,c,mu)
McCullochParametersEstim(x+y)
vtest=McCullochParametersEstim(a*b)
# Kogon regression with the McCulloch initialisation bad and long for the reason of Kogon regression
IGParametersEstim(x,pm=0)

#GMM bad need to tuning
regularization="cut-off"
WeightingMatrix="OptAsym"
alphaReg=0.005
t_seq=seq(0.1,2,length.out=12)

## If you are just interested by the value
## of the 4 estimated parameters
t_scheme="free"
algo = "2SGMM"

suppressWarnings(GMMParametersEstim(x=x,
                                    algo=algo,alphaReg=alphaReg,
                                    regularization=regularization,
                                    WeightingMatrix=WeightingMatrix,
                                    t_scheme=t_scheme,
                                    pm=pm,PrintTime=TRUE,t_free=t_seq))
#Cgmm good
alphaReg=0.01
subdivisions=20
randomIntegrationLaw="unif"
IntegrationMethod="Uniform"

## Estimation
twoS <- CgmmParametersEstim(x=x,type="2S",alphaReg=alphaReg,
                            subdivisions=subdivisions,
                            IntegrationMethod=IntegrationMethod,
                            randomIntegrationLaw=randomIntegrationLaw,
                            s_min=0,s_max=1,theta0=NULL,
                            pm=pm,PrintTime=TRUE)
twoS$Estim$par


#svd
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
X <- hilbert(9)[, 1:6]
s <- svd(X)
D <- diag(s$d)
s$u %*% D %*% t(s$v)
#x is the vector of parameters(alpha(a),beta(b),gamma(c),delta(mu))
tail.amplitude=function(x){
  amph=x[3]*(sin(pi*x[1]/2)*gamma(x[1])/pi)^(1/x[1])
  return(as.numeric(amph))
}
tail.amplitude(vtest)

#simulation 
pm = 0
theta <- c(1.45,0,1,14) 
set.seed(2345) 
nb=1200;M=2;nv=3 
set.seed(2345) 
a=rstable(nb,theta[1],theta[2],theta[3],theta[4],pm) 
theta <- c(1.45,0,2,1) 
set.seed(1234) 
b=rstable(nb,theta[1],theta[2],theta[3],theta[4],pm) 
theta <- c(1.45,-0.8,6,-14) 
set.seed(1234)
c=rstable(nb,theta[1],theta[2],theta[3],theta[4],pm)
# rab=rstable(nb,theta[1],theta[2],theta[3],theta[4],pm)
# ind=seq(1,nb,nb/M)
# a =rab[ind[1]:(nb/M)]
# b=rab[ind[2]:(2*nb/M)]
x<-0.5*a+0.5*b+0.3*c;y=-2*a+b+0.6*c;z=0.4*a+0.7*b+0.6*c
d <- density(x+y+z)
plot(d, main="Kernel Density of stable distributions")
polygon(d, col="red", border="blue") 



C.matrix=matrix(0,nrow=nv,ncol=nv)
data=rbind(x,y,z)
eta=matrix(0,nrow=(nv*(nv+1)/2),ncol=(nb))
idx=0
for(i in 1:nv){
  for(j in i:nv){
    idx=idx+1
    eta[idx,]=data[i,]*data[j,]
  }
}
#we can change sign of s1$u according to s$u,vice versa
spara=apply(eta,1,McCullochParametersEstim)
tcor=apply(spara,2,function(x){tail.amplitude(x)})
tcor=spara[2,]*tcor^spara[1,]
C.matrix[lower.tri(C.matrix, diag=TRUE)] <- tcor
C.matrix[upper.tri(C.matrix)] <- t(C.matrix)[upper.tri(C.matrix)]

A.matrix=matrix(0,nrow=nv,ncol=nv)
A.matrix[lower.tri(A.matrix, diag=TRUE)]=1/spara[2,]*tcor
A.matrix[upper.tri(A.matrix)] <- t(A.matrix)[upper.tri(A.matrix)]
s1=eigen(A.matrix)
A1=diag(s1$values)
s1$vectors
s <- eigen(C.matrix)
A <- diag(s$values)
s$vectors %*% A %*% t(s$vectors)
# tail.amplitude(c(1.45,0,1,0))
# tail.amplitude(c(1.45,0,0.2,1))
# tail.amplitude(c(1.45,0,0.3,4))
#regulariation

C.mod=nearPD(C.matrix, corr = FALSE,ensureSymmetry=TRUE,conv.norm.type = "F",maxit=200)
s.mod=eigen(C.mod$mat)
s.v=s.mod$vectors
sign.v=s.v/abs(s.v)

# abs(v)%*%diag(s.mod$values)%*%abs(v)
#parameters of factors a,b,c
alpha=1.96#fixed before
v=sign.v*abs(s.v)^(2/alpha)
rpara=apply(rbind(x,y,z),1,McCullochParametersEstim)
mu=solve(v)%*%rpara[4,]
# c=(solve(abs(v)^(alpha))%*%(rpara[3,]^alpha))^(1/alpha)#c may be calculated negative so change it to constrained optimisation
require(limSolve)#optimization using limsolve
ineq=as.matrix(rbind(diag(3)))
c=lsei(abs(v)^(alpha),rpara[3,]^alpha,G=ineq,H=rep(0,3),type=2)$X
cinv=t(matrix(rep(c,3),nrow=3))
cinv=abs(v)*cinv
cinv=cinv^alpha
cinv=cinv/rowSums(cinv)
# b=solve(cinv)%*%rpara[2,] normally it should be constrained b in -1 to 1
ineq=as.matrix(rbind(diag(3),-diag(3)))
b=lsei(cinv,rpara[2,],G=ineq,H=rep(-1,6),type=2)$X
b[b>1]=1
b[b<-1]=1
fpara=cbind(rep(alpha,3),b,c,mu)

# #constrOptim
# d = as.vector(rpara[2,])
# D= cinv
# start <- rep(1/(ncol(cinv)+1), ncol(cinv))
# min_fn <- function(b, dvec, Dmat){t(Dmat%*%b-dvec)%*%(Dmat%*%b-dvec)/2}
# grad_min_fn <- function(b, dvec, Dmat){(Dmat%*%b-dvec)%*%t(b)}
# cond=rep(-1,6)
# 
# constrOptim(theta=start, f=min_fn, grad=grad_min_fn, ui=ineq, ci=cond, control=list(reltol=10*.Machine$double.eps),
#             dvec=d, Dmat=D )


random=matrix(0,nrow=3,ncol=6000)
for(i in 1:3){
  set.seed(2345)
  random[i,]=rstable(dim(random)[2],fpara[i,1],fpara[i,2],fpara[i,3],fpara[i,4],pm)
}
rport=v%*%random
# rport=v[,1]%*%t(random[1,])
# rport=v[,c(1,3)]%*%random[c(1,3),]
d <- density(colSums(rport))
plot(d, main="Kernel Density of stable distributions")
polygon(d, col="red", border="blue") 

#plot of reconstructed individuals
library(ggplot2)
#Sample data
sample=as.vector(cbind(x,y,z))
dat <- data.frame(dens = sample
                  , individuals = rep(c("x", "y","z"), each = 1200))
#
ggplot(dat, aes(x = dens,fill = individuals,colour = individuals)) + geom_density(adjust = 3,alpha = 0.5)+
  xlim(-100, 50)

#Sample data
sample=as.vector(t(rport))
dat <- data.frame(dens =sample
                  , individuals = rep(c("x", "y","z"), each = 6000))
#
ggplot(dat, aes(x = dens,fill = individuals,colour = individuals)) + geom_density(adjust =3,alpha = 0.5)+
  xlim(-100, 50)

#plot of factors
#Sample data
sample=as.vector(t(random))
dat <- data.frame(dens = sample
                  , factors = rep(c("a", "b","c"), each = 6000))
#Plot of factors
ggplot(dat, aes(x = dens,fill = factors,colour = factors)) + geom_density(adjust = 10,alpha = 0.5)+
  xlim(-100, 50)


#reconstruire A matrix
fcor=apply(fpara,1,function(x){tail.amplitude(x)})
A.left=abs(v)^(alpha/2)%*%diag(fcor)^alpha%*%t(abs(v)^(alpha/2))

###use last ten quantile to identify best alpha
ms=x+y+z
mr=colSums(rport)
plot(sort(ms)[seq(2,120,1)],sort(mr)[seq(10,600,5)])
abline(0, 1)
plot(sort(z)[seq(12,120,12)],sort(rport[3,])[seq(60,600,60)])
abline(0, 1)


max(colSums(rport));min(colSums(rport))
max(x+y+z);min(x+y+z)
A.mod=nearPD(A.matrix, corr = FALSE,conv.norm.type = "F")
eigen(A.mod$mat)$vectors

#shrinkage
A=A[1:2,1:2]
u=s$u[,1:2]
v=s$v[,1:2]
u %*% A %*% t(v)
tt=matrix(c(-0.2223282,0.9550149,-0.1962567,-0.97479174,-0.21386781,0.06357381,0.01874094,0.20544367,0.97848949),nrow=3)


#normal distribution
xtest=x+y+z
exp=function(x) dnorm(x,mean=mean(xtest),sd=sd(xtest))
curve(exp, -200, 200, xname = "t")
qnorm(0.001,mean=mean(xtest),sd=sd(xtest))
#whether quantile instead of sample follow same distribution
#put more weight in tail
require(kSamples)
ad.test(sort(ms)[seq(12,1200,12)],sort(mr)[seq(60,6000,60)])
ks.test(sort(ms)[seq(12,1200,12)],sort(mr)[seq(60,6000,60)])
ks.test(sort(ms),sort(mr))

ad.test(sort(z)[seq(12,1200,12)],sort(rport[3,])[seq(60,6000,60)])
ks.test(sort(z)[seq(12,1200,12)],sort(rport[3,])[seq(60,6000,60)])
ks.test(sort(x),sort(rport[1,]))

ad.test(rnorm(30),runif(30))
ad.test(sort(rnorm(30))[seq(3,30,3)],sort(runif(30))[seq(3,30,3)])
ad.test(sort(rnorm(30,0,1)),sort(rnorm(30,0,2)))
ad.test(sort(rnorm(300,0,1))[seq(10,300,10)],sort(rnorm(300,0,2))[seq(10,300,10)])
