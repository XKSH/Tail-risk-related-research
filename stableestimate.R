require(StableEstim)
theta <- c(1.45,0.55,1,0)
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
McCullochParametersEstim(x+y)
McCullochParametersEstim(a*b)
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