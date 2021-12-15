library(RSienaTest)

mynet1 <- sienaDependent(array(c(s501,s502,s503), dim=c(50, 50, 3)))
mynet2 <- sienaDependent(s50a, type='behavior')
mydata <- sienaDataCreate(mynet1, mynet2)
#mynet1 <- sienaDependent(array(c(s501,s502), dim=c(50, 50,2)))
#mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, recip, include=TRUE)
myeff <- includeEffects(myeff, transTrip, inPopSqrt)
myeff <- includeEffects(myeff, altX, interaction1='mynet2')
myeff <- includeEffects(myeff, outRateInv, type="rate", include=TRUE)

mynetm1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32,32,2)))
mydatam <- sienaDataCreate(mynetm1)
myeffm <- getEffects(mydatam)
myeffm <- includeEffects(myeffm, transTrip, inPopSqrt)
MOMmodel <- sienaModelCreate(cond=FALSE, maxlike=FALSE, n3=100, seed=1)
MLmodel <- sienaModelCreate(cond=FALSE,  maxlike=TRUE, n3=100, seed=1)

ans <- siena07(MOMmodel, data=mydata, effects=myeff, useCluster=TRUE,
                verbose=TRUE)
ans50ml <- siena07(MLmodel, data=mydata, effects=myeff, useCluster=FALSE,
                verbose=TRUE)
ansm <- siena07(MOMmodel, data=mydatam, effects=myeffm, useCluster=TRUE,
               verbose=TRUE)
ansml <- siena07(MLmodel, data=mydatam, effects=myeffm, useCluster=FALSE,
                verbose=TRUE)

## MOM statistics / Gelman few/ few many/all
print("resp1")
system.time(resp1 <- algorithms(mydata, myeff, MOMmodel, nIter=20, numiter=20,
                    nbrNodes=2))
myeff <- includeEffects(myeff, outRateInv, type="rate", include=FALSE)

## ML scores / Gelman few/ few many/all
print("resp2")
system.time(resp2 <- algorithms(mydata, myeff, MLmodel, nIter=20, numiter=20, nbrNodes=2))

## SEM + Gelman
print("resp3")
system.time(resp3 <- algorithms(mydata, myeff, MLmodel, finalIter=50,
                    useOptim=TRUE, optimFinal=1, optimSchedule=rep(50, 10),
                    nbrNodes=2))
## EM varied sampling rate
print("resp4")
system.time(resp4 <- algorithms(mydata, myeff, MLmodel,  useOptim=TRUE,
                    optimFinal=1,
                    scale=1, nbrNodes=2))
## EM reuse history
print("resp5")
system.time(resp5 <- algorithms(mydata, myeff, MLmodel,  useOptim=TRUE,
                    optimFinal=1, useHistory=TRUE, optimWeight=0.25,
                    scale=1, nbrNodes=2))

## Geyer in final loop
print("resp6")
system.time(resp6<-algorithms(mydata, myeff, MLmodel,  useOptim=TRUE,
                    optimFinal=2, scale=0.5, nbrNodes=2))

## Geyer with varied theta
print("resp7")
system.time(resp7 <- algorithms(mydata, myeff, MLmodel,  useOptim=TRUE,
                    optimFinal=2, variedThetas=TRUE, scale=0.5,
                    finalIter=100, optimSchedule=c(10, 20, 20, 20, 50,50),
                    nbrNodes=2))


## response surface
print("resp8")
system.time(resp8 <- algorithms(mydata, myeff, MOMmodel, nIter=100, numiter=10,
                    finalIter=100, responseSurface=TRUE, scale=1,
                    optimWeight=0.75, nbrNodes=2))

## profileLikelihoods

profileLikelihoods(list(theta=c(5.2,5.2,-2.5, 2, .7, 0.05, .1,1, 1, .5, .5)),
system.time(RSiena:::profileLikelihoods(resp8, MLmodel, mydata, myeff,
                   1,3, gridl=c(0.9,1.1),nIter=100, nbrNodes=2))
system.time(RSienaTest:::profileLikelihoods(resp8, MLmodel, mydata, myeff,
                   1,3, gridl=c(0.9,1.1),nIter=100, nbrNodes=2))
RSienaTest:::profileLikelihoods(list(theta=c(5.2,5,-2.5, 2, .7, 0.05, .1, 1, 1, .5, .5)),
MLmodel, mydata, myeff,
                   1,3, gridl=c(0.9,1.1),nIter=100, nbrNodes=1)

##bayes
print("resp9")
system.time(resp9 <- sienaBayes(mydata, myeff, MLmodel, nwarm=100,nmain=100, nbrNodes=2))


system.time(resp1 <- algorithms(mydatam, myeffm, MOMmodel, nIter=20, numiter=20,
                    nbrNodes=1))

## ML scores / Gelman few/ few many/all
print("resp2")
system.time(resp2 <- algorithms(mydatam, myeffm, MLmodel, nIter=20, numiter=20,
                    nbrNodes=1))

## SEM + Gelman
print("resp3")
system.time(resp3 <- algorithms(mydatam, myeffm, MLmodel,  finalIter=50,
                    useOptim=TRUE, optimFinal=1, optimSchedule=rep(50, 10),
                    nbrNodes=1))
## EM varied sampling rate
print("resp4")
system.time(resp4 <- algorithms(mydatam, myeffm, MLmodel,  useOptim=TRUE,
                    optimFinal=1, optimSchedule=c(10, 20, 20, 20, 50),
                    scale=1, nbrNodes=1))
## EM reuse history
print("resp5")
system.time(resp5 <- algorithms(mydatam, myeffm, MLmodel,  useOptim=TRUE,
								optimFinal=1, useHistory=TRUE, scale=0.75,
								optimSchedule=c(10,20,20,20,50)))
## Geyer in final loop
print("resp6")
system.time(resp6 <- algorithms(mydatam, myeffm, MLmodel, useOptim=TRUE,
                    optimFinal=2, scale=0.5, nbrNodes=1,
								optimSchedule=c(10,20,20,20,50)))

## Geyer with varied theta
print("resp7")
system.time(resp7 <- algorithms(mydatam, myeffm, MLmodel, useOptim=TRUE,
                    optimFinal=2, variedThetas=TRUE, scale=0.5,
                    finalIter=100, optimSchedule=c(10, 20, 20, 20, 50,50),
                    nbrNodes=1))


## response surface
print("resp8")
system.time(resp8 <- algorithms(mydatam, myeffm, MOMmodel, nIter=100, numiter=10,
                    finalIter=100, responseSurface=TRUE, scale=1,
                    optimWeight=0.75, nbrNodes=1))

## profileLikelihoods
system.time(RSiena:::profileLikelihoods(resp8,
                   MLmodel, mydatam, myeffm,
                   1,3, gridl=c(0.9,1.1),nIter=100, nbrNodes=1))
system.time(RSienaTest:::profileLikelihoods(resp8,
                   MLmodel, mydatam, myeffm,
                   1,3, gridl=c(0.9,1.1),nIter=100, nbrNodes=1))

##bayes
print("resp9")
system.time(resp9 <- sienaBayes(mydatam, myeffm, MLmodel,nwarm=10,nmain=10, nbrNodes=1, plot=TRUE))

