



####  Example  Nuclear Energy adjacent cat scaled

library("ordinal")
library("VGAM")
library("xtable")
library("ordinalgmifs")
library("ordDisp")



load("./GLES17angst.rda")
dat <- GLES
summary(GLES)
### 

#set.seed(1)
#sam<-sample(dim(GLES)[1], 800, replace = FALSE, prob = NULL)



#### use all observations
sam<-1:2036
GLESred<-GLES[sam,]
pred <- as.matrix(GLES[sam,c(7,8,9,11)]) #age,gender,eastwest unemployment
resp<-as.matrix(GLES[sam,6]) ### nuclear energy

#pred <- as.matrix(GLES[sam,7:11]) #all
#resp<-as.matrix(GLES[sam,2]) ### climate change
#resp<-as.matrix(GLES[sam,4]) ### globalization

k<-7

##### Nuclear energy

formula <-"NuclearEnergy ~Age+Gender+EastWest+Unemployment"
formuladisp<-"~Age+Gender+EastWest+Unemployment"

#### fit model

#### vglm fit adjacent categories model

GLESred<-GLES[sam,]
fitvglm <- vglm(formula,family=acat(parallel=TRUE),data=GLESred)
fitvglm

#### vglm fit cumulative model

fitvcum <- vglm(formula,family=cumulative(parallel=TRUE),data=GLESred)
fitvcum

### without dispersion own functions

k<-7
fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = TRUE,maxit=500,start=0)
fit
round(fit$location,digits=3) 

### with dispersion start from fit

disp<-as.matrix(pred)
fitdispder<-fitadj(resp,k,pred,disp=disp,der='der',hessian = TRUE,maxit=500,start=c(fit$parameter,rep(.1,dim(disp)[2])))  ## with der
fitdispder

round(fitdispder$location,digits=3) 
round(fitdispder$parunc,digits=3)



#############  CUB models

library("CUB")
install.packages("FastCUB")
library("FastCUB")
#library("ordDisp")

#######################GESIS
#citation("CUB")


CUB  <-GEM(Formula(NuclearEnergy~Age+Gender+EastWest+Unemployment|Age+Gender+EastWest+Unemployment|0), 
           family="cub",data =GLES)
summary(CUB)

### without covariates in mixture
CUB0  <-GEM(Formula(NuclearEnergy~0|Age+Gender+Unemployment+EastWest|0), 
            family="cub",data =GLES)
summary(CUB0)

