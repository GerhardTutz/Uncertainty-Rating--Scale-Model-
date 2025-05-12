

### Example attidudes: role in society

## data from Maurer und Tutz (2024): An Ordinal Item Response Model for Understanding Attitudes, SSRN


#readRDS("C:/Users/tutz/LRZ Sync+Share/TuRegAdjCat/R/socgfsav.rds")

socgf<-readRDS("C:/Users/tutz/LRZ Sync+Share/TuRegAdjCat/R/socgfsav.rds")
summary(socgf)
dim(socgf)



##### v74  as response: family life suffers when woman has full-time job 


### CUB models

library("CUB")
CUB  <-GEM(Formula(v74~age+male+v261_ppp+townsize|age+male+v261_ppp+townsize|0), 
           family="cub",data =socgf)
summary(CUB)


CUB0  <-GEM(Formula(v74~0|age+male+v261_ppp+townsize|0), 
            family="cub",data =socgf)
summary(CUB0)

### own software 
source("ProgramsAdjacent")

pred <- as.matrix(socgf[,c(34,46,39,45)]) #age,male.income, townsize
#pred <- as.matrix(socgf[,c(34,46)]) #age,male.income, townsize
resp<-(as.matrix(as.numeric(socgf[,8]))) ### 74


### without dispersion

k<-4
fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = TRUE,maxit=500,start=0)
fit
fit$location
round(fit$location,digits=3)

#with dispersion

disp<-pred
fitdispder<-fitadj(resp,k,pred,disp=disp,der='der',hessian = TRUE,maxit=500,start=c(fit$parameter,rep(.1,dim(disp)[2])))  ## with der
fitdispder

round(fitdispder$location,digits=3)
round(fitdispder$parunc,digits=3)



