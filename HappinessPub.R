

##### Example happiness: adjacent categories models 

library("CUB")
data(relgoods)
summary(relgoods)


### preparation data

####
myvars <- c("Gender", "Safety", "BirthYear", "EducationDegree", "WalkAlone","RelNeighbours", "Environment","Happiness")
newdata <- relgoods[myvars]

newdata$age <- 2014 -newdata$BirthYear

summary(newdata)
newdata <- newdata[newdata$age>=18,]
newdata <- newdata[newdata$age<=80,]
summary(newdata)

newdata <- na.omit(newdata) 
names(newdata)<-c("Gender", "Safety","Birthyear","Education","Walk","Neighbours","Environment","Happiness","age")

relgood2 <- newdata
relgood2$happicat<-floor((relgood2$Happiness/112)*10)+1 
relgood2$happicat<-as.factor(relgood2$happicat) 
#hcat<-round((relgood2$Happiness/110)*9,digits=0) 
summary(relgood2)

plot(relgood2$Happiness)
x <- relgood2$Happiness
h<-hist(x, breaks=9, col="blue", xlab="",
        main="")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="red", lwd=2)

summary(relgood2$happicat)



#### fitting with vgam

fitvglmh <- vglm(happicat~age+Gender+Education+Walk+Neighbours,family=acat(parallel=TRUE),data=relgood2)
summary(fitvglmh)

#fitcumh <- vglm(happicat~age+Gender+Education+Walk+Neighbours,family=cumulative(parallel=TRUE),data=relgood2)
#fitcumh
#summary(fitcumh)

## fitting with own function

source("ProgramsAdjacent")

pred <- as.matrix(relgood2[,c(9,1,4,5,6)]) #age,gender,education,walk, neighbors
resp<-(as.matrix(as.numeric(relgood2[,10]))) ### happiness


##  model without dispersion

k<-10
fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = TRUE,maxit=500,start=0)
fit
fit$location
round(fit$location,digits=3)

### model with dispersion

disp<-pred
fitdispder<-fitadj(resp,k,pred,disp=disp,der='der',hessian = TRUE,maxit=500,start=c(fit$parameter,rep(.1,dim(disp)[2])))  ## with der
fitdispder

round(fitdispder$location,digits=3) 
round(fitdispder$parunc,digits=3)


#### CUB models

CUB  <-GEM(Formula(happicat~age+Gender+Education+Walk+Neighbours|age+Gender+Education+Walk+Neighbours|0), 
           family="cub",data =relgood2)
summary(CUB)

### without covariates in mixture
CUB0  <-GEM(Formula(happicat~0|age+Gender+Education+Walk+Neighbours|0), 
            family="cub",data =relgood2)
summary(CUB0)



