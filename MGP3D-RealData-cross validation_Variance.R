rm(list=ls()) 
library(Matrix)
library(nloptr)
library(minqa)
library(optimx)
library(rootSolve)
library(readxl)
########################################### System Input ##############################################
########################################### Specify Input #################################################
n=3 #Number of outputs (measurement)
p=3  #Number of Dimensions

SS316L <- read_excel("C:/Cloud Drives/Box Sync/Projects/Machine Learning/Data/dataset_Rready.xlsx",sheet = "SS316L", range = "B1:F28")

SS316L_Cu <- read_excel("C:/Cloud Drives/Box Sync/Projects/Machine Learning/Data/dataset_Rready.xlsx",sheet = "316-Cu-75%", range = "B1:F28")
SS316L_Cu=SS316L_Cu[complete.cases(SS316L_Cu[,4]),]

Cu <- read_excel("C:/Cloud Drives/Box Sync/Projects/Machine Learning/Data/dataset_Rready.xlsx",sheet = "Cu", range = "B1:F28")

data=rbind(cbind(id=1,SS316L),cbind(id=2,Cu),cbind(id=3,SS316L_Cu))
sdata1=scale(data[,-1],center=T,scale=T)
sdata=as.data.frame(cbind(id=data$id,sdata1))

iter=0 #Iteration
ae=c() #Absolute Error

data_iter=data.frame()

while(iter<5){
  iter=iter+1
  
  testy=as.vector(sdata$Density[sdata$id==3])
  
  cons=c(which(testy==max(testy)),which(testy==min(testy))) #Pick min and max from training data
  s_temp=sample((1:length(testy))[! 1:length(testy) %in%  cons],3,replace=F) #Select the rest randomly
  s_index=c(cons,s_temp)
  testy=testy[s_index]
  
  test=as.matrix(sdata[sdata$id==3,2:4])
  tests=test[s_index,]
  
  trains=list()
  trains[[1]]=as.matrix(sdata[sdata$id==1,2:4])
  trains[[2]]=as.matrix(sdata[sdata$id==2,2:4])
  
  trainy=c(as.vector(sdata$Density[sdata$id==1]),sdata$Density[sdata$id==2])
  m1=nrow(tests);
  ############################################# Given/ Matrix Shape ######################################## 
  index=function(n,len,m)
  {
    p1=c();p2=c();p3=c();p4=c();p5=c();p6=c()
    pp=sum(len)
    for(j in 1:(n-1))
    {
      i1=1 + sum(len[0:(j-1)])
      for(i in i1:(i1+len[j]-1))
      {
        p1=c(p1,i1:i)
        p2=c(p2,rep(i,length(i1:i)))
      }
    }
    p3=rep(1:pp,m)
    for(i in 1:m)
    {
      p4=c(p4,rep(pp+i,pp))
    }
    i2=pp+1
    for(i in i2:(i2+m-1))
    {
      p5=c(p5,i2:i)
      p6=c(p6,rep(i,length(i2:i)))
    }
    
    return(list(pfi=c(p1,p3,p5),pfj=c(p2,p4,p6)))
  }
  pf=index(n,lengths(trains)/p,m1)
  pfi=pf$pfi;pfj=pf$pfj
  ###################################### Covariance Functions #############################################
  cyii=function(a,b,L)
  {
    p=ncol(a)
    dd=lapply(1:p, function(i){outer(a[,i],b[,i],`-`)})
    II=lapply(1:p, function(i){outer(a[,i],b[,i],`==`)})
    I=Reduce("*",II)
    d=lapply(1:p, function(i){dd[[i]][upper.tri(dd[[i]],diag=T)]})
    I=I[upper.tri(I,diag=T)]
    L[1]^2/sqrt( Reduce("*",lapply(1:p,function(i){1/L[i+1]^2})))  *exp(-0.25*(Reduce("+",lapply(1:p,function(i){d[[i]]^2*L[i+1]^2})))) + I*L[p+2]^2
  }
  cyip=function(a,b,L)
  {
    p=ncol(a)
    d=lapply(1:p, function(i){outer(a[,i],b[,i],`-`)})
    2*L[1]*L[p+2]/sqrt((1/L[2]^2+1/L[6]^2)*(1/L[3]^2+1/L[7]^2)*(1/L[4]^2+1/L[8]^2))*exp(-0.5*(d[[1]]^2*(L[2]^2*L[6]^2)/(L[2]^2+L[6]^2)+d[[2]]^2*(L[3]^2*L[7]^2)/(L[3]^2+L[7]^2)+d[[3]]^2*(L[4]^2*L[8]^2)/(L[4]^2+L[8]^2)))
  }
  ###################################### Covariance Matrix #############################################
  y=c(trainy,testy)
  dd=lapply(1:p, function(i){outer(tests[,i],tests[,i],`-`)})
  II=lapply(1:p, function(i){outer(tests[,i],tests[,i],`==`)})
  d=lapply(1:p, function(i){dd[[i]][upper.tri(dd[[i]],diag=T)]})
  P=Reduce("*",II)
  P=P[upper.tri(P,diag=T)]
  leny=length(y)
  ######################################################################################################
  C=function(strain,H)
  {
    zii=list();zip=list();zpp=c()
    zii = lapply(1:(n-1), function(i){cyii(strain[[i]],strain[[i]],H[c(4*i-3,4*i-2,4*i-1,4*i,8*n-3)])})
    zip = lapply(1:(n-1), function(i){cyip(strain[[i]],tests,H[c(4*i-3,4*i-2,4*i-1,4*i,4*n+4*i-3,4*n+4*i-2,4*n+4*i-1,4*n+4*i)])})
    K=H[(4*n-3):(8*n-3)]
    zpp=Reduce("+",lapply(1:n, function(i){K[4*i-3]^2/sqrt((1/K[4*i-2]^2)*(1/K[4*i-1]^2)*(1/K[4*i]^2))*exp(-0.25*(Reduce("+",lapply(1:p,function(j){d[[j]]^2*K[4*i-(p-j)]^2}))))})) + P*K[length(K)]^2
    b1=unlist(zii);b2=as.vector(do.call("rbind",zip));
    return(sparseMatrix(i=pfi,j=pfj,x=c(b1,b2,zpp),symmetric=T))
  }
  ##################################### likelihood #####################################################
  logL=function(H,fn)
  {
    B=C(trains,H)
    deter=det(B)
    if(deter>0) {a=0.5*(log(deter)+t(y)%*%solve(B,y)+log(2*pi)*leny)
    } else {
      ch=chol(B)
      logdeter=2*(sum(log(diag(ch))))
      a=0.5*(logdeter+t(y)%*%solve(B,y)+log(2*pi)*leny)
    }
    return(as.numeric(a))
  }
  logL_grad=function(H,fn)
  {
    return(nl.grad(H,fn))
  }
  ################################## Optimize #############################################
  x0=c(rep(1,8*n-4),0.05)
  opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 10000,print_level=3) #print_level=3
  #t=proc.time()
  one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,opts= opts,fn= logL ), error = function(e) e)
  # if(any(class(one) == "error")==T){cat("nloptr Problem");iit=iit-1;next}
  #proc.time()-t
  H1=one$solution
  H0=H1
  ############################################ Results ######################################################
  zs=as.vector(sdata$Density[sdata$id==3])
  zs=zs[-(s_index)]
  xstar=test[-(s_index),]
  
  ### Prediction Variance and Mean
  zip_pred=list()
  zip_pred =lapply(1:(n-1), function(i){cyip(trains[[i]],xstar,H0[c(4*i-3,4*i-2,4*i-1,4*i,4*n+4*i-3,4*n+4*i-2,4*n+4*i-1,4*n+4*i)])})
  
  d1=lapply(1:p, function(i){outer(xstar[,i],tests[,i],`-`)})
  K1=H0[(4*n-3):(8*n-3)]
  zip_pred[[n]]=t(Reduce("+",lapply(1:n, function(i){K1[4*i-3]^2/sqrt((1/K1[4*i-2]^2)*(1/K1[4*i-1]^2)*(1/K1[4*i]^2))*exp(-0.25*(Reduce("+",lapply(1:p,function(j){d1[[j]]^2*K1[4*i-(p-j)]^2}))))})))
  
  Pk=t(do.call("rbind",zip_pred))
  
  d2=lapply(1:p, function(i){outer(xstar[,i],xstar[,i],`-`)})
  II=lapply(1:p, function(i){outer(xstar[,i],xstar[,i],`==`)})
  P2=Reduce("*",II)
  
  sk=Reduce("+",lapply(1:n, function(i){K1[4*i-3]^2/sqrt((1/K1[4*i-2]^2)*(1/K1[4*i-1]^2)*(1/K1[4*i]^2))*exp(-0.25*(Reduce("+",lapply(1:p,function(j){d2[[j]]^2*K1[4*i-(p-j)]^2}))))})) + P2*K1[length(K1)]^2
  
  covM=C(trains,H0)
  salman=solve(covM,y)
  ypred=Pk%*%salman
  yvar=diag(sk-Pk%*%solve(covM,t(Pk)))
  
  ae[iter]=sum(abs(zs-ypred))/length(zs)*100
  
  ypred2=ypred * attr(sdata1, 'scaled:scale')[4] + attr(sdata1, 'scaled:center')[4]  #Unscaled prediction
  zs2=zs * attr(sdata1, 'scaled:scale')[4] + attr(sdata1, 'scaled:center')[4]  #Unscaled original data

  yvar2=yvar * attr(sdata1, 'scaled:scale')[4] + attr(sdata1, 'scaled:center')[4]  #Unscaled prediction
  se=sqrt(yvar2)
  
  temp=sweep(xstar,MARGIN = 2,attr(sdata1, 'scaled:scale')[1:3],"*")
  predictionpoints=sweep(temp,MARGIN = 2,attr(sdata1, 'scaled:center')[1:3],"+")
  data_iter=rbind(data_iter,data.frame(Iteration=iter,Power=predictionpoints[,1],Velocity=predictionpoints[,2],Hatch=predictionpoints[,3],Density_Prediction=as.numeric(ypred2),Density_True=as.numeric(zs2),SE=se))
  }
# ############################################ Results - Write to file ######################################################
#write.csv(ae,file="C:/Users/rankouhi/Box/Projects/Machine Learning/Data/AbsoluteError.csv")
write.csv(data_iter,file="C:/Cloud Drives/Box Sync/Projects/Machine Learning/Data/Results/Varaince/Predictions.csv")


# plot(data_iter$Power,data_iter$Density_Prediction,ylim=c(1,10))
# par(new=T)
# plot(data_iter$Power,data_iter$Density_Prediction+se,ylim=c(1,10),col='red')
