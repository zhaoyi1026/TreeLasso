###############################
# Hierarchical Tree Data in Regularized Regression
# Zhao et al.

# Example code
###############################

###############################

library("genlasso")        # generalized lasso pacakge
library("glmnet")
library("elasticnet")      # elastic net package
library("glmgraph")        # glmgraph package (Li and Li, 2008), graph-constrained regularization (GCR)
library("MLGL")            # overlapped grop lasso package

rm(list=ls())
###############################

###############################
# model setting parameters

p<-127        # # of predictors
L<-7          # # of levels of a binary tree design

root.idx<-c(1)                    # root node index
rest.idx<-(1:p)[-root.idx]        # rest node

level.idx<-vector("list",length=L)            # node index in each level
names(level.idx)<-paste0("Level",1:L)
for(l in 1:L)
{
  if(l==1)
  {
    level.idx[[l]]<-1
  }else
  {
    level.idx[[l]]<-((2^(l-1)-1)+1):(2^l-1)
  }
}

A.base<-1                                    # weights in the adjacency matrix
A<-matrix(0,p,p)                             # adjacency matrix
colnames(A)=rownames(A)<-paste0("node",1:p)
for(l in 1:(L-1))
{
  for(i in 1:length(level.idx[[l]]))
  {
    pidx.tmp<-level.idx[[l]][i]
    cidx.tmp<-level.idx[[l+1]][((i-1)*2+1):(i*2)]
    A[pidx.tmp,cidx.tmp]<-A.base
  }
}
D<-solve(diag(1,p)-A)                        # influence matrix

beta.base<-1                                 # direct effect value (beta)
# Simulation (1) Model (a)-(c)
beta.mat1<-matrix(0,p,3)                     
rownames(beta.mat1)<-rownames(A)
colnames(beta.mat1)<-paste0("model",1:ncol(beta.mat1))
s<-0
for(l in c(1,3,7))
{
  s<-s+1
  beta.mat1[level.idx[[l]][1],s]<-beta.base
}
# Simulation (1) Model (d)
beta.mat2<-matrix(0,p,1)                   
rownames(beta.mat2)<-rownames(A)
colnames(beta.mat2)<-"model4"
l<-6
beta.mat2[level.idx[[l]][1],1]<--beta.base
beta.mat2[level.idx[[L]][1],1]<-beta.base
# Simulation (1) Model (e)
beta.mat3<-matrix(0,p,1)                   
rownames(beta.mat3)<-rownames(A)
colnames(beta.mat3)<-"model5"
for(l in seq(2,L-1,by=2))
{
  beta.mat3[level.idx[[l]][1],1]<-beta.base*(rep(c(-1,1),ceiling(length(seq(2,L-1,by=2))/2))[which(l==seq(2,L-1,by=2))])
}
beta.mat3[level.idx[[L]][1],1]<-beta.base
# beta coefficient of 5 models
beta.mat<-cbind(beta.mat1,beta.mat2,beta.mat3)

gamma.mat<-matrix(NA,nrow(beta.mat),ncol(beta.mat))      # total effect
rownames(gamma.mat)<-rownames(beta.mat)
colnames(gamma.mat)<-colnames(beta.mat)
for(ss in 1:ncol(beta.mat))
{
  gamma.mat[,ss]<-D%*%beta.mat[,ss]
}

eps.sd<-1
Y.sd<-1
###############################

###############################
# generate simulation data

n<-100          # sample size

dir.out<-paste0(getwd(),"/n",n)
if(file.exists(dir.out)==FALSE)
{
  dir.create(dir.out)
}

runonce<-function(l)
{
  #########################################
  # generate data
  data.file<-paste0(dir.out,"/RUN_",l,"_Data.RData")
  if(file.exists(data.file)==FALSE)
  {
    set.seed(100+l)
    
    epsilon<-matrix(rnorm(n*p,mean=0,sd=eps.sd),n,p)
    X<-epsilon%*%D
    colnames(X)<-colnames(A)
    
    Y<-matrix(NA,n,ncol(beta.mat))
    colnames(Y)<-colnames(beta.mat)
    for(ss in 1:ncol(beta.mat))
    {
      set.seed(500+l+ss)
      Y[,ss]<-X%*%beta.mat[,ss]+rnorm(n,mean=0,sd=Y.sd)
    }
    
    save(list=c("X","epsilon","Y"),file=data.file)
  }else
  {
    load(data.file)
  }
  #########################################
  
  #########################################
  run.file<-paste0(dir.out,"/RUN_",l,"_",method,".RData")
  
  if(file.exists(run.file)==FALSE)
  {
    if(method=="R1")
    {
      # R1 only
      re<-vector("list",length=ncol(beta.mat))
      names(re)<-colnames(beta.mat)
      for(ss in 1:ncol(beta.mat))
      {
        try(re[[ss]]<-genlasso(Y[,ss],X,D))
      }
    }
    if(method=="R2")
    {
      # R2 only
      re<-vector("list",length=ncol(beta.mat))
      names(re)<-colnames(beta.mat)
      for(ss in 1:ncol(beta.mat))
      {
        try(re[[ss]]<-genlasso(Y[,ss],X,diag(rep(1,p))))
      }
    }
    if(length(grep("R1R2-",method))>0)
    {
      # R1 and R2
      alpha.idx<-as.numeric(sub("R1R2-","",method))
      re<-vector("list",length=ncol(beta.mat))
      names(re)<-colnames(beta.mat)
      D.t<-rbind(D,alpha[alpha.idx]*diag(rep(1,p)))
      for(ss in 1:ncol(beta.mat))
      {
        try(re[[ss]]<-genlasso(Y[,ss],X,D.t))
      }
    }
    if(length(grep("EN-",method))>0)
    {
      # elastic net
      alpha.EN.idx<-as.numeric(sub("EN-","",method))
      re<-vector("list",length=ncol(beta.mat))
      names(re)<-colnames(beta.mat)
      for(ss in 1:ncol(beta.mat))
      {
        try(re[[ss]]<-enet(X,Y[,ss],lambda=alpha.EN[alpha.EN.idx],intercept=FALSE))
        if(is.null(re[[ss]])==FALSE)
        {
          re[[ss]]$beta<-t(re[[ss]]$beta.pure)
          re[[ss]]$lambda<-re[[ss]]$penalty
          re[[ss]]$y<-Y[,ss]
          re[[ss]]$fit<-X%*%re[[ss]]$beta
        }
      }
    }
    if(method=="GCR")
    {
      # calculate the Laplacian matrix
      A.sym<-A+t(A)          # symmetric adjacency matrix
      dg.vec<-apply(A.sym,1,sum)     # weighted degree
      L.mat<-diag(rep(1,ncol(A.sym)))-diag(1/sqrt(dg.vec))%*%A.sym%*%diag(1/sqrt(dg.vec))     # normalized Laplacian matrix
      
      re<-vector("list",length=ncol(beta.mat))
      names(re)<-colnames(beta.mat)
      for(ss in 1:ncol(beta.mat))
      {
        try(re[[ss]]<-cv.glmgraph(X,Y[,ss],L=L.mat,type.measure="mse",standardize=TRUE,intercept=FALSE,trace=FALSE))
      }
    }
    if(method=="GL")
    {
      group.var<-NULL                        # variable index included
      group<-NULL                            # group indicators
      s<-0
      for(i in 1:nrow(A))
      {
        if(sum(A[i,])!=0)
        {
          s<-s+1
          itmp<-which(A[i,]!=0)
          group.var<-c(group.var,itmp)
          group<-c(group,rep(s,length(itmp)))
        }
      }
      
      re<-vector("list",length=ncol(beta.mat))
      names(re)<-colnames(beta.mat)
      for(ss in 1:ncol(beta.mat))
      {
        try(re[[ss]]<-overlapgglasso(X=X,y=Y[,ss],var=group.var,group=group,intercept=FALSE))
        
        re[[ss]]$y<-Y[,ss]
        re[[ss]]$fit<-matrix(NA,n,length(re[[ss]]$lambda))
        out.beta<-matrix(NA,p,length(re[[ss]]$lambda))
        re[[ss]]$df<-rep(NA,length(re[[ss]]$lambda))
        for(i in 1:length(re[[ss]]$lambda))
        {
          re[[ss]]$df[i]<-length(re[[ss]]$var[[i]])
          beta.est.tmp<-rep(0,ncol(X))
          otmp<-re[[ss]]$beta[[i]]
          if(length(otmp)>0)
          {
            beta.est.tmp[re[[ss]]$var[[i]]]<-otmp
          }
          
          out.beta[,i]<-beta.est.tmp
          re[[ss]]$fit[,i]<-X%*%beta.est.tmp
        }
        re[[ss]]$beta<-out.beta
      }
    }
    
    save(list=c("re"),file=run.file)
  }
  #########################################
}
###############################

###############################
# method parameters

# R2/R1 ratio
alpha<-c(0.1,0.5,0.9,1.5)

# elastic net ratio
alpha.EN<-c(0.01,0.05,0.1,0.5,1,10)
###############################

###############################
# Methods:
# EN: elastic net with different alpha values
# R1
# R1R2: with different alpha values
# R2: lasso
# GCR: graph-constrained regularization
# GL: group lasso
method0<-c(paste0("EN-",1:length(alpha.EN)),"R1",paste0("R1R2-",1:length(alpha)),"R2","GCR","GL")

for(ii in 1:length(method0))
{
  method<-method0[ii]
  
  for(l in 1:L)
  {
    runonce(l)

    print(paste0("Run ",l," done!"))
  }
  
  print(paste0("Method ",method," done!"))
}
###############################










