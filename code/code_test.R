### 数据预处理####


data <- read.csv("FirstData.csv")

data$sex <- as.numeric(data$sex==1)
data=data[order(data$NO,decreasing = F),]


#########################################################################  height 

data <- data[c("NO","height","time","sex",
               "birthweight","fheight","mheight")] #for height
#防止complete.case多删除数据！


#数据特殊处理，删除缺失值，删除重复测试点不到2次的人
data <- data[complete.cases(data),]
for(i in 1:431){
  if (length(data$NO[data$NO==i])<2){data=data[data$NO!=i,]}
}
n=0
for(i in 1:431){if (length(data$NO[data$NO==i])!=0){n=n+1;data$iden[data$NO==i]=n}}
ID <- data$iden
N1=length(ID)
R <- c(1:N1)
T <- data$time/730
############################自变量转化
sex <- data$sex
birthweight <- data$birthweight/1000 #把单位克换成千克。保持统一。
fheight <- data$fheight
mheight <- data$mheight
X_fix <- cbind(birthweight,fheight,mheight)
nfix <- ncol(X_fix)   #不随时间变化的变量个数

Y <- data$height          #因变量转化





################### 需要的函数 ###########

#knots function
default.knots=function(x,num.knots)
{  if (missing(num.knots)) num.knots=max(5,min(floor(length(unique(x))/4),40))
return(quantile(unique(x),seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))]))}

#variance smoother
newtonr=function(lambda2, eta, W2, R, tol=10^(-2), maxiter=25){
  dif=1; k=0
  eta1 <- eta                                               #GAILE
  while(dif>tol & k<maxiter){
    k=k+1
    sigma=exp(W2%*%eta)
    nR=R/sqrt(sigma)
    H1=0
    for(i in 1:n){ H1=H1+t(W2[ID==i,]*nR[ID==i])%*%(nR[ID==i]*W2[ID==i,])}
    u=drop(apply(W2, 2, sum)-t(nR^2)%*%W2)+drop(2*lambda2*D2%*%eta)
    H=H1+2*lambda2*D2
    eta1=eta-solve(H)%*%u
    lambda2=1/(mean((diag(solve(H)))[(p2+1):(p2+K2)])+mean((eta1[(p2+1):(p2+K2)])^2))
    dif=sum((eta1-eta)^2)
    eta=eta1}
  print(lambda2)
  return(list(eta=eta1,lambda=lambda2))}        #GAILE

#update smoothing parameter for variance smoother
eta_f=function(lambda2, eta, W2, R, tol=10^(-2), maxiter=2){
  k=0; dif=1
  while(dif>tol && k < maxiter){
    k=k+1
    temp=newtonr(lambda2,eta, W2, R)
    ERROR=temp$ERROR                                        #GAILE
    nlambda2=temp$lambda
    dif=abs(nlambda2-lambda2)
    lambda2=nlambda2
  }
  return(list(eta=temp$eta, lambda2=lambda2))}  #GAILE


#EM algorithm
EM_f=function(beta.hat, sigma, lambda1, lambda3, lambdanu, D,tol=10^(-2), maxiter=25){
  Dnu=solve(solve(D)+lambdanu*D3)
  dif=1; k=0
  while(dif>tol & k<maxiter){
    k=k+1
    a=d=0
    for(i in 1:n){
      varm=diag(sigma[ID==i])
      invcovm=solve(varm+ Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]))
      a=a+t(M[ID==i,])%*%invcovm%*%M[ID==i,]
      d=d+t(M[ID==i,])%*%invcovm%*%Y[ID==i]}
    nbeta.hat=solve(a+ c(rep(lambda1, nfix+p+K1), rep(lambda3,p+K1))*D12)%*%d   #gaile  #gaile2
    
    b=matrix(rep(0,n*(pnu+Knu)),nrow=n)
    R=Y-M%*%nbeta.hat                          #####gaile           
    for(i in 1:n){
      varm=diag(sigma[ID==i])
      invcovm=solve(varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]))
      b[i,]=Dnu%*%t(Znu[ID==i,])%*%invcovm%*%(R[ID==i])}
    
    D0=0
    for(i in 1:n){
      varm=diag(sigma[ID==i])
      invcovm=solve(varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]))
      nQ=invcovm-invcovm%*%M[ID==i,]%*%solve(a+lambda1*D12)%*%t(M[ID==i,])%*%invcovm
      D0=D0+b[i,]%*%t(b[i,])+Dnu-Dnu%*%t(Znu[ID==i,])%*%nQ%*%Znu[ID==i,]%*%Dnu}
    
    dif=mean((beta.hat-nbeta.hat)^2)
    beta.hat=nbeta.hat}
  
  return(list(D=D0/n, nbeta.hat=nbeta.hat, b.hat=b))}

#choice of smoothing parameter  for mean function
lambda1f=function(lambda, lambda3, lambdanu, sigma, beta, maxiter=3){
  r=Y-M1%*%c(beta[1:(p+nfix)], beta[(p+K1+1+nfix):(2*p+K1+nfix)])     ##gaile  #gaile2 除了第一项，全家加1
  n.row=0
  V1=matrix(rep(0,N1^2),nrow=N1)
  Dnu=solve(solve(D)+lambdanu*D3)
  for(i in 1:n){
    temp1=length(Y[ID==i])
    varm=diag(sigma[ID==i])
    V1[(n.row+1):(n.row+temp1),(n.row+1):(n.row+temp1)]=varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,])
    n.row=n.row+temp1}
  
  dif=1; k=0
  while(dif>.01 & k<maxiter){
    V=V1+Z1%*%t(Z1)*c(lambda)+GZ1%*%t(GZ1)/lambda3   ##maybe wrong
    inV=solve(V)
    dV=Z1%*%t(Z1)  
    dinV=-inV%*%dV%*%inV
    tempm2=solve(t(M1)%*%inV%*%M1)
    u=sum(diag(inV%*%dV))+sum(diag(tempm2%*%t(M1)%*%dinV%*%M1))+t(r)%*%dinV%*%r
    H=sum(diag(dinV%*%dV))-2*t(r)%*%dinV%*%dV%*%inV%*%r-sum(diag(tempm2%*%t(M1)%*%dinV%*%M1%*%tempm2%*%t(M1)%*%dinV%*%M1))-2*sum(diag(tempm2%*%t(M1)%*%dinV%*%dV%*%inV%*%M1))
    nlambda=lambda-u/H;
    dif=abs(nlambda-lambda)
    lambda=nlambda
    k=k+1}
  
  return(c(1/nlambda))}

#choice of smoothing parameter  for varying coefficient
lambda3f=function(nlambda=25, lambda1, lambdanu, sigma, beta){
  REML=rep(0, nlambda) 
  r=Y-M1%*%c(beta[1:(p+nfix)], beta[(p+K1+1+nfix):(2*p+K1+nfix)])     ##gaile  #gaile2 除了第一项，全家加1
  V1=matrix(rep(0,N1^2),nrow=N1) 
  Dnu=solve(solve(D)+lambdanu*D3)
  for(j in 1:nlambda){ lambda3=1.5^(-1/2*nlambda+j) 
  n.row=0 
  for(i in 1:n){temp1=length(Y[ID==i]) 
  varm=diag(sigma[ID==i])
  V1[(n.row+1):(n.row+temp1),(n.row+1):(n.row+temp1)]=varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]) 
  n.row=n.row+temp1}
  V=V1+Z1%*%t(Z1)/lambda1++GZ1%*%t(GZ1)/lambda3   ##maybe wrong
  inV=solve(V) 
  a=det(V);  if (a<10^(-300)) loga=-300*log(10) else loga=log(a)
  REML[j]=loga+t(r)%*%inV%*%r+log(det(t(M1)%*%inV%*%M1))
  if(is.na(REML[j]) | REML[j]>10^5) REML[j]=10^5 }
  vec=cbind(REML, (1:nlambda))[order(REML),]
  return(1.5^(-1/2*nlambda+vec[1,2]))}


#choice of smoothing parameter for subject specific curves
lambdanuf=function(nlambda=25, sigma, lambda1){
  REML=rep(0,nlambda)
  for(j in 1: nlambda){
    lambdanu=1.5^(-1/2*nlambda+j) 
    r=Y-M%*%nbeta.hat
    Dnu=solve(solve(D)+lambdanu*D3) 
    
    temp1=0
    for(i in 1:n){
      varm=diag(sigma[ID==i])
      covm=(varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]))  
      invcovm=solve(covm)
      REML[j]=REML[j]+log(det(covm))+t(r[ID==i])%*%invcovm%*%r[ID==i]
      temp1=temp1+t(M[ID==i,])%*%invcovm%*%M[ID==i,]}
    
    REML[j]=REML[j]+log(det(temp1))
    if(is.na(REML[j])) REML[j]=10^5 }
  vec=cbind(REML, (1:nlambda))[order(REML),]
  return(1.5^(-1/2*nlambda+vec[1,2]))}

############需要的函数 ending ###########



###开始分析###

K1=40;p=4 
Knu=40;pnu=2
K2=40;p2=4

#mean function

#T是时间，K1是mu的knots的个数。Knots1是切割好的时间点。
#N1是subject*meansurement.每个个体的knots都一样.

knots1=abc
knots1_esti=knots1

X1=cbind(rep(1,N1), T, T^2,T^3) 
Z1=outer(T, knots1,"-") 
Z1=(Z1*(Z1>0))^3
W1=cbind(X1,Z1) 
D1=diag(c(rep(0,p),rep(1,K1)))  #就是penalty matrix。

#varying coefficient

#p=the order of basis function? 和pnu指的是什么？#############改了 

GX1=X1*sex
GZ1=Z1*sex
G=cbind(X1,Z1)*sex


D12=diag(c(rep(0,nfix),rep(0,p),rep(1,K1),rep(0,p),rep(1,K1)))  #gaile2
M1=cbind(X_fix,X1,GX1)                                             #gaile2
M=cbind(X_fix,W1,G) #rewrite后的Xi矩阵。####################改完     #gaile2

#random subject specific curves

knotsnu=default.knots(T, Knu)

Xnu=cbind(rep(1,N1), T) 
Znu=outer(T, knotsnu,"-") 
Znu=(Znu*(Znu>0))
Znu=cbind(Xnu,Znu) 
D3=diag(c(rep(0,pnu),rep(1,Knu)))

#variance function
#D是penalty matrix。p是order
knots2=default.knots(T, K2) 

X2=cbind(rep(1,N1), T,T^2,T^3) 
Z2=outer(T, knots2,"-") 
Z2=(Z2*(Z2>0))^3
W2=cbind(X2,Z2)
D2=diag(c(rep(0,p2),rep(1,K2)))


library(nlme)
group=rep(1,N1)
fit=lme(Y~-1+X_fix+X1+GX1, random=list(group=pdIdent(~-1+Z1+GZ1)))  ######gaile   gaile2
beta.hat=c(fit$coef$fixed[1:(p+nfix)],unlist(fit$coef$random)[1:K1],
           fit$coef$fixed[(p+1+nfix):(2*p+nfix)],unlist(fit$coef$random)[(1+K1):(2*K1)]) #gaile 2: fix effect +1
#分别是截距项固定和trt项固定。
lambda1=as.numeric(VarCorr(fit)[(2*K1+1),1])/as.numeric(VarCorr(fit)[K1,1])  #maybe wrong
lambdanu=1;  lambda3=lambda1
D=diag(rep(1,(pnu+Knu)))
sigma.hat=rep(as.numeric(VarCorr(fit)[(2*K1+1),1]), N1)   #######gaiwan

temp=EM_f(beta.hat, sigma.hat, lambda1, lambda3, lambdanu, D)
D=temp$D
nbeta.hat=temp$nbeta.hat
b.hat=temp$b.hat
for(i in 1:n){R[ID==i]=Y[ID==i]-M[ID==i,]%*%beta.hat-Znu[ID==i,]%*%b.hat[i,]}

fit2=lme(log(R^2)~-1+X2, random=list(group=pdIdent(~-1+Z2)))
eta=c(fit2$coef$fixed, unlist(fit2$coef$random))
lambda2=as.numeric(VarCorr(fit2)[(K2+1),1])/as.numeric(VarCorr(fit2)[K2,1])
variance=eta_f(lambda2, eta, W2,R)
eta1=variance$eta
lambda2=variance$lambda2
sigma.hat=exp(W2%*%eta1) 

count=0;   maxiter=5;  dif=1
while(dif>.01 && count<maxiter){
  lambda1=lambda1f(1/lambda1, lambda3, lambdanu, sigma.hat, beta.hat,maxiter=2) 
  lambda3=lambda3f(nlambda=15, lambda1, lambdanu, sigma.hat, beta.hat)
  lambdanu=lambdanuf(nlambda=15, sigma.hat, lambda1) 
  
  temp=EM_f(beta.hat, sigma.hat, lambda1, lambda3, lambdanu, D)
  nD=temp$D
  nbeta.hat=temp$nbeta.hat
  b.hat=temp$b.hat
  
  for(i in 1:n){R[ID==i]=Y[ID==i]-M[ID==i,]%*%nbeta.hat-Znu[ID==i,]%*%b.hat[i,]}
  variance=eta_f(lambda2, eta1, W2,R)
  eta1=variance$eta
  lambda2=variance$lambda2
  sigma.hat=exp(W2%*%eta1)
  dif=mean(diag(D-nD)^2)
  D=nD
  count=count+1
}

R.hat <- R
beta.hat0=nbeta.hat



beta_l <- length(beta.hat0)-nfix  #gaile2
mu_t <- function(t)(c(1,t,t^2,t^3,(t-knots1_esti)^3*(t>knots1_esti))%*%beta.hat0[(1+nfix):(beta_l/2+nfix)])#gaile2
x <- seq(0,1,length.out = 100)
y <- sapply(x,mu_t)



beta_t <- function(t)(c(1,t,t^2,t^3,(t-knots1)^3*(t>knots1))%*%beta.hat0[(beta_l/2+1+nfix):(2*beta_l/2+nfix)])#gaile2
yy <- sapply(x,beta_t)


beta.hat0[1:nfix]





#postscript(file='height3.eps',horizontal=T,width=12,height=6)

par(mfrow=c(1,2)) 
plot(2*x,y,xlab = "age",col= "red",ylab=expression(mu(t)),type="l",lwd=2,ylim=c(0,65),
     main = "mean function for height") #之前将T除以730，即两年，故这里用2*x


plot(2*x,yy,xlab = "age",ylab=expression(beta(t)),type="l",lwd=2,ylim=c(0,2),
     main = "beta(t) of sex for height",col="blue")

abline(h=0,lty=4,col=6,lwd=2)
par(mfrow=c(1,1))

#dev.off()