NH.Monte.Carlo <- function(n, R, alpha, lambda){ 
  
  z=c()
  set.seed(369)
  alpha=alpha  #forma
  lambda=lambda #escala
  real_par=c(alpha,lambda)
  R=R
  n=n
  
  fdp <- function(param,x){
    alpha=param[1]
    lambda=param[2]
    alpha*lambda*(1+lambda*x)^(alpha-1) * exp(1-(1+lambda*x)^alpha)
  }
  
  loglik <- function(param,x){
    alpha=param[1]
    lambda=param[2]
    log(prod(fdp(param,x)))
  }
  
  score <- function(param,x){
    c( (n/alpha) +sum(log(1+lambda*x))- sum((1+lambda*x)^alpha * log(1+lambda*x)),
       (n/lambda)+(alpha-1) * sum(x*(1+lambda*x)^(-1)) - alpha*sum(x*(1+lambda*x)^(alpha - 1)) )
  }	
  
  dist_ks<-function(param,x)
  {
    alpha=param[1]
    lambda=param[2]
    F_acum<-(1-exp(1-(1+lambda*x)^alpha))
    Fn<-ecdf(x)
    max(abs(F_acum-knots(Fn)))
  }
  
  result=c(); resultS=c(); resultKS=c()
  i=0
  
  while(i < R){
    u<-runif(n)
    xi <- ((1-log(1-u))^(1/alpha)-1)/lambda 
    
    opt <- try(optim(real_par, loglik, method="BFGS", x=xi, control=list(fnscale=-1)),silent = T)
    optScore <- try(optim(real_par, loglik, score, method = "BFGS", x=xi, control=list(fnscale=-1)),silent=T)
    
    opt_ks<-optim(real_par,dist_ks,x=xi)
    resultKS<-rbind(resultKS,opt_ks$par)
    
    if((class(opt) != "try-error") && (class(optScore) != "try-error")){ i<-i+1
    result <-rbind(result,  opt$par)
    resultS<-rbind(resultS, optScore$par) 
    colnames(result) <- c("alpha","lambda") ; colnames(resultS)<- c("alphaScore","lambdaScore")
    if((100*i/R)%%25==0) print(c((100*i/R),"% ---------- R ----------"), quote = T)
    }
  }
  
  z$result <- result
  z$resultScore <- resultS
  z$resultKS <- resultKS
  
  #media
  medias<-colMeans(result)
  mediasS<-colMeans(resultS)
  mediasKS<-colMeans(resultKS)
  
  #vies
  vies<-medias-real_par
  viesS<-mediasS-real_par
  viesKS<-mediasKS-real_par
  
  #vies relativo
  vr<- (vies/real_par)*100
  vrS<- (viesS/real_par)*100
  vrKS<- (viesKS/real_par)*100
  
  #erro padrao
  sd <- apply(result,2,sd)  
  sdS <- apply(resultS,2,sd) 
  sdKS <- apply(resultKS,2,sd) 
  
  #EQM 
  eqm <- apply(result,2,var)+vies^2
  eqmS <- apply(resultS,2,var)+viesS^2 
  eqmKS <- apply(resultKS,2,var)+viesS^2 
  
  # final
  all<- cbind(real_par,medias, vies, vr, sd, eqm)
  allS<-cbind(real_par,mediasS,viesS,vrS,sdS,eqmS)
  allKS<-cbind(real_par,mediasKS,viesKS,vrKS,sdKS,eqmKS)
  
  z$medidas <- all
  z$medidasScore <- allS
  z$medidasKS <- allKS
  
  print(c(paste("Tamanho amostral =",n),paste("Replicas Monte Carlo =",R)),q=F)
  print(round(all,3))
  cat("\n")
  cat("Considerando escore analitico\n")
  print(round(allS,3))
  cat("\n")
  cat("Minima Distancia via KS\n")
  print(round(allKS,3))
  
  return(z)
}

alpha05_30 <-NH.Monte.Carlo(n=30, R=50000,alpha=0.5,lambda=1)
alpha05_100<-NH.Monte.Carlo(n=100,R=50000,alpha=0.5,lambda=1)
alpha05_300<-NH.Monte.Carlo(n=300,R=50000,alpha=0.5,lambda=1)

alpha1_30 <-NH.Monte.Carlo(n=30, R=50000,alpha=1,lambda=1)
alpha1_100<-NH.Monte.Carlo(n=100,R=50000,alpha=1,lambda=1)
alpha1_300<-NH.Monte.Carlo(n=300,R=50000,alpha=1,lambda=1)

alpha2_30 <-NH.Monte.Carlo(n=30, R=50000,alpha=2,lambda=1)
alpha2_100<-NH.Monte.Carlo(n=100,R=50000,alpha=2,lambda=1)
alpha2_300<-NH.Monte.Carlo(n=300,R=50000,alpha=2,lambda=1)

alpha3_30 <-NH.Monte.Carlo(n=30, R=50000,alpha=3,lambda=1)
alpha3_100<-NH.Monte.Carlo(n=100,R=50000,alpha=3,lambda=1)
alpha3_300<-NH.Monte.Carlo(n=300,R=50000,alpha=3,lambda=1)

# Plots	
fdpp <- function(x){
  alpha*lambda*(1+lambda*x)^(alpha-1) * exp(1-(1+lambda*x)^alpha)
}
integrate(fdpp,lower=0,upper=Inf) # OK. Confirmado

par(mar=c(2.5, 2.6, 0.5, 0.5)) # margens c(baixo,esq,cima,direia)
par(mgp=c(1.6, 0.5, 0))
alpha=0.5  #forma
lambda=1 #escala #fixa
curve(fdpp,from=fromx, to=tox, add = FALSE, lty=1, type = "l", xlim = c(0,2),
      ylab = expression("Densidade NH"), ylim =c(yliminf,1.8), col =1, lwd = 1.3)
alpha=1
curve(fdpp,from=fromx, to=tox, add = T, lty=2, type = "l", lwd = 1.3,col=1)
alpha=2
curve(fdpp,from=fromx, to=tox, add = T, lty=3, type = "l", lwd = 1.3,col=1)
alpha=3
curve(fdpp,from=fromx, to=tox, add = T, lty=5, type = "l", lwd = 1.3,col=1)
legend(1.6,1.8, legend = c(expression(paste(alpha,"=0.5"),paste(alpha,"=1"),paste(alpha,"=2"),paste(alpha,"=3"))),
       lty = c(1,2,3,5) , bty = "n", y.intersp = 1.3)

# HRF hazard
F_acum<-function(x) (1-exp(1-(1+lambda*x)^alpha))
surv<-function(x) 1-F_acum(x)
hrf <- function(x) fdpp(x)/surv(x)

alpha=0.5  #forma
lambda=1 #escala #fixa
curve(hrf,from=fromx, to=tox, add = FALSE, lty=1, type = "l", 
      ylab = expression("Taxa de risco (hrf)"), ylim =c(yliminf,15), col =1, lwd = 1.3)
alpha=1
curve(hrf,from=fromx, to=tox, add = T, lty=2, type = "l", lwd = 1.3,col=1)
alpha=2
curve(hrf,from=fromx, to=tox, add = T, lty=3, type = "l", lwd = 1.3,col=1)
alpha=3
curve(hrf,from=fromx, to=tox, add = T, lty=5, type = "l", lwd = 1.3,col=1)
legend(2,15.5, legend = c(expression(paste(alpha,"=0.5"),paste(alpha,"=1"),paste(alpha,"=2"),paste(alpha,"=3"))),
       lty = c(1,2,3,5) , bty = "n", y.intersp = 1.3)