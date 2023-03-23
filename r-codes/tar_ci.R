###########################################################################
##
##  TAR_CI.R
##  To estimate and test a threshold bi-variate VECM
## 
##  written by:
##  
##  Bruce E. Hansen
##  Department of Economics
##  Social Science Building
##  University of Wisconsin
##  Madison, WI 53706-1393
##  behansen@wisc.edu
##  http://www.ssc.wisc.edu/~bhansen/
##  
##  and 
##  
##  Byeongseon Seo
##  Department of Economics
##  Soongsil University
##  Seoul, 156-743
##  Korea
##  seo@saint.soongsil.ac.kr
##  
##  
##  This R program estimates a bi-variate VECM, a threshold bi-variate VECM, and
##  tests for the presence of a threshold.  The methods are those described in
##  "Testing for Threshold Cointegration" by Bruce E. Hansen and Byeongseon Seo.
##  
##  The program is set up to replicate the empirical application to the 3-month 
##  and 6-month interest rates series.  For your own application, load your data 
##  into the matrix "dat", and change the controls listed below
##  
###########################################################################

# Controls #

k <- 1	      # Lags in VAR beyond EC 
gn <- 300	      # number of gridpoints for gamma 
bn <- 300		# number of gridpoints for beta      
trim <- .05		# trimming percentage for threshold 
boot <- 5000	# number of bootstrap replications 
		      # set equal to zero to not do testing  
coint <- 1		# set to 1 to estimate cointegrating vector
			# set to 0 to fix cointegrating vector at _cvalue 
cvalue <- 1		# cointegrating vector, if coint=0 
cov <- 1		# covariance matrix estimation method
		      # set to 1 for Eicker-White
		      # set to 0 for conventional homoskedastic estimator 
p_ests <- 1       # set to 1 to print estimates, else 0 
graph <- 1        # set to 1 to generate graph of nonlinear ECM, else 0 
graph_rotate <- 0	# set to 1 to generate rotated graphs 
		      # (useful for some print jobs, but ackward for screen viewing) 

# Load your own data into matrix "dat" #
# Here we load in the 3-month and 6-month T-Bill series #

dat <- read.table("zeroyld.dat")
dat <- dat[1:nrow(dat),(7:62)]
rs <- rbind(as.matrix(seq(0,18,1)),21,24,30,as.matrix(seq(36,(36+7*12),12)))
short <- 12
long <- 120
short_i <- which.max(rs==short)
long_i <- which.max(rs==long)
dat <- dat[,cbind(long_i,short_i)]   

#*************************************************************#
#   MAIN PROGRAM                                              #
#*************************************************************#

for (main in 1:1){
if (p_ests==1){
cat ("**********************************************************", "\n")
cat ("\n")
cat ("Number of Bootstrap Replications    ", boot, "\n")
cat ("Number of Gridpoints for threshold  ", gn, "\n")
if (coint==1){
cat ("Estimated Cointegrating Vector", "\n")
cat ("Number of Gridpoints for ci vector ", bn, "\n")
}else{
cat ("Cointegrating Vector fixed at      ", cvalue, "\n")
}
cat ("\n")
cat ("\n")
cat ("Long Rate  (month):     ", long, "\n")
cat ("Short Rate (month):     ", short, "\n")
cat ("Number of VAR lags:     ", k, "\n")
cat ("\n")
}

# Organize Data #
n <- nrow(dat)
y <- as.matrix(dat[(2+k):n,]-dat[(1+k):(n-1),])
t <- nrow(y)
xlag <- as.matrix(dat[(1+k):(n-1),])
x <- matrix(1,t,1)
for (j in 1:k)  x <- cbind(x,(dat[(2+k-j):(n-j),]-dat[(1+k-j):(n-1-j),]))
x <- as.matrix(x)

# Compute Linear Model #
xx <- solve(t(x)%*%x)
xxx <- x%*%xx
u <- y-x%*%xx%*%(t(x)%*%y)
if (coint==1){
  v <- xlag-x%*%xx%*%(t(x)%*%xlag)
  uu <- t(u)%*%u
  m <- solve(t(v)%*%v)%*%(t(v)%*%u)%*%solve(uu)%*%(t(u)%*%v)
  ev <- eigen(m)
  va <- ev$val
  ve <- ev$vec
  h <- ve[,which.max(va)]
  b0 <- -h[2]/h[1]
}else{
  b0 <- cvalue
}
w0 <- as.matrix(xlag%*%rbind(1,-b0))
z0 <- cbind(w0,x)
kk <- ncol(z0)
zz0 <- solve(t(z0)%*%z0)
zzz0 <- z0%*%zz0
beta0 <- t(zzz0)%*%y 
e <- y - z0%*%beta0 
sige <- t(e)%*%e/t 
nlike <- (t/2)*log(det(sige))
bic <- nlike+log10(t)*4*(1+k)
aic <- nlike+2*4*(1+k)
b_like <- function(b){
  z <- cbind((xlag%*%rbind(1,-b)),x) 
  yz <- y - z%*%qr.solve(z,y)
  sigma <- (t(yz)%*%yz)/t  
  like <- (t/2)*log(det(sigma))
  like
}
nlike1 <- b_like(b0+.001)
nlike2 <- b_like(b0-.001)
hp <- (nlike1+nlike2-2*nlike)/(.001^2)
seb <- 1/sqrt(hp)
k_product <- function(ma,mb){
    mat <- matrix(0,nrow(ma)*nrow(mb),ncol(ma)*ncol(mb))
    for (i in 1:nrow(ma)){
        for (j in 1:ncol(ma)){
            mat[(1:nrow(mb))+(i-1)*nrow(mb),(1:ncol(mb))+(j-1)*ncol(mb)] <- mb*ma[i,j]
        }
    }
    mat
}
if (cov ==1){
  xe <- cbind((z0*(e[,1]%*%matrix(1,1,ncol(z0)))),(z0*(e[,2]%*%matrix(1,1,ncol(z0)))))
  m0 <- k_product(diag(2),zz0)
  v <- m0%*%t(xe)%*%xe%*%m0
}else{
  v <- k_product(sige,zz0)
}
se <- as.matrix(sqrt(diag(v)))

if (p_ests==1){
cat ("Linear VECM Estimates", "\n")
cat ("\n")
if (coint==1){
cat ("Cointegrating Vector ",b0," ",seb,"\n")
}else{
cat ("Cointegrating Vector ",b0,"\n")
}
cat ("Negative Log-Like    ",nlike,"\n")
cat ("AIC                  ",aic,"\n")
cat ("BIC                  ",bic,"\n")
cat ("\n")
cat ("  Equation 1","\n")
tbeta0 <- format(beta0,digits=4)
tse <- format(se,digits=4)
for (j in 1:kk) cat(" ",tbeta0[j,1],"  ",tse[j],"\n")
cat ("\n")
cat ("  Equation 2","\n")
for (j in 1:kk) cat(" ",tbeta0[j,2],"  ",tse[kk+j],"\n")
cat ("\n")
cat ("\n")
cat ("\n")
}

# Set-up Grids #
q <- unique(w0)
q <- as.matrix(sort(q))
gamma1 <- q[round(seq(trim,trim+(gn-1)*(1-2*trim)/gn,(1-2*trim)/gn)*nrow(q))]
gamma2 <- q[round(seq(1/(gn+1),1/(gn+1)*gn,1/(gn+1))*nrow(q))]
if (coint==1){
  if (bn==1){ 
      betas <- b0
  }else{
      betas <- seq((b0-6*seb),(b0-6*seb+(bn-1)*12*seb/(bn-1)),12*seb/(bn-1))
  }
}else{
  betas <- b0
  bn <- 1
}

# Estimate Threshold Model #

store <- matrix(0,gn,bn)+nlike
for (j in 1:gn){
  gam <- gamma2[j]
  for (bj in 1:bn){
    beta <- betas[bj]
    w <- xlag%*%rbind(1,(-beta))
    z <- cbind(w,x)
    d1 <- (w <= gam)
    n1 <- sum(d1)
    if (min(rbind(n1,(t-n1)))/t > trim){
      zj <- cbind((z*(d1%*%matrix(1,1,ncol(z)))),w)
      zzj <- as.matrix(zj-x%*%xx%*%(t(x)%*%zj))
      if (qr(zzj)$rank==ncol(zzj)) bz <- qr.solve(zzj,u)
      if (qr(zzj)$rank<ncol(zzj)) bz <- (qr(t(zzj)%*%zzj)$qr)%*%(t(zzj)%*%u)
      store[j,bj] <- (t/2)*log(det(t(u-zzj%*%bz)%*%(u-zzj%*%bz)/t))
    }
  }
}
m <- apply(store,2,which.min)
c <- which.min(diag(as.matrix(store[m,])))
r <- m[c]
gammahat <- gamma2[r]
b1 <- betas[c]
nlike <- store[r,c]
bic <- nlike+log10(t)*8*(1+k)
aic <- nlike+2*8*(1+k)

  w <- xlag%*%rbind(1,-b1)
  z <- cbind(w,x)
  d1 <- as.matrix((w <= gammahat))
  d2 <- 1-d1
  y1 <- as.matrix(y[d1%*%matrix(1,1,ncol(y))>0])
  y1 <- matrix(y1,nrow(y1)/ncol(y),ncol(y))
  z1 <- as.matrix(z[d1%*%matrix(1,1,ncol(z))>0])
  z1 <- matrix(z1,nrow(z1)/ncol(z),ncol(z))
  y2 <- as.matrix(y[d2%*%matrix(1,1,ncol(y))>0])
  y2 <- matrix(y2,nrow(y2)/ncol(y),ncol(y))
  z2 <- as.matrix(z[d2%*%matrix(1,1,ncol(z))>0])
  z2 <- matrix(z2,nrow(z2)/ncol(z),ncol(z))
  zz1 <- solve(t(z1)%*%z1)
  zz2 <- solve(t(z2)%*%z2)
  beta1 <- zz1%*%(t(z1)%*%y1) 
  beta2 <- zz2%*%(t(z2)%*%y2)
  e1 <- y1 - z1%*%beta1
  e2 <- y2 - z2%*%beta2
  if (cov==1){
    xe1 <- cbind((z1*(e1[,1]%*%matrix(1,1,ncol(z1)))),(z1*(e1[,2]%*%matrix(1,1,ncol(z1)))))
    m1 <- k_product(diag(2),zz1)
    v1 <- m1%*%t(xe1)%*%xe1%*%m1  
    xe2 <- cbind((z2*(e2[,1]%*%matrix(1,1,ncol(z2)))),(z2*(e2[,2]%*%matrix(1,1,ncol(z2)))))
    m2 <- k_product(diag(2),zz2)
    v2 <- m2%*%t(xe2)%*%xe2%*%m2  
  }else{
    sig1 <- t(e1)%*%e1/nrow(e1)
    v1 <- k_product(sig1,zz1)
    sig2 <- t(e2)%*%e2/nrow(e2)
    v2 <- k_product(sig2,zz2)
  }
  se1 <- as.matrix(sqrt(diag(v1)))
  se2 <- as.matrix(sqrt(diag(v2))) 
  rb <- nrow(beta1)
  bb <- beta1-beta2
  if (k>0){
    bw <- bb[3:rb,]
    bw <- matrix(bw,ncol(bw)*nrow(bw),1)
    vr <- rbind(as.matrix(seq(3,rb,1)),as.matrix(seq(rb+3,2*rb,1)))
    vv <- v1[vr,vr]+v2[vr,vr]
    ww <- t(bw)%*%solve(vv)%*%bw
  }
  indx <- rbind(1,rb+1)
  wecm <- bb[1,]%*%solve(v1[indx,indx]+v2[indx,indx])%*%(as.matrix(bb[1,]))

if (p_ests==1){
cat ("Threshold VECM Estimates", "\n")
cat ("\n")
cat ("Threshold Estimate            ", gammahat, "\n")
cat ("Cointegrating Vector Estimate ", b1, "\n")
cat ("Negative Log-Like             ", nlike, "\n")
cat ("AIC                           ", aic, "\n")
cat ("BIC                           ", bic, "\n")
cat ("\n")
cat ("First Regime", "\n")
cat ("Percentage of Obs", mean(d1), "\n")
cat ("\n")
cat ("  Equation 1","\n")
tbeta1 <- format(beta1,digits=4)
tse1 <- format(se1,digits=4)
for (j in 1:kk) cat(" ",tbeta1[j,1],"  ",tse1[j],"\n")
cat ("\n")
cat ("  Equation 2","\n")
for (j in 1:kk) cat(" ",tbeta1[j,2],"  ",tse1[kk+j],"\n")
cat ("\n")
cat ("Second Regime", "\n")
cat ("Percentage of Obs", mean(d2),"\n")
cat ("\n")
cat ("  Equation 1","\n")
tbeta2 <- format(beta2,digits=4)
tse2 <- format(se2,digits=4)
for (j in 1:kk) cat(" ",tbeta2[j,1],"  ",tse2[j],"\n")
cat ("\n")
cat ("  Equation 2","\n")
for (j in 1:kk) cat(" ",tbeta2[j,2],"  ",tse2[kk+j],"\n")
cat ("\n")
if (k>0){
cat ("Wald Test for Equality of Dynamic Coefs ", ww," ",1-pchisq(ww,(rb-2)*2),"\n")
}
cat ("Wald Test for Equality of ECM Coef      ", wecm," ",1-pchisq(wecm,2),"\n")
cat ("\n")
cat ("\n")
}

#***************************************************************#

if (graph==1 | graph_rotate==1){

if (coint==1){
  mtit <- rbind("Figure 2","Concentrated Negative Log Likelihood")
  xtit <- "Beta"
  ytit <- "Negative Log-Likelihood"
  loglike <- apply(store,2,min)   
  x11()
  plot(betas,loglike,type="l",ann=0)
  title(main=mtit,ylab=ytit,xlab=xtit)
  if (graph_rotate==1){
      x11()
      plot(loglike,betas,type="l",ann=0)
      title(main=mtit,ylab=xtit,xlab=ytit)
  }
}

mtit <- rbind("Figure 1","Concentrated Negative Log Likelihood")
xtit <- "Gamma"
ytit <- "Negative Log-Likelihood"
loglike <- apply(t(store),2,min) 
x11()
plot(gamma2,loglike,type="l",ann=0)
title(main=mtit,ylab=ytit,xlab=xtit)
if (graph_rotate==1){
    x11()
    plot(loglike,gamma2,type="l",ann=0)
    title(main=mtit,ylab=xtit,xlab=ytit)
}

yhat2 <- (z[,1:2]%*%beta1[1:2,])*(d1%*%matrix(1,1,2))+(z[,1:2]%*%beta2[1:2,])*(d2%*%matrix(1,1,2))
wi <- order(w)
yhat2 <- yhat2[wi,]
wi <- w[wi]
ylong <- paste(c("R_"),long,sep="")
yshort <- paste(c("R_"),short,sep="")
ylong_t <- paste(ylong,c("(t-1)"),sep="") 
yshort_t <- paste(yshort,c("(t-1)"),sep="") 
if (coint==0) xtit <- paste(c("Error Correction: "),ylong_t,c("-"),yshort_t,sep="")  
if (coint==1) xtit <- paste(c("Error Correction: "),ylong_t,c("-beta*"),yshort_t,sep="") 
ytit <- "Interest Rate Response"
mtit <- rbind("Figure 3","Interest Rate Response to Error-Correction")
x11()
xxlim <- range(wi)
yylim <- range(yhat2)
plot(wi,yhat2[,1],lty=1,col=1,xlim=xxlim,ylim=yylim,type="l",ann=0)
lines(wi,yhat2[,2],lty=2,col=2)
title(main=mtit,ylab=ytit,xlab=xtit)
tit1 <- paste(c("Long Rate ("),ylong,c(")"),sep="")
tit2 <- paste(c("Short Rate ("),yshort,c(")"),sep="")
legend("bottomright",c(tit1,tit2),lty=c(1,2),col=c(1,2))

if (graph_rotate==1){
  x11()
  plot(yhat2[,1],wi,lty=1,col=1,xlim=yylim,ylim=xxlim,type="l",ann=0)
  lines(yhat2[,2],wi,lty=2,col=2)
  title(main=mtit,ylab=xtit,xlab=ytit)
  legend("topleft",c(tit1,tit2),lty=c(1,2),col=c(1,2))
}
}

#********************************************************#

# Function for LM test #

lmtest01 <- function(y,x,w0,gammas,plot){
z0 <- cbind(w0,x)
zz0 <- t(z0)%*%z0
if (qr(zz0)$rank==ncol(zz0)) z0zz <- t(solve(zz0)%*%t(z0))
if (qr(zz0)$rank<ncol(zz0)) z0zz <- z0%*%(qr(zz0)$qr)
e <- y-z0zz%*%(t(z0)%*%y) 
e1 <- e[,1]
e2 <- e[,2]
store <- matrix(0,gn,1)
for (j in 1:gn){
    d1 <- (w0 <= gammas[j])
    n1 <- sum(d1)
    if (min(rbind(n1,(t-n1)))/t > trim){
      z1 <- z0*(d1%*%matrix(1,1,ncol(z0)))
      z11 <- z1-z0zz%*%(t(z0)%*%z1)
      ze <- cbind((z11*(e1%*%matrix(1,1,ncol(z11)))),(z11*(e2%*%matrix(1,1,ncol(z11)))))
      v <- t(ze)%*%ze 
      s <- t(z11)%*%y
      s <- matrix(s,ncol(s)*nrow(s),1)
      if (qr(v)$rank==ncol(v)) sv <- solve(v)%*%s
      if (qr(v)$rank<ncol(v)) sv <- (qr(v)$qr)%*%s
      store[j] <- t(s)%*%sv
    }
}
lm01 <- max(store)
if (plot==1){
    xtit <- "Gamma"
    mtit <- "LM Statistic as function of Gamma"
    x11()
    plot(gammas,store,type="l",ann=0)
    title(main=mtit,xlab=xtit)
    if (graph_rotate==1){
        plot(store,gammas,type="l",ann=0)
        title(main=mtit,ylab=xtit)
    }
}
lm01
}

# Testing #

lm01 <- lmtest01(y,x,w0,gamma1,graph)
if (p_ests==1){
cat ("Lagrange Multipler Threshold Test", "\n")
cat ("\n")
cat ("Test Statistic                   ", lm01,"\n")
cat ("\n")
}

if (boot>0){
# Fixed Regressor Bootstrap #
fix01 <- matrix(0,boot,1)
for (r in 1:boot){
  yr <- matrix(rnorm(2*t),t,2)*e 
  fix01[r] <- lmtest01(yr,x,w0,gamma1,0)
}
pfix01 <- mean(fix01>lm01)
fix01 <- sort(fix01)
cfix01 <- fix01[round(.95*boot)]

# Parametric Bootstrap #
x0 <- dat[1:(1+k),]
mu <- beta0[2,]
ab <- rbind(1,(-b0))%*%beta0[1,]
boot01 <- matrix(0,boot,1)
if (k>0){
  capgamma <- beta0[3:(2+2*k),]
  dx0 <- (t(apply((x0[2:(k+1),]-x0[(1:k),]),2,rev)))
  dx0 <- matrix(dx0,ncol(dx0)*nrow(dx0),1)
}
for (r in 1:boot){
  datb <- x0
  if (k>0) dx <- dx0   
  datbi <- as.matrix(datb[(k+1),])
  ei <- as.matrix(ceiling(runif(t)*t))  
  eb <- e[ei,]
  for (i in (k+2):n){
    u <- mu+eb[(i-k-1),]+datbi%*%ab
    if (k>0) u <- u+t(dx)%*%capgamma
    datbi <- datbi+u
    datb <- rbind(datb,datbi)
    if (k==1) dx <- t(u)
    if (k>1) dx <- rbind(t(u),as.matrix(dx[1:(2*k-2)]))
  }
  yr <- as.matrix(datb[(2+k):n,]-datb[(1+k):(n-1),])
  xlagr <- as.matrix(datb[(1+k):(n-1),])
  xr <- matrix(1,t,1)
  for (j in 1:k) xr <- cbind(xr,(datb[(2+k-j):(n-j),]-datb[(1+k-j):(n-1-j),]))
  xr <- as.matrix(xr) 
  if (coint==1){
    xxr <- solve(t(xr)%*%xr) 
    u <- yr-xr%*%xxr%*%(t(xr)%*%yr)
    v <- xlagr-xr%*%xxr%*%(t(xr)%*%xlagr)
    m <- solve(t(v)%*%v)%*%(t(v)%*%u)%*%solve(t(u)%*%u)%*%(t(u)%*%v)
    ev <- eigen(m)
    va <- ev$val
    ve <- ev$vec
    h <- ve[,which.max(va)]
    wr <- as.matrix(xlagr%*%rbind(1,(h[2]/h[1])))
  }else{
    wr <- as.matrix(xlagr%*%rbind(1,(-cvalue)))
  }
  boot01[r] <- lmtest01(yr,xr,wr,gamma1,0)
}
boot01 <- sort(boot01)
cb01 <- boot01[round(.95*boot)]
pb01 <- mean(boot01>lm01)

if (p_ests==1){

cat ("Fixed Regressor (Asymptotic) 5% Critical Value ", cfix01,"\n")
cat ("Bootstrap 5% Critical Value                    ", cb01,"\n")
cat ("Fixed Regressor (Asymptotic) P-Value           ", pfix01,"\n")
cat ("Bootstrap P-Value                              ", pb01,"\n")
cat ("\n")
cat ("\n")
cat ("\n")
  }
 }
}









