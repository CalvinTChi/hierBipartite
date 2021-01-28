#' Sparse canonical covariance analysis
#'
#' 'scca' is used to perform sparse canonical covariance analysis (SCCA)
#'
#' @usage scca(X,Y,penalty="HL",lamx=c(1,2,3),lamy=c(1,2,3),nc=1,
#' tuning="CV.alt",K=5,seed=NULL,center=TRUE,scale=FALSE)
#'
#' @param X n-by-p data matrix, where n is the number of subjects and p is the
#' number of variables
#' @param Y n-by-q data matrix, where q is the number of variables
#' @param penalty "HL" is the unbounded penalty proposed by Lee and Oh (2009).
#' "LASSO" (Tibshirani, 1996), "SCAD" (Fan and Li, 2001) and "SOFT"
#' (soft thresholding) are also available as other penalty options.
#' Default is "HL".
#' @param lamx A vector specifying grid points of the tuning parameter for X.
#' Default is (1,2,3).
#' @param lamy A vector specifying grid points of the tuning parameter for Y.
#' Default is (1,2,3).
#' @param nc Number of components (canonical vectors). Default is 1.
#' @param tuning How to find optimal tuning parameters for the sparsity.
#' If tuning="CV.full", then the tuning parameters are selected
#' automatically via K-fold cross-validation by using 2-dim'l grid search.
#' If "CV.alt", then a sequential 1-dim'l search method is applied instead of
#' the 2-dim'l grid search. Default is "CV.alt".
#' @param K Perform K-fold cross-validation.
#' @param seed Seed number for initialization. A random initial point is
#' generated for tuning="CV.alt".
#' @param center The columns of the data matrix are centered to have mean zero.
#' Default is TRUE.
#' @param scale The columns of the data matrix are scaled to have variance 1.
#' Default is FALSE.
#' @details Sparse CCA uses a random-effect model approach to obtain sparse
#' regression. This model gives unbounded gains for zero loadings at the
#' origin. Various penalty functions can be adapted as well.
#' @return
#' \itemize{
#'   \item A: p-by-nc matrix, k-th colum of A corresponds to k-th pattern
#'   \item B: q-by-nc matrix, k-th colum of B corresponds to k-th pattern
#'   (canonical vector) for Y
#'   \item U: n-by-nc matrix. k-th column of U corresponds to k-th score
#'   associated with k-th pattern for X
#'   \item V: n-by-nc matrix. k-th column of V corresponds to k-th score
#'   associated with k-th pattern for Y
#'   \item lambda: nc-by-2 matrix. k-th row of lambda corresponds to the optimal
#'   tuning parameters for k-th pattern pairs
#'   \item CR: average cross-validated sample covariance
#' }
#' @author Woojoo Lee, Donghwan Lee, Youngjo Lee and Yudi Pawitan
#' @references Lee, W., Lee, D., Lee, Y. and Pawitan, Y. (2011) Sparse Canonical
#' Covariance Analysis for High-throughput Data
#' @examples
#' ## Example 1
#' ## A very simple simulation example
#' n<-10; p<-50; q<-20
#' X = matrix(rnorm(n*p),ncol=p)
#' Y = matrix(rnorm(n*q),ncol=q)
#' scca(X,Y)
#'
#' @export
scca <-
  function(X, Y, penalty="HL", lamx=c(1,2,3),lamy=c(1,2,3), nc=1, tuning="CV.alt",K=5, seed=NULL, center=TRUE, scale=FALSE){

    if (penalty=="SOFT" && min(lamx)>=1) {stop("Range of lamx for SOFT should be (0,1)")}
    if (penalty=="SOFT" && min(lamy)>=1) {stop("Range of lamy for SOFT should be (0,1)")}

    X<-scale(X,center=center,scale=scale);Y<-scale(Y,center=center,scale=scale)

    if (is.null(seed)) {seed<-Sys.time()}
    if (nrow(X) != nrow(Y)) {stop("X and Y should have same number of rows")}

    svd.X <- svd(X)

    n<-dim(X)[1]; p<-dim(X)[2]; q<-dim(Y)[2]

    U<-matrix(0,n,nc)
    V<-matrix(0,n,nc)
    A<-matrix(0,p,nc)
    B<-matrix(0,q,nc)
    CR<-rep(0,nc)

    L<-matrix(0,nc,2)

    X.new<-X; Y.new<-Y

    for (cnt in 1:nc){

      if (tuning[1]=="CV.full") {cv.full<-opt.cv.full(X.new,Y.new,lamx,lamy,K=K,penalty=penalty,seed=seed); olamx<-cv.full$optx; olamy<-cv.full$opty; CR1<-cv.full$avecvcov}
      if (tuning[1]=="CV.alt") {cv.alt<-opt.cv.alt(X.new,Y.new,lamx,lamy,K=K,penalty=penalty,seed=seed); olamx<-cv.alt$optx; olamy<-cv.alt$opty; CR1<-cv.alt$avecvcov}

      if (is.numeric(tuning)) {olamx<-tuning[cnt,1];olamy<-tuning[cnt,2]}


      if (penalty!="SOFT") {sscca.rst<-NIPALS.sparse(X.new,Y.new,olamx, olamy, penalty=penalty)}
      if (penalty=="SOFT") {sscca.rst<-NIPALS.soft(X.new,Y.new,olamx, olamy)}


      X.new<-X.new- sscca.rst$u1 %*%t(sscca.rst$a1)
      Y.new<-Y.new- sscca.rst$v1 %*%t(sscca.rst$b1)


      ytx.new<-t(Y.new)%*%X.new
      L[cnt,1]<-olamx; L[cnt,2]<-olamy

      cat("Computing Component number ",cnt,"\n")

      U[,cnt]<-sscca.rst$u1; V[,cnt]<-sscca.rst$v1
      A[,cnt]<-sscca.rst$a1; B[,cnt]<-sscca.rst$b1
      CR[cnt]<-CR1

    }

    return(list(U=U, V=V, A=A, B=B, lambda=L,CR=CR))
  }

opt.cv.alt <-
  function(X,Y,K, lamx, lamy, penalty,seed){

    row.x<-dim(X)[1]
    seq.x<-lamx
    seq.y<-lamy

    i.x<-length(seq.x)
    i.y<-length(seq.y)
    cvcov<-rep(0,K)
    cvcovdiff<-rep(0,K)
    avecvcov<-matrix(0,i.x,i.y)
    avecvcovdiff<-matrix(0,i.x,i.y)

    set.seed(seed)

    cv.set<-split(sample(1:row.x),rep(1:K,length=row.x))

    opt.i<-sample(c(1:i.x),1)

    opt.conv<-1;max.old<-0;opt.iter<-0

    while(opt.conv>1e-02 & opt.iter<20){

      for (jj in 1:i.y){

        for (kk in 1:K){
          XX<-X[-cv.set[[kk]],]
          YY<-Y[-cv.set[[kk]],]
          if (penalty=="LASSO") {resXY<-NIPALS.sparse(XX,YY,seq.x[opt.i],seq.y[jj], "LASSO")}
          if (penalty=="SCAD") {resXY<-NIPALS.sparse(XX,YY,seq.x[opt.i],seq.y[jj], "SCAD")}
          if (penalty=="HL") {resXY<-NIPALS.sparse(XX,YY,seq.x[opt.i],seq.y[jj], "HL")}
          if (penalty=="SOFT") {resXY<-NIPALS.soft(XX,YY,seq.x[opt.i],seq.y[jj])}

          ra1<-resXY$a1
          rb1<-resXY$b1

          cvcov[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))
          cvcovdiff[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))-abs(resXY$rho1)
        }

        avecvcov[opt.i,jj]<-sum(abs(cvcov))/K
        avecvcovdiff[opt.i,jj]<-sum(abs(cvcovdiff))/K
      }

      opt.j<- which.max(avecvcov[opt.i,])

      for (ii in 1:i.x){

        for (kk in 1:K){
          XX<-X[-cv.set[[kk]],]
          YY<-Y[-cv.set[[kk]],]
          if (penalty=="LASSO") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[opt.j], "LASSO")}
          if (penalty=="SCAD") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[opt.j], "SCAD")}
          if (penalty=="HL") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[opt.j], "HL")}
          if (penalty=="SOFT") {resXY<-NIPALS.soft(XX,YY,seq.x[ii],seq.y[opt.j])}

          ra1<-resXY$a1
          rb1<-resXY$b1

          cvcov[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))
          cvcovdiff[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))-abs(resXY$rho1)
        }

        avecvcov[ii,opt.j]<-sum(abs(cvcov))/K
        avecvcovdiff[ii,opt.j]<-sum(abs(cvcovdiff))/K
      }
      opt.i<- which.max(avecvcov[,opt.j])

      max.new<-max(avecvcov[opt.i,opt.j])
      opt.conv<-abs(max.new-max.old)
      max.old<-max.new
      opt.iter<-opt.iter+1
    }


    r.index<-opt.i
    c.index<-opt.j

    return(list(optx=seq.x[r.index],opty=seq.y[c.index],avecvcov=avecvcov[r.index,c.index]))
  }

opt.cv.full <-
  function(X,Y,K, lamx,lamy,penalty,seed){

    row.x<-dim(X)[1]
    seq.x<-lamx;  seq.y<-lamy

    i.x<-length(seq.x)
    i.y<-length(seq.y)
    cvcov<-rep(0,K)
    cvcovdiff<-rep(0,K)
    avecvcov<-matrix(0,i.x,i.y)
    avecvcovdiff<-matrix(0,i.x,i.y)

    set.seed(seed)

    cv.set<-split(sample(1:row.x),rep(1:K,length=row.x))


    for (ii in 1:i.x){
      for (jj in 1:i.y){

        for (kk in 1:K){
          XX<-X[-cv.set[[kk]],]
          YY<-Y[-cv.set[[kk]],]
          if (penalty=="LASSO") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[jj], "LASSO")}
          if (penalty=="SCAD") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[jj], "SCAD")}
          if (penalty=="HL") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[jj], "HL")}
          if (penalty=="SOFT") {resXY<-NIPALS.soft(XX,YY,seq.x[ii],seq.y[jj])}

          ra1<-resXY$a1
          rb1<-resXY$b1

          cvcov[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))
          cvcovdiff[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))-abs(resXY$rho1)
        }

        avecvcov[ii,jj]<-sum(abs(cvcov))/K
        avecvcovdiff[ii,jj]<-sum(abs(cvcovdiff))/K
      }
    }
    r.index<-(which.max(avecvcov)%%i.x)*((which.max(avecvcov)%%i.x)>0)+i.x*((which.max(avecvcov)%%i.x)==0)
    c.index<-(floor(which.max(avecvcov)/i.x)+1)*((which.max(avecvcov)%%i.x)>0)+(floor(which.max(avecvcov)/i.x))*((which.max(avecvcov)%%i.x)==0)

    r.indexw<-(which.min(avecvcovdiff)%%i.x)*((which.min(avecvcovdiff)%%i.x)>0)+i.x*((which.min(avecvcovdiff)%%i.x)==0)
    c.indexw<-(floor(which.min(avecvcovdiff)/i.x)+1)*((which.min(avecvcovdiff)%%i.x)>0)+(floor(which.min(avecvcovdiff)/i.x))*((which.min(avecvcovdiff)%%i.x)==0)
    return(list(optx=seq.x[r.index],opty=seq.y[c.index],optxw=seq.x[r.indexw],optyw=seq.y[c.indexw],avecvcov=avecvcov[r.index,c.index],avecvcovdiff=avecvcovdiff))
  }

NIPALS.sparse <-
  function(X,Y,lamx,lamy, penalty){
    eps.cov<-1
    p<-dim(X)[2]
    q<-dim(Y)[2]
    v1<-Y[,1]

    res0<-NIPALS(X,Y)
    a1<-res0$a1
    b1<-res0$b1
    niter<-0
    a<-3.7
    w<-30
    if (lamx==0) {sigx<-median(sd(X))}
    if (lamy==0) {sigy<-median(sd(Y))}

    if (lamx!=0) {sigx<-sd(res0$a1)/sqrt(2)}
    if (lamy!=0) {sigy<-sd(res0$b1)/sqrt(2)}

    if (lamx==0 & penalty=="SCAD") {lamx<-1e-10}
    if (lamy==0 & penalty=="SCAD") {lamy<-1e-10}

    while ((eps.cov>1e-03)&(niter<100)){
      eps.cov.a1<-1
      nitera<-0
      while ((eps.cov.a1>1e-03)&(nitera<50)){

        if (penalty == "LASSO") {uux<-abs(a1)+1e-08}
        if (penalty == "SCAD") {uux<-(abs(a1)+1e-08)/ (as.numeric(abs(a1)<=lamx)+as.numeric(abs(a1)>lamx)*(a*lamx-abs(a1))*as.numeric(a*lamx>abs(a1))/((a-1)*lamx))}
        if (penalty == "HL"){    kax <- sqrt(4*a1^2/(w*sigx^2)+((2/w)-1)^2) ;    uux <- 0.25*w*(((2/w)-1)+kax)+1e-08 }

        WWx<-as.vector(1/uux)
        a1.old<-a1

        v1.norm<-sum(v1^2)
        Xtv1<-crossprod(X, v1)

        for (i in 1:p){
          a1[i]<-(1/(v1.norm+lamx*WWx[i]))*Xtv1[i]
        }
        eps.cov.a1<-max(abs(a1-a1.old))
      }
      a1.norm<-sum(a1^2)
      a1<-a1/sqrt(a1.norm)
      u1<-X%*%a1


      eps.cov.b1<-1
      niterb<-0
      while ((eps.cov.b1>1e-03)&(niterb<50) ){

        if (penalty=="LASSO") {uuy<-abs(b1)+1e-08}
        if (penalty=="SCAD") {uuy<-(abs(b1)+1e-08)/ (as.numeric(abs(b1)<=lamy)+as.numeric(abs(b1)>lamy)*(a*lamy-abs(b1))*as.numeric(a*lamy>abs(b1))/((a-1)*lamy))}
        if (penalty=="HL"){    kay <- sqrt(4*b1^2/(w*sigy^2)+((2/w)-1)^2) ;    uuy <- 0.25*w*(((2/w)-1)+kay)+1e-08 }

        WWy<-as.vector(1/uuy)
        b1.old<-b1

        u1.norm<-sum(u1^2)
        Ytu1<-crossprod(Y, u1)
        for (i in 1:q){
          b1[i]<-(1/(u1.norm+lamy*WWy[i]))*Ytu1[i]
        }
        eps.cov.b1<-max(abs(b1-b1.old))
      }


      b1.norm<-sum(b1^2)
      b1<-b1/sqrt(b1.norm)

      v1.old<-v1
      v1<-Y%*%b1
      eps.cov<-max(abs(v1-v1.old))
      niter<-niter+1
      if ((sum(round(a1,digits=4)==0)>=(p-2))&(sum(round(b1,digits=4)==0)>=(q-2))) {break}
    }

    a1<-a1/sqrt(sum(a1^2))
    b1<-b1/sqrt(sum(b1^2))
    u1<-X%*%a1
    v1<-Y%*%b1
    rho1<-crossprod(u1, v1)
    if (niter>=100) {print("No Convergence")}
    return(list(rho1=rho1, u1=u1, v1=v1, a1=round(a1, digits=4),
                b1=round(b1, digits=4), WWx=WWx, WWy=WWy, niter=niter))

  }

NIPALS.soft <-
  function(X,Y,lamx,lamy){
    eps.cov <- 1
    p <- dim(X)[2]
    q <- dim(Y)[2]
    n <- dim(X)[1]
    v1 <- Y[,1]
    res.NIPALS <- NIPALS(X,Y)
    a1 <- res.NIPALS$a1
    b1 <- res.NIPALS$b1
    niter <- 0

    while ((eps.cov > 10^-3) & (niter < 100)){

      a1.old <- a1

      v1.norm <- sum(v1^2)
      Xtv1 <- crossprod(X, v1)

      for (i in 1:p){
        a1[i] <- soft((Xtv1[i] / v1.norm),lamx)
      }

      if (sum(a1==0) < p) {a1 <- a1 / sqrt(sum(a1^2))}
      if (sum(a1==0) == p) {
        u1 <- rep(0,n);
        b1 <- rep(0,q);
        v1 <- rep(0,n);
        break
      }
      u1 <- X%*%a1

      b1.old <- b1

      u1.norm <- sum(u1^2)
      Ytu1 <- crossprod(Y, u1)
      for (i in 1:q) {
        b1[i] <- soft(Ytu1[i] / u1.norm, lamy)
      }

      if (sum(b1==0) < q) {
        b1 <- b1 / sqrt(sum(b1^2))
      }
      if (sum(b1 == 0) == q) {
        u1 <- rep(0,n);
        a1 <- rep(0,p);
        v1 <- rep(0,n);
        break
      }
      v1.old <- v1
      v1 <- Y%*%b1
      eps.cov <- max(abs(v1-v1.old))
      niter <- niter+1
    }

    rho1 <- t(u1)%*%v1
    return(list(rho1=rho1, u1=u1, v1=v1, a1=round(a1, digits=4),
                b1=round(b1, digits=4), niter=niter))
  }

NIPALS <-
  function(X,Y){
    eps.cov<-1
    v1<-Y[,1]

    while (eps.cov>10^-6){
      a1<-crossprod(X,v1)/sum(v1^2)
      a1<-a1/sqrt(sum(a1^2))
      u1<-X%*%a1

      b1<-crossprod(Y,u1)/sum(u1^2)
      b1<-b1/sqrt(sum(b1^2))
      v1.old<-v1
      v1<-Y%*%b1
      eps.cov<-max(abs(v1-v1.old))
    }
    rho1<-crossprod(u1, v1)
    return(list(rho1=rho1,u1=u1,v1=v1,a1=round(a1,digits=4),b1=round(b1,digits=4)))
  }

soft <-
  function(y,thr){
    sign(y)*(abs(y)-thr)*(abs(y)>thr)
  }

