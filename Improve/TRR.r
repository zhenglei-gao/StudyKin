## a collection of the trust region methods

STIR_LS <- function()
{
  ## Steps 1: Consider the problem of finding the trial step size s_i between the ith point x_i and the next 
  ## point x_i+1 that minimize \Psi_i()s
  
  
  ## Step 2: Compute the 2 dimentional subspace V with 2 spanning vectors: one vector in the direction of the 
  ## gradient at x_i, the other in the approximate Gauss-Newto direction
  
  ## Step 3: Solve for the trial step s_i of the 2 dimentional subproblem
  
  ## Step 4: if for a predefined constant \tau in (0,1), u(x_i+s_i) < u(x_i), then x_i+1 becomes the current point,
  ## otherwise x_i+1=x_i
  
  ## Step 5: update \Delat_i
  
  ## Step 6, if \partial u(s_i) is below a chosen tolerance, end, o.w, repeat and increment i
  
  
}


trr.options <- function(nopar=1,gradflag=0,pcflags=Inf,MaxPCGIter=1,TolPCG=0.1,MaxFunEvals=100,MaxIter=500,TolFun=1e-6,TolX=1e-6)
{
  MaxPCGIter <- max(MaxPCGIter,nopar/2)
  MaxFunEvals <- max(MaxFunEvals,100*nopar)
  return(list(gradflag=gradflag,pcflags=pcflags,kmax=MaxPCGIter,
              pcgtol=TolPCG,maxfunevals=MaxFunEvals,itb=MaxIter,
              tol1=TolFun,tol2=TolX))
}

seceqn <- function(lambda,eigval,alpha,delta)
{
  m <- length(lambda)
  n <- length(eigval)
  unn <- rep(1,n)
  unm <- rep(1,m)
  M <- eigval%*%t(unm)+unn%*%t(lambda)
  MC <- M
  MM <- alpha%*%t(unm)
  M[M!=0] <- MM[M!=0]/M[M!=0]
  M[MC==0] <- Inf
  M <- M*M
  value <- sqrt(unm/(t(M)%*%unn))
  value[is.nan(value)] <- 0
  value <- (1/delta)*unm-value
  return(value)
}
rfzero <- function(f,...,lower=0,tol=1e-6,maxiter=100)
{
  stop("cannot be used still! Error inside!!")
  itfun <- 0
  x <- lower
  if(x!=0) dx <- abs(x)/2 else dx <- 0.5
  a <- x
  c <- a
  fa <- f(a,...)
  itfun <- itfun+1
  b <- x + dx
  ##b <- x+1
  fb <- f(b,...)
  itfun <- itfun +1
  while((fa>0) * (fb>0) ==1){
    dx <- 2*dx
    if((fa>0) != (fb>0)){
      break
    } else{
      b <- x + dx
      ##b <- x+1
      fb <- f(b,...)
      itfun <- itfun +1
      if(itfun>maxiter) break
    }
    
      
  }
  fc <- fb
  while(fb!=0){
    if((fb>0)* (fc>0) ==1){
      c <- a
      fc <- fa
      d <- b-a
      e <-d
    }
    if(abs(fc)<abs(fb)){
      a<-b
      b<-c
      c<-a
      fa <- fb
      fb <- fc
      fc <- fa
    }
    
    if(itfun>itbnd) break
    m <- 0.5*(c-b)
    toler <- 2.0*tol*max(abs(b),1)
    if(abs(m)<=toler || fb==0) break
    if(abs(e)<toler || (abs(fa)<=abs(fb))){
      d <- m
      e <- m
    }else{
      s <- fb/fa
      if(a==c){
        p <- 2.0*m*s
        q <- 1.0-s
      }else{
        q <- fa/fc
        r <- fb/fc
        p <- s*(2*m*q*(q-r)-(b-a)*(r-1))
        q <- (q - 1.0)*(r - 1.0)*(s - 1.0)
      }
      if(p>0) q <- -q else p <- -p  
      if((2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))){
        e <- d
        d <- p/q
      }else{
        d <- m
        e <- m
      }
    }
    a <- b
    fa <- fb
    if(abs(d)>toler) b <- b+d else{
      if(b>c) b <- b-toler else b <- b+toler
    }
    fb <- f(b,...)
    itfun <- itfun+1
    
  }
  #browser()
  return(list(b=b,c=c,itfun=itfun))
}

trust.exact <- function(g,H,delta)
{
  #   %TRUST  Exact soln of trust region problem
  #   %
  #   % [s,val,posdef,count,lambda] = TRUST(g,H,delta) Solves the trust region
  #   % problem: min{g^Ts + 1/2 s^THs: ||s|| <= delta}. The full
  #   % eigen-decomposition is used; based on the secular equation,
  #   % 1/delta - 1/(||s||) = 0. The solution is s, the value
  #   % of the quadratic at the solution is val; posdef = 1 if H
  #   % is pos. definite; otherwise posdef = 0. The number of
  #   % evaluations of the secular equation is count, lambda
  #   % is the value of the corresponding Lagrange multiplier.
  #   %
  #   %
  #   % TRUST is meant to be applied to very small dimensional problems.
  norm <- function(x) sqrt(sum(x^2))
  tol <- 1e-12
  tol2 <- 1e-8
  key <- 0
  itbnd <- 50
  lambda <- 0
  n <- length(g)
  coeff <- rep(0,n)
  tmp <- eigen(H,symmetric=TRUE)
  V <- tmp$vectors
  D <- tmp$values
  count <- 0
  eigval <- D
  jmin <- which.min(eigval)
  mineig <- eigval[jmin]
  alpha <- -t(V)%*%g
  sig <- sign(alpha[jmin])+(alpha[jmin]==0)
  
  ## positive definite case
  if(mineig > 0){
    coeff <- alpha / eigval; 
    lambda <- 0;
    s <- V*coeff; 
    posdef <- 1; 
    nrms <- norm(s);
    if(nrms <= 1.2*delta) key <-1 else laminit <- 0
  }else{
    laminit <- -mineig
    posdef <- 0
  }
  
  if(key==0){
    tmp <- seceqn(laminit,eigval,alpha,delta)
    if(tmp>0){
      ##tmp1 <- rfzero(f=seceqn,eigval=eigval,alpha=alpha,delta=delta,lower=laminit,tol=tol,maxiter=itbnd)
      ##tmp1 <- optim(par=laminit,fn=seceqn^2,eigval=eigval,alpha=alpha,delta=delta,lower=laminit,control=list(tol=tol,maxiter=itbnd))
      tmp1 <- uniroot(f=seceqn,eigval=eigval,alpha=alpha,delta=delta,lower=laminit,upper=1e100,tol=tol,maxiter=itbnd)
      b <- tmp1$root
      count <- tmp1$iter
      vval <- abs(seceqn(b,eigval,alpha,delta))
      if(vval <= tol2){
        lambda <- b
        key <- 2
        lam <- lambda*rep(1,n)
        w <- eigval + lam
        arg1 <- (w==0) & (alpha == 0); 
        arg2 <- (w==0) & (alpha != 0);
        coeff[w !=0] <- alpha[w !=0] / w[w!=0];
        coeff[arg1] = 0;
        coeff[arg2] = Inf;
        coeff[is.nan(coeff)]=0;
        s <- V%*%coeff; 
        nrms <- norm(s);
        if( (nrms > 1.2*delta) || (nrms < .8*delta)){
          key <- 5
          lambda <- -mineig
        }
      }
    }else{
      lambda <- -mineig
      key <- 4
    }
    lam <- lambda*rep(1,n)
    if(key>2){
      arg <- abs(eigval+lam) < 10*eps*max(abs(eigval),1)
      alpha[arg]<-0
    }
    w <- eigval+lam
    arg1 <- w==0 & alpha==0
    arg2 <- w==0 & alpha!=0
    coeff[w!=0]<-alpha[w!=0]/w[w!=0]
    coeff[arg1] <- 0
    coeff[arg2] <- Inf
    coeff[is.nan(coeff)] <- 0
    s <- V%*%coeff
    nrms <- norm(s)
    if(key>2 && (nrms < 0.8*delta)){
      beta <- sqrt(delta^2-nrms^2)
      s <- s+ beta*sig*V[,jmin]
    }
    if((key>2) && (nrms>1.2*delta)){
      tmp1 <- uniroot(f=seceqn,eigval=eigval,alpha=alpha,delta=delta,lower=laminit,upper=1e100,tol=tol,maxiter=itbnd)
      b <- tmp1$root
      count <- tmp1$iter
      lambda <- b
      lam <- lambda*rep(1,n)
      w <- eigval+lam
      arg1 <- w==0 & alpha==0
      arg2 <- w==0 & alpha!=0
      coeff[w!=0] <- alpha[w!=0]/w[w!=0]
      coeff[arg1] <- 0
      coeff[arg2] <- Inf
      coeff[is.nan(coeff)] <- 0
      s <- V%*%coeff
      nrms <- norm(s)
    }
  }
  val <- t(g)%*%s+t(0.5*s)%*%(H%*%s)
  return(list(s=s,val=val,posdef=posdef,count=count,lambda=lambda))
}
snls <- function()
{
  n <- length(xstart)
  xcurr <- xstart
  msgData <- NULL
  dnewt <- vector()
  iter <- 0
  numFunEvals <- 1 ## done in calling function lsqnonlin
  numGradEvals <- 1
  tol2 <- options$TolX
  tol1 <- options$TolFun
  itb <- options$MaxIter
  maxfunevals <- options$MaxFunEvals
  pcgtol <- options$TolPCG
  kmax <- options$MaxPCGIter
  ############
  ex <-0
  posdef <- 1
  npcg <- 0
  pcgit <-0
  
  delta <- 10
  nrmsx <- 1
  ratio <- 0
  
  degen <- Inf
  dv <- rep(1,n)
  x <- xstart
  oval <- Inf
  nbnds <- 1
  Z <- vector()
  if(){
    nbnds <- 0
    degen <- -1
  }
  
  xcurr <- xstart
  fvec <- full(fval) ????????????
  ## %   Evaluate F and J
  
  
  
  
  numFunEvals <- numFunEvals + findiffevals
  delbnd <- max(100*norm(xstart),1)
  size_fvec <- dim(fvec)
  mm <- size_fvec[1]
  pp <- size_fvec[2]
  if(mm < n) stop("Note enough equations for the parameters,M<N")
  if(pp==2) dnewt <- fvec[,2]
  ## %   Determine gradient of the nonlinear least squares function
  g <- ?????????????
  ## %   Evaluate F (initial point)
  val <- t(fvec[,1])%*%fvec[,1]
  ## Create the fault tolerance structure: no need in our case
  
  ## Main Loop: Generates feasible sequence x(iter), s.t. ||F(x(iter))|| is decreasing.
  while(!ex){
    ## Update:
    vres <- definev(g,x,l,u)
    v <- vres$v
    dv <- vres$dv
    
    gopt <- v*g
    optnrm <- norm(matrix(gopt),type="I") ## get the infinity norm
    r <- abs(min(u-x,x-l))
    degen <- min(r+abs(g))
    if(!nbnds) degen <--1
    
    ## Display:(ignore first)
    
    ## Test for Convergence
    diff <- abs(oval-val)
    oval <- val
    if((optnorm<tol1) && (posdef==1)){
      ex <- 1
      EXITFLAG <- 1
      if(iter==1) msgFlag <-100 else msgFlag <- EXITFLAG
    }else{
      if ((nrmsx < .9*delta) && (ratio > .25) && (diff < tol1*(1+abs(oval)))){
        ex <- 2
        EXITFLAG <- 3
        msgData <- paste('snls')
      }else{
        if((iter > 1) && (nrmsx < tol2)){
          ex <- 3
          EXITFLAG <- 2
          msgData <- paste('snls')
        }else{
          if(iter > itb){
            ex <- 4
            EXITFLAG <-0
            msgData <- paste('snls')
          }else{
            if(numFunEvals > maxfunevals){
              ex <- 4
              EXITFLAG <- 0
            }
          }
        }
      }
    }
    ## % Reset fault tolerance structure I think no need in R
    ## -------------
    ## %     Continue if ex = 0 (i.e., not done yet)
    if(ex != 0){
      ## Determine the trust region correction ******** using trdog()
      dd <- abs(v)
      D <- ??????
      sx <- rep(0,n)
      theta <- max(.95,1-optnrm)
      oposdef <- posdef
      trdogRes <- trdog(x,g,A,D,delta,dv,)
      posdef <- ifelse(is.null(trdogRes$posdef),oposdef,trdogRes$posdef)
      nrmsx <- norm(trdogRes$snod)
      npcg <- npcg + trdogRes$pcgit
      newx <- x+trdogRes$sx
      ## Newton Reflection:
      pertRes <- perturbTrustRegionReflective(newx,l,u)
      pert <- pertRes$pert
      newx <- pertRes$newx
    ## -------------
      xcurr <- newx
    ## Evaluate F and J
    
    ## update the fault tolerance structure: In R, there is no need
    ## to check whether it is NaN or Inf or  complex.
    
    ## update the trust region radius:
    
    ## Accept or reject the trial point:
    if( 1 && (newval<val)){
      x <- newx
      val <- newval
      g <- newgrad
      A <- newA
      Z <- ?
      fvec <- newfvec
      if(pp==2) dnewt <- newfvec(,2) ## Extract the newton direction
      
    }
    iter <- iter+1
    } ## end  if(ex != 0)
  }## end  while(!ex)
}

definv <- function(g,x,l,u)
{
  ## Scaling vector and derivative
#   [v,dv]= DEFINEV(g,x,l,u) returns v, distances to the
#   %   bounds corresponding to the sign of the gradient g, where
#   %   l is the vector of lower bounds, u is the vector of upper 
#   %   bounds. Vector dv is 0-1 sign vector
  
  n <- length(x); 
  v <- rep(0,n);
  dv <- rep(0,n) 
  arg1 <- (g < 0)  & (u <  Inf ); 
  arg2 <- (g >= 0) & (l > -Inf);
  arg3 <- (g < 0)  & (u == Inf); 
  arg4 <- (g >= 0) & (l == -Inf);
  v[arg1]= (x[arg1] - u[arg1]); 
  dv[arg1] = 1;
  v[arg2]  = (x[arg2] - l[arg2]); 
  dv[arg2] = 1;
  v[arg3]  = -1;
  dv[arg3] = 0;
  v[arg4]  = 1;
  dv[arg4] = 0;
  return(list(v=v,dv=dv))
}

atamult <- function(A,Y,flag=0)
{
  if(flag==0) V <- t(A)%*%(A%*%Y) else{
    if(flag < 0) V <- t(A)%*%Y else V <- A%*%Y
  }
  return(V)
}

perturbTrustRegionReflective <- function(x,l,u,del=100*.Machine$double.eps,y=NULL,sigma=NULL)
{
#   perturbs the current 
#   %   point X slightly to shake it loose from tight (less than DEL away)
#   %   bounds U and L to be strictly feasible.  also perturbs 
#  %   the reflected point Y with respect to SIGMA,
  if(min(abs(u-x))<del | min(abs(x-l))<del){
    upperi <- (u-x) < del; 
    loweri <- (x-l) < del;
    x[upperi] = x[upperi] - del;
    x[loweri] = x[loweri] + del;
    if(!is.null(y)){
      y[upperi] = y[upperi] - del*sigma[upperi];
      y[loweri] = y[loweri] + del*sigma[loweri];
    }
    pert<-TRUE
  }else pert <- FALSE
  return(list(x=x,y=y,pert=pert))
}
##' Banded preconditioner function for least-squares problems.
##' the sparse nonsingular upper triangular matrix RPCMTX such that
##' RPCMTX'*RPCMTX is a preconditioner for the
##' matrix M = DM*(A'*A)*DM + DG, where DM is a positive
##' diagonal matrix, DG is a non-negative diagonal matrix,
##' and A is sparse rectangular matrix with more rows than columns.
##' PPVEC is the associated permutation (row) vector and UPPERBW
##' specifies the upperbandwidth of RPCMTX.
aprecon <- function(A=NULL,upperbandw=NULL,DM,DG,...)
{
  tmp <- dim(A)
  n <- length(DM) ??????
  RPCMTX <- Matrix(diag(n)) 
  ppvec <- 1:n
  ## % In case "A" isn't really A, but something else to use with JacobMult function.
  if((!(is.numeric(A)))||(is.null(A))||(tmp[1]<tmp[2])||(tmp[2]!=n)){
  
    d1 <- diag(DM)
    d2 <- diag(DG)
    dd <- sqrt(d1*d1+abs(d2))
    ## RPCMTX <- matrix(0,n,n) ????????
    ## diag(RPCMTX) <- dd
    RPCMTX <- Matrix(diag(dd))
    
  }else{
    epsi <- sqrt(.Machine$double.eps)
    if(is.null(upperbandw)) upperbandw <-0
    ## Form Matrix M
    TM <- A%*%DM
    if(upperbandw==0){
      ## M <- t(TM)%*%TM+DG
      M <- crossprod(TM)+DG
      dnrms <- sqrt(sum(M*M)) ??????
      d <- max(sqrt(dnrms),epsi)
      RPCMTX <- Matrix(diag(d))
      # ppvec <- 1:n ## no need to do it
      ## But what if RPCMTX is singular?
    }else{
      if(upperbandw >= (n-1)) {## Attempt sparse QR
        dgg <- sqrt(diag(DG))
        DDG <- diag(dgg)
        TM <- rbind(TM,DDG)???
        ##p <- colamd(TM)????
        ##RPCMTX <- qr(TM[,p])
        ##RPCMTX <- RPCMTX[1:n,1:n]
        ##ppvec <- p
        tmp <- qr(TM,pivot=TRUE)
        RPCMTX <- tmp$qr
        ppvec <- tmp$pivot
        ## %    Modify for singularity?
        mdiag <- min(abs(diag(RPCMTX)))
        lambda <- 1
        while(mdiag <- epsi){
          TM <- rbind(A%*%M,DDG+lambda*diag(n))##???
          lambda <- 4*lambda
#           p <- colamd(TM)
#           RPCMTX <- qr(TM[,p])
#           RPCMTX <- RPCMTX[1:n,1:n]
#           ppvec <- p
          tmp <- qr(TM,pivot=TRUE)
          RPCMTX <- tmp$qr
          ppvec <- tmp$pivot            
          mdiag <- min(abs(diag(RPCMTX)))
        }
      }else{
        if(upperbandw>0 && (upperbandw < (n-1))){
          ## % Band approximation.
          M <- crossprod(TM)+DG
          p <- 1:n
          M <- tril(triu(M(p,p),-upperbandw),upperbandw)??????
          RPCMTX <- Matrix(0,n,n)
          tmp <- try(chol(M),silent=TRUE)
          lambda <- 1
          while(attr(tmp,"class")=="try-error"){
            M <- M+lambda*Matrix(diag(n))
            tmp <- try(chol(M),silent=TRUE)
            lambda <- 4*lambda
          }
          RPCMTX <- tmp
          ppvec <- p
        }else{
          stop("Invalid upperbandw")
        }
      }
    }

    
  }
  
  
  return(list(RPCMTX=RPCMTX,ppvec=ppvec))
}
pcgr <- function(DM,DG,g,kmax,tol,mtxmpy=atamult,
                 H,R,pR,callerflag,
                 pcoptions,...)
{## Preconditioned conjugate gradients
  n <- nrow(DG)
  r <- -g
  p <- rep(0,n)
  
  ## precondition
  z <- preproj(r,R,pR)
  ynrm <- nrom(Z)
  stoptol <- tol*znrm
  inner2 <- 0
  inner1 <- t(r)%*%z
  posdef <- 1
  
  kmax <- max(kmax,1)
  
  for(k in 1:kmax){
    if(k==1){
      d <- z
    }else{
      beta <- inner1/inner2
      d <- z+beta*d
    }
    ww <- DM %*% d
    if(callerflag=="jacobprecon") w <- mtxmpy(H,ww,...)
    ww <- DM%*%w+DG%*%d
    denom <- t(d)%*%ww
    if(denom <=0){
      if(norm(d)==0){
        p <- d
      }else p <- d/norm(d)
      posdef <-0
      break
    }else{
      alpha <- inner1/denom
      p <- p+alpha*d
      r <- r-alpha*ww
    }
    z <- preproj(r,R,pR)
    if(norm(z)<=stoptol) break
    inner2 <- inner1
    inner1 <- t(r)%*%z
  }
  
  if(pcoptions==Inf) k <- 0
  return(list(p=p,posdef=posdef,k=k))
  
  
}

preproj <- function(r,RPCMTX=NULL,ppvec=NULL)
{
  ## %PREPROJ Apply preconditioner
  n <- length(r)
  if(is.null(ppvec)) ppvec <- 1:n
  if(is.null(RPCMTX)) RPCMTX <- Matrix(diag(n))
  ## wbar <- t(RPCMTX)\r[ppvec]
  wbar <- solve(t(RPCMTX),r[ppvec])
  ## w[ppvec,1] <- RPCMTX\wbar
  w <- rep(NA,n)
  w[ppvec] <- solve(RPCMTX,wbar)
  return(w)
}
trdog <- function(x,g,H,D,delta,dv,
                  mtxmpy=atamult,pcmtx=aprecon,pcoptions=Inf,
                  tol=0.1,kmax=2,theta=0.95,
                  l=-Inf,u=Inf,Z=NULL,
                  dnewt=NULL,preconflag="jacobprecon",...)
{
  ## a clone of the trdog.m function in Matlab
#   %TRDOG Reflected (2-D) trust region trial step (box constraints)
#   %
#   % [s,snod,qpval,posdef,pcgit,Z] = TRDOG(x,g,H,D,delta,dv,...
#      %                 mtxmpy,pcmtx,pcoptions,tol,theta,l,u,Z,dnewt,preconflag);
#   %
#   %   Determine the trial step `s', an approx. trust region solution.
# %   `s' is chosen as the best of 3 steps: the scaled gradient
# %   (truncated to  maintain strict feasibility),
#   %   a 2-D trust region solution (truncated to remain strictly feas.),
#   %   and the reflection of the 2-D trust region solution,
#   %   (truncated to remain strictly feasible).
#   %
#   %   The 2-D subspace (defining the trust region problem) is defined 
#   %   by the scaled gradient direction and a CG process (returning
#    %   either an approximate Newton step of a direction of negative curvature.
#      Driver functions are: SNLS, SFMINBX
#        SNLS actually calls TRDOG with the Jacobian matrix (and a special 
#      Jacobian-matrix multiply function in MTXMPY).
#                                                          
#     Copyright 1990-2011 The MathWorks, Inc.
#      $Revision: 1.1.6.3 $  $Date: 2011/05/09 01:16:31 $
#      
  
  
#       Initialization
  ##require(Matrix)
  n = length(g);  # g is the gradient
  pcgit = 0; 
  grad = D%*%g;    ##
  DM = D;  
  #%DG = sparse(1:n,1:n,full(abs(g).*dv)); ## Matlab implementation
  # DG <- diag(abs(g)*dv)  ## R implementation without sparse matrix definition
  DG <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  diag(DG) <- abs(g)*dv
  posdef = 1; 
  pcgit = 0; 
  tol2 = sqrt(eps);
  v1 = dnewt; 
  qpval1 = Inf; 
  qpval2 = Inf; 
  qpval3 = Inf;
  
  ## % DETERMINE A 2-DIMENSIONAL SUBSPACE
  if(is.null(Z)){
    if(is.null(v1)){
      if(preconflag=="jacobprecon") {
        tmp <- pcmtx(H,pcoptions,DM,DG,...)
        R <- tmp$R
        permR <- tmp$permR
      }else{
        stop("Not implemented yet!!")
      }
      if(tol <= 0) tol <- 0.1
      tmp <- pcgr(DM,DG,grad,kmax,tol,mtxmpy,H,R,permR,preconflag,pcoptions,...)
      v1 <- tmp$v1
      posdef <- tmp$posdef
      pcgit$tmp$pcgit
    }
    
    nrmv1 <- norm(v1)
    if(nrmv1>0) v1 <- v1/nrmv1
    Z <- matrix(NA,n,2)
    Z[,1] <- v1
    if(n>1){
      if(posdef<1){
        v2 <- D%*%sign(grad)
        nrmv2 <- norm(v2)
        if(nrmv2>0) v2 <- v2/nrmv2
        v2 <- v2 - v1%*%(t(v1)%*%v2)
        nrmv2 <- norm(v2)
        if(nrmv2>tol2) {
          v2 <- v2/nrmv2
          Z[,2] <- v2
        }
      }else{
        nrmgrad <- norm(grad)
        if(nrmgrad>0) v2 <- grad/nrmgrad else v2 <- grad
        v2 <- v2 - v1%*%(t(v1)%*%v2)
        nrmv2 <- norm(v2)
        if(nrmv2>tol2) {
          v2 <- v2/nrmv2
          Z[,2] <- v2
        }
      }
    }
  }## end  if(is.null(Z))
  W <- DM%*%Z
  if(preconflag=="jacobprecon") {
    WW <- mtxmpy(H,W,0)
  }else stop("not implemented yet!!")
  WW <- DM%*%WW
  MM <- t(Z)%*%W+t(Z)%*%DG%*%Z
  rhs <- t(Z)%*%grad
  tmp <- trust.exact(rhs,MM,delta)
  st <- tmp$s
  ss <- Z%*%st
  s <- abs(diag(D))*ss
  ssave <- s
  sssave <- ss
  stsave <- st
  
  ## check direction for NaNs
  if(any(is.nan(s))) stop("erros in s, NaNs")
  ## Truncate the TR solution
  arg <- abs(s)>0
  if(!any(arg)){
    alpha <- 1
    mmdis <- 1
    ipt <- vector()
  }else{
    dis <- pmax((u(arg)-x(arg))/s(arg), (l(arg)-x(arg))/s(arg))
    ipt <- which.min(dis)
    mmdis <- dis[ipt]
    mdis <- theta * mmdis
    alpha <- min(1,mdis)
  }
  s <- alpha*s
  st <- alpha*st
  ss <- alpha*ss
  qpval1 <- t(rhs)*st+t(.5*st)%*%MM%*%st
  if(n>1){
    ##%   Evaluate along the reflected direction?
    qpval3 <- Inf
    ssssave <- mmdis*sssave
    if(norm(ssssave)<.9*delta){
      r <- mmdis*ssave
      nx <- x+r
      stsave <- mmdis*stsave
      qpval0 <- t(rhs)%*%stsave+t(0.5*stsave)%*%MM%*%stsave
      if(preconflag=="jacobprecon") ng <- mtxmpy(H,r,0)
      ng <- ng+g
      ngrad <- D%*%ng
      ngrad <- ngrad +DG%*%ssssave
      ## %      nss is the reflected direction
      nss <- sssave
      nss[ipt]<- -nss[ipt]
      ZZ[,1]<-nss/(norm(nss))
      W<-DM%*%ZZ
      if(preconflag=="jacobprecon") WW <- mtxmpy(H,W,0)
      W <- DM%*%WW
      MM <- t(ZZ)%*%W+t(ZZ)%*%DG%*%ZZ
      nrhs <- t(ZZ)%*%ngrad
      tmp1 <- quad1d(nss,ssssave,delta)
      nss <- tmp1$nss
      tau <- tmp1$tau
      nst <- tau/norm(nss)
      ns <- abs(diag(D))*nss
      ## check direction for NaNs
      if(any(is.nan(s))) stop("erros in s, NaNs")
      ## % Truncate the reflected direction?
      arg <- (abs(ns) > 0)
      if(!(any(arg))) alpha <- 1 else{
        dis <- pmax((u(arg)-nx(arg))/ns(arg), (l(arg)-nx(arg))/ns(arg));
        mdis <- min(dis);  
        mdis <- theta*mdis;
        alpha <- min(1,mdis)
      }
      ns <- alpha*ns
      nst <- alpha*nst
      nss <- alpha*nss
      qpval3 <- qpval0+t(nrhs)%*%nst+t(0.5*nst)%*%MM%*%nst
    }
    gnorm <- norm(grad)
    ZZ[,1]<- grad/(gnorm+(gnorm==0))## % Protect against norm of 0
    W <- DM%*%ZZ
    if(preconflag=="jacobprecon") WW <- mtxmpy(H,W,0)
    W <- DM%*%WW
    MM <- t(ZZ)%*%W+t(ZZ)%*%DG%*%ZZ
    rhs <- t(ZZ)%*%grad
    tmp <- trust.exact(rhs,MM,delta)
    st <- tmp$s
    ssg <- ZZ%*%st
    sg <- abs(diag(D))*ssg
    if(any(is.nan(sg))) stop("erros in sg, NaNs")
    arg <- (abs(sg) > 0)
    if(!(any(arg))) alpha <- 1 else{
      dis = pmax((u(arg)-x(arg))/sg(arg), (l(arg)-x(arg))/sg(arg))
      mdis <- min(dis)
      mdis <- theta*mdis;
      alpha <- min(1,mdis)
    }
    sg <- alpha*sg
    st <- alpha*st
    ssg <- alpha*ssg
    qpval2 <- t(rhs)%*%st+t(0.5*st)%*%MM%*%st
  }
  
  ## Choose the best of s, sg, ns.
  if(qpval2 <0 min(qpval1,qpval3)){
    qpval <- qpval2
    s <- sg
    snod <- ssg
  }else{
    if(qpval1 <= min(qpval2,qpval3)){
      qpval <- qpval1
      snod <- ss
    }else{
      qpval <- qpval3
      s <- ns+r
      snod <- nss+ssssave
    }
          
  }
  return(list(s=s,snod=snod,qpval=qpval,posdef=posdef,pcgit=pcgit,Z=Z))
                                                         
}

## utilities
gradfun_ls <- function(modelfun,par,obs,...)
{
  if(is.null(names(par))) names(par) <- pnames
  np <- length(par)
  resfun <- function(par,obs, ...){
    modelfun(par,...)-obs
  }
  res <- resfun(par,obs=obs)
  Jmat <- jacobian(modelfun,x=par)
  g <- t(Jmat)%*%res
  return(g)
}

##' @examples
##' \dontrun{
##' farb<- function(x=c(3,1),t=1:10,...){
##'   x1 <- x[1]
##'   x2 <- x[2]
##'   x[1]+x[2]*t
##'   }
##'   y <- farb(x=c(3,1),t=1:10)+rnorm(10,0,2)
##'   t <- 1:10
##'   x <- c(3,1)
##'   objfun <- objfun_ls(farb,c(3,1),obs=y)
##'   }
##'   
objfun_ls <- function(modelfun,par,obs,Q=FALSE,pnames=NULL,...)
{
  np <- length(par)
  if(is.null(names(par))) names(par) <- pnames
  resfun <- function(par,obs, ...){
    modelfun(par,...)-obs
  }
  objfun <- function(par,obs,...)
  {
    0.5*sum((modelfun(par,...)-obs)^2)
  }

  res <- resfun(par,obs=obs)
  ## calculate the jacobian of model function(returning a vector)
  id <- which(is.na(res))
  
  if(Q==TRUE)
  {
    D <- genD(modelfun,x=par)$D
    Jmat <- D[,1:np]
    g <- t(Jmat[-id,])%*%res[-id]
    Hmat1 <- D[,(np+1):(np+np*(np+1)/2)]
    LT <- res[-id]%*%Hmat1[-id,]
    M <- matrix(NA,np,np)
    k <- 0
    for(i in 1:np)
      for(j in 1:i)
      {
        k <- k+1
        M[i,j] <- LT[k]
      }
    ## M is symmetric
    ##browser()
    M[upper.tri(M)] <- t(M)[upper.tri(M)]
    B <- t(Jmat[-id,])%*%Jmat[-id,]+M
  }else {
    Jmat <- jacobian(modelfun,x=par)
    g <- t(Jmat[-id,])%*%res[-id]
    
    B <- t(Jmat[-id,])%*%Jmat[-id,]
  }
  ## using numerical approximation : grad(objfun,x=c(3,1),obs=y)
  ## using numerical approximation : Hmat <- hessian(objfun,x=par,obs=obs)
  
  list(value=0.5*sum(res^2,na.rm=TRUE),gradient=g,hessian=B)  
}

