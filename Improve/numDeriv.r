############################################################################

#    functions for gradient calculation

############################################################################

grad <- function (func, x, method="Richardson", method.args=list(), ...)
  UseMethod("grad")

grad.default <- function(func, x, method="Richardson",
                         method.args=list(), ...){
  # modified by Paul Gilbert from code by Xingqiao Liu.
  # case 1/ scalar arg, scalar result (case 2/ or 3/ code should work)
  # case 2/ vector arg, scalar result (same as special case jacobian)
  # case 3/ vector arg, vector result (of same length, really 1/ applied multiple times))
  f <- func(x, ...)
  n <- length(x)   #number of variables in argument
  case1or3 <- n == length(f)
  if((1 != length(f)) & !case1or3)
    stop("grad assumes a scalar valued function.")
  if(method=="simple"){
    #  very simple numerical approximation
    args <- list(eps=1e-4) # default
    args[names(method.args)] <- method.args
    eps <- args$eps
    if(case1or3) return((func(x+eps, ...)-f)/eps) 
    # now case 2
    df <- rep(NA,n)
    for (i in 1:n) {
      dx <- x
      dx[i] <- dx[i] +eps 
      df[i] <- (func(dx, ...)-f)/eps
    }
    return(df)
  } else
    if(method=="complex"){ # Complex step gradient
      eps <- .Machine$double.eps
      v <- try(func(x + eps * 1i, ...))
      if(inherits(v, "try-error")) 
        stop("function does not accept complex argument as required by method 'complex'.")
      if(!is.complex(v)) 
        stop("function does not return a complex value as required by method 'complex'.")
      
      if(case1or3) return(Im(v)/eps) 
      # now case 2
      h0 <- rep(0, n)
      g  <- rep(NA, n)
      for (i in 1:n) {
        h0[i] <- eps * 1i
        g[i] <- Im(func(x+h0, ...))/eps 
        h0[i]  <- 0
      }
      return(g)
    } else
      if(method=="Richardson"){
        args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE) # default
        args[names(method.args)] <- method.args
        eps <- args$eps
        d <- args$d
        r <- args$r
        v <- args$v
        show.details <- args$show.details
        a <- matrix(NA, r, n) 
        #b <- matrix(NA, (r - 1), n)
        
        #  first order derivatives are stored in the matrix a[k,i], 
        #  where the indexing variables k for rows(1 to r), i for columns (1 to n),
        #  r is the number of iterations, and n is the number of variables.
        
        h <- abs(d*x)+eps*(abs(x) < args$zero.tol)
        for(k in 1:r)  { # successively reduce h		    
          if(case1or3)  a[k,] <- (func(x + h, ...) -  func(x - h, ...))/(2*h)
          else for(i in 1:n)  {
            if((k != 1) && (abs(a[(k-1),i]) < 1e-20)) a[k,i] <- 0 #some func are unstable near zero
            else  a[k,i] <- (func(x + h*(i==seq(n)), ...) - 
              func(x - h*(i==seq(n)), ...))/(2*h[i])
          }
          if (any(is.na(a[k,]))) stop("function returns NA at ", h," distance from x.")
          h <- h/v     # Reduced h by 1/v.
        }	
        if(show.details)  {
          cat("\n","first order approximations", "\n")		
          print(a, 12)
        }
        
        #------------------------------------------------------------------------
        # 1 Applying Richardson Extrapolation to improve the accuracy of 
        #   the first and second order derivatives. The algorithm as follows:
        #
        #   --  For each column of the derivative matrix a,
        #	  say, A1, A2, ..., Ar, by Richardson Extrapolation, to calculate a
        #	  new sequence of approximations B1, B2, ..., Br used the formula
        #
        #	     B(i) =( A(i+1)*4^m - A(i) ) / (4^m - 1) ,  i=1,2,...,r-m
        #
        #		N.B. This formula assumes v=2.
        #
        #   -- Initially m is taken as 1  and then the process is repeated 
        #	 restarting with the latest improved values and increasing the 
        #	 value of m by one each until m equals r-1
        #
        # 2 Display the improved derivatives for each
        #   m from 1 to r-1 if the argument show.details=T.
        #
        # 3 Return the final improved  derivative vector.
        #-------------------------------------------------------------------------
        
        for(m in 1:(r - 1)) {	  
          a <- (a[2:(r+1-m),,drop=FALSE]*(4^m)-a[1:(r-m),,drop=FALSE])/(4^m-1)
          if(show.details & m!=(r-1) )  {
            cat("\n","Richarson improvement group No. ", m, "\n") 	  
            print(a[1:(r-m),,drop=FALSE], 12)
          }
        }
        return(c(a))
      } else stop("indicated method ", method, "not supported.")
}


jacobian <- function (func, x, method="Richardson",
                      method.args=list(), ...) UseMethod("jacobian")

jacobian.default <- function(func, x, method="Richardson",
                             method.args=list(), ...){
  f <- func(x, ...)
  n <- length(x)	 #number of variables.
  if(method=="simple"){
    #  very simple numerical approximation
    args <- list(eps=1e-4) # default
    args[names(method.args)] <- method.args
    eps <- args$eps
    df <-matrix(NA, length(f), n)
    for (i in 1:n) {
      dx <- x
      dx[i] <- dx[i] +eps 
      df[,i] <- (func(dx, ...)-f)/eps
    }
    return(df)
  } else
    if(method=="complex"){ # Complex step gradient
      # Complex step Jacobian
      eps <- .Machine$double.eps
      h0  <-  rep(0, n)
      h0[1] <- eps * 1i
      v <- try(func(x+h0, ...))
      if(inherits(v, "try-error")) 
        stop("function does not accept complex argument as required by method 'complex'.")
      if(!is.complex(v)) 
        stop("function does not return a complex value as required by method 'complex'.")
      
      h0[1]  <- 0
      jac <- matrix(NA, length(v), n)
      jac[, 1] <- Im(v)/eps
      if (n == 1) return(jac)
      for (i in 2:n) {
        h0[i] <- eps * 1i
        jac[, i] <- Im(func(x+h0, ...))/eps 
        h0[i]  <- 0
      }
      return(jac)
    } else
      if(method=="Richardson"){
        args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), 
                     r=4, v=2, show.details=FALSE) # default
        args[names(method.args)] <- method.args
        eps <- args$eps
        d <- args$d
        r <- args$r
        v <- args$v 	  
        a <- array(NA, c(length(f),r, n) )
        
        h <- abs(d*x)+eps*(abs(x) < args$zero.tol)
        for(k in 1:r)  { # successively reduce h		    
          for(i in 1:n)  {
            a[,k,i] <- (func(x + h*(i==seq(n)), ...) -  
              func(x - h*(i==seq(n)), ...))/(2*h[i])
            #if((k != 1)) a[,(abs(a[,(k-1),i]) < 1e-20)] <- 0 #some func are unstable near zero
          }
          h <- h/v     # Reduced h by 1/v.
        }	
        for(m in 1:(r - 1)) {	  
          a <- (a[,2:(r+1-m),,drop=FALSE]*(4^m)-a[,1:(r-m),,drop=FALSE])/(4^m-1)
        }
        # drop second dim of a, which is now 1 (but not other dim's even if they are 1
        return(array(a, dim(a)[c(1,3)]))  
      } else stop("indicated method ", method, "not supported.")
}


hessian <- function (func, x, method="Richardson",
                     method.args=list(), ...) UseMethod("hessian")

hessian.default <- function(func, x, method="Richardson",
                            method.args=list(), ...){  
  
  if(1!=length(func(x, ...)))
    stop("Richardson method for hessian assumes a scalar valued function.")
  
  if(method=="complex"){ # Complex step hessian
    args <- list(eps=1e-4, d=0.1, 
                 zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)
    args[names(method.args)] <- method.args
    # the CSD part of this uses eps=.Machine$double.eps
    # but the jacobian is Richardson and uses method.args
    return(jacobian(func=function(fn, x, ...){grad(func=fn, x=x, 
                                                   method="complex", method.args=list(eps=.Machine$double.eps), ...)}, 
                    x=x, fn=func, method.args=args, ...))
  } else 
    if(method != "Richardson")  stop("method not implemented.")
  args <- list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), 
               r=4, v=2, show.details=FALSE) # default
  args[names(method.args)] <- method.args
  D <- genD(func, x, method=method, method.args=args, ...)$D
  if(1!=nrow(D)) stop("BUG! should not get here.")
  H <- diag(NA,length(x))
  u <- length(x)
  for(i in 1:length(x))
  {for(j in 1:i) 
  {u <- u + 1
   H[i,j] <- D[,u]
  }  }
  H <- H +t(H)
  diag(H) <- diag(H)/2
  H
}


#######################################################################

#               Bates & Watts   D matrix calculation

#######################################################################

genD <- function(func, x, method="Richardson",
                 method.args=list(), ...)UseMethod("genD")

genD.default <- function(func, x, method="Richardson",
                         method.args=list(), ...){
  #   additional cleanup by Paul Gilbert (March, 2006)
  #   modified substantially by Paul Gilbert (May, 1992)
  #    from original code by Xingqiao Liu,   May, 1991.
  
  #  This function is not optimized for S speed, but is organized in 
  # the same way it could be (was) implemented in C, to facilitate checking.
  
  #  v  reduction factor for Richardson iterations. This could
  #	 be a parameter but the way the formula is coded it is assumed to be 2.
  
  if(method != "Richardson")  stop("method not implemented.")
  args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7),
               r=4, v=2) # default
  args[names(method.args)] <- method.args
  eps <- args$eps
  d <- args$d
  r <- args$r
  v <- args$v
  if (v!=2) stop("The current code assumes v is 2 (the default).")	 
  #func.args <- list(...)
  
  #f0 <- do.call("func",append(list(x), func.args))
  f0 <- func(x, ...)
  #  f0 is the value of the function at x.
  
  p <- length(x)  #  number of parameters (theta)
  h0 <- abs(d*x)+eps*(abs(x) < args$zero.tol)
  D <- matrix(0, length(f0),(p*(p + 3))/2)
  #length(f0) is the dim of the sample space
  #(p*(p + 3))/2 is the number of columns of matrix D.( first
  #   der. & lower triangle of Hessian)
  Daprox <- matrix(0, length(f0),r) 
  Hdiag  <-  matrix(0,length(f0),p)
  Haprox <-  matrix(0,length(f0),r)
  for(i in 1:p)    # each parameter  - first deriv. & hessian diagonal
  {h <-h0
   for(k in 1:r)  # successively reduce h 
   {f1 <- func(x+(i==(1:p))*h, ...)
    f2 <- func(x-(i==(1:p))*h, ...) 
    #f1 <- do.call("func",append(list(x+(i==(1:p))*h), func.args))
    #f2 <- do.call("func",append(list(x-(i==(1:p))*h), func.args))
    Daprox[,k] <- (f1 - f2)  / (2*h[i])    # F'(i) 
    Haprox[,k] <- (f1-2*f0+f2)/ h[i]^2     # F''(i,i) hessian diagonal
    h <- h/v     # Reduced h by 1/v.
    NULL
   }
   for(m in 1:(r - 1))
     for ( k in 1:(r-m))
     {Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1)
      Haprox[,k]<-(Haprox[,k+1]*(4^m)-Haprox[,k])/(4^m-1)
      NULL
     }
   D[,i] <- Daprox[,1]
   Hdiag[,i] <- Haprox[,1]
   NULL
  }	  
  u <- p
  
  for(i in 1:p)   # 2nd derivative  - do lower half of hessian only
  {for(j in 1:i) 
  {u <- u + 1
   if (i==j) { D[,u] <- Hdiag[,i]; NULL}
   else 
   {h <-h0
    for(k in 1:r)  # successively reduce h 
    {f1 <- func(x+(i==(1:p))*h + (j==(1:p))*h, ...)
     f2 <- func(x-(i==(1:p))*h - (j==(1:p))*h, ...)
     #f1 <- do.call("func", append(
     #	  list(x+(i==(1:p))*h + (j==(1:p))*h), func.args))
     #f2 <- do.call("func",append(
     #	  list(x-(i==(1:p))*h - (j==(1:p))*h), func.args))  
     Daprox[,k]<- (f1 - 2*f0 + f2 -
       Hdiag[,i]*h[i]^2 - 
       Hdiag[,j]*h[j]^2)/(2*h[i]*h[j])  # F''(i,j)  
     h <- h/v	# Reduced h by 1/v.
    }
    for(m in 1:(r - 1))
      for ( k in 1:(r-m))
      {Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1); NULL}
    D[,u] <- Daprox[,1]
    NULL
   }
  }  
  }
  D <- list(D=D, p=length(x), f0=f0, func=func, x=x, d=d,
            method=method, method.args=args)# Darray constructor (genD.default)
  class(D) <- "Darray"
  invisible(D)
}

